def read_hires(path,library_id):
	import scanpy as sc
	from h5py import File
	import cv2
	import pandas as pd
	import json
	adata = sc.read_10x_h5(path +"filtered_feature_bc_matrix.h5")
	adata.uns["spatial"] = dict()
	with File(path +"filtered_feature_bc_matrix.h5", mode="r") as f:
		attrs = dict(f.attrs)
	adata.uns["spatial"][library_id] = dict()
	adata.uns["spatial"][library_id]["images"] = dict()
	adata.uns["spatial"][library_id]["images"]["hires"] = cv2.imread(path +"/spatial/tissue_hires_image.png")
	adata.uns["spatial"][library_id]['use_quality']="hires"
	with open(path +"spatial/scalefactors_json.json",'r') as jsonfile:
		adata.uns["spatial"][library_id]["scalefactors"] = json.load(jsonfile)
	positions = pd.read_csv(path +"spatial/tissue_positions.csv",header=None,index_col=0)
	positions.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
	adata.obs = adata.obs.join(positions, how="left")

	adata.obs['imagecol']=adata.obs['pxl_col_in_fullres'].astype(float).astype(int)*adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
	adata.obs['imagerow']=adata.obs['pxl_row_in_fullres'].astype(float).astype(int)*adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
	adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()

	return adata


def sTcellinteraction(sTRNAFold, SCTFold, OutputFold):
    import stlearn as st
    import numpy as np
    import pandas as pd
    import scanpy as sc
   
    sp_adata = read_hires(sTRNAFold, library_id="spatial")
    sp_adata.var_names_make_unique()
    SCTData = sc.read_10x_mtx(SCTFold, cache=False)
    genes = np.intersect1d(SCTData.var_names.to_list(), sp_adata.var_names.to_list())
    metadata = pd.read_csv(SCTFold + "/metadata.txt", sep="\t")
    sp_adata = sp_adata[metadata.barcode.to_list(), genes]
    SCTData = SCTData[metadata.barcode.to_list(), genes]
    sp_adata.X = SCTData.X
    st.pp.filter_genes(sp_adata, min_cells=3)

    lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
    st.tl.cci.run(sp_adata, lrs, min_spots=20, distance=None, n_pairs=2000, n_cpus=None)

    if 'lr_summary' in sp_adata.uns.keys():
        st.tl.cci.adj_pvals(sp_adata, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')
        lr_interaction_score = pd.DataFrame(sp_adata.obsm["lr_scores"], index=sp_adata.obs.index,
                                            columns=sp_adata.uns["lr_summary"].index)
        lr_interaction_pval = pd.DataFrame(sp_adata.obsm["p_adjs"], index=sp_adata.obs.index,
                                           columns=sp_adata.uns["lr_summary"].index)
        lr_interaction_score.to_csv(OutputFold + "/lr_interaction_score.txt")
        lr_interaction_pval.to_csv(OutputFold + "/lr_interaction_pval.txt")
        sp_adata.uns['lr_summary'].to_csv(OutputFold + "/lr_summary.txt")
        EpiSpots = metadata.iloc[np.where(metadata.loc[:, "cnvcluster"] == "Normoxia")[0], :].loc[:, "barcode"].to_list()
        MesSpots = metadata.iloc[np.where(metadata.loc[:, "cnvcluster"] == "Hypoxia")[0], :].loc[:, "barcode"].to_list()
        LREpiMes = pd.concat(
            [lr_interaction_score.loc[EpiSpots, :].mean(0), lr_interaction_score.loc[MesSpots, :].mean(0)], axis=1)
        LREpiMes.columns = ["Normoxia", "Hypoxia"]
        LREpiMes.to_csv(OutputFold + "/lr_interaction_score_compare.txt")
    else:
        lr_interaction_score = pd.DataFrame(index=sp_adata.obs.index)
        lr_interaction_score.to_csv(OutputFold + "/lr_interaction_score.txt")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-FunctionFlag', dest='FunctionFlag', default="sTcellinteraction", type=str, help='Function aim')
    parser.add_argument('-sTcellinteraction_sTRNAFold', dest='sTcellinteraction_sTRNAFold', default="./NSCLC_STData/spaceRangerOutput/OTAR_LNGsp9476039/outs/", type=str,
                        help='Raw spatial folder')
    parser.add_argument('-sTcellinteraction_SCTFold', dest='sTcellinteraction_SCTFold', default="./NSCLC_STData/sTinteraction/preData/Interaction(OTAR_LNGsp9476039)", type=str,
                        help='Spatial data analyzed by SCT')
    parser.add_argument('-sTcellinteraction_OutputFold', dest='sTcellinteraction_OutputFold', default="./NSCLC_STData/sTinteraction/stlearnOutputCellChat/OTAR_LNGsp9476039", type=str,
                        help='output folder')
    args = parser.parse_args()

    if args.FunctionFlag == "sTcellinteraction":
        sTcellinteraction(args.sTcellinteraction_sTRNAFold, args.sTcellinteraction_SCTFold,
                          args.sTcellinteraction_OutputFold)


if __name__ == '__main__':
    main()
