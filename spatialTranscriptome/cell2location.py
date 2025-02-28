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
    adata.uns["spatial"][library_id]["images"]["hires"] = cv2.imread(
    path + "/spatial/tissue_hires_image.png")
    adata.uns["spatial"][library_id]['use_quality'] = "hires"
    with open(path + "spatial/scalefactors_json.json", 'r') as jsonfile:
        adata.uns["spatial"][library_id]["scalefactors"] = json.load(jsonfile)
    positions = pd.read_csv(path +"spatial/tissue_positions.csv",header=None,index_col=0)
    positions.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
    adata.obs = adata.obs.join(positions, how="left")
    adata.obs['imagecol'] = adata.obs['pxl_col_in_fullres'].astype(float).astype(int)*adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    adata.obs['imagerow'] = adata.obs['pxl_row_in_fullres'].astype(float).astype(int)*adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
    return adata

def sTcell2location(scRNAFold,sTRNAFold,OutputFile):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import cell2location
    import random
    from cell2location.utils.filtering import filter_genes
    from cell2location.models import RegressionModel
    from cell2location.plt import plot_spatial
    from cell2location import run_colocation

    sc_adata = sc.read_mtx(scRNAFold + "/matrix.mtx")
    sc_adata = sc_adata.T
    sc_adata.obs["cell_id"] = pd.read_csv(scRNAFold + "/barcodes.tsv", header=None)[0].tolist()
    sc_adata.obs.index = sc_adata.obs["cell_id"]
    sc_adata.var["gene_name"] = pd.read_csv(scRNAFold + "/genes.tsv", sep="\t", header=None)[1].tolist()
    sc_adata.var.index = sc_adata.var["gene_name"]

    metadata = pd.read_csv(scRNAFold + "/MetaData.txt", sep="\t",index_col=0)
    sc_adata = sc_adata[metadata.index, :]
    metadata['SampleID'] = metadata.index
    sc_adata.obs['Sample'] = metadata['SampleID']
    sc_adata.obs['batch'] = metadata['SampleID']
    sc_adata.obs['celltype'] = metadata['ident']
    sc_adata.obs['cluster'] = sc_adata.obs['celltype'].astype("category")
    sc_adata.var['symbol'] = sc_adata.var['gene_name']
    sc_adata.var_names_make_unique()

    # filter: remove MT genes, genes expressed in less than 5 cells,
    sc_adata.var['MT_gene'] = [gene.startswith('MT-') for gene in sc_adata.var['symbol']]
    sc_adata.obsm['MT'] = sc_adata[:, sc_adata.var['MT_gene'].values].X.toarray()
    sc_adata = sc_adata[:, ~sc_adata.var['MT_gene'].values]
    sc.pp.filter_genes(sc_adata, min_cells=5)
    sc_adata.var['n_cells'] = (sc_adata.X.toarray() > 0).sum(0)
    sc_adata.var['nonz_mean'] = sc_adata.X.toarray().sum(0) / sc_adata.var['n_cells']
    nonz_mean_cutoff = np.log10(2)
    sc_adata = sc_adata[:, np.array(np.log10(sc_adata.var['nonz_mean']) > nonz_mean_cutoff)]
    sc_adata = sc_adata.copy()

    cell2location.models.RegressionModel.setup_anndata(adata=sc_adata,batch_key='batch',labels_key='cluster')
    mod = RegressionModel(sc_adata)
    mod.train(max_epochs=500,train_size=1,lr=0.002,accelerator="cpu")  # epoch: more than 500; only accept count matrix
    sc_adata = mod.export_posterior(sc_adata,sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
    inf_aver = sc_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in sc_adata.uns['mod']['factor_names']]].copy()
    inf_aver.columns = sc_adata.uns['mod']['factor_names']

    sp_adata = read_hires(sTRNAFold, library_id="spatial")
    sp_adata.var['symbol'] = sp_adata.var.index
    sp_adata.var_names_make_unique()

    sp_adata.var['MT_gene'] = [gene.startswith('MT-') for gene in sp_adata.var['symbol']]
    sp_adata.obsm['MT'] = sp_adata[:, sp_adata.var['MT_gene'].values].X.toarray()
    sp_adata = sp_adata[:, ~sp_adata.var['MT_gene'].values]

    intersect = np.intersect1d(sp_adata.var_names, inf_aver.index)
    sp_adata = sp_adata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=sp_adata)
    mod = cell2location.models.Cell2location(sp_adata,cell_state_df=inf_aver,N_cells_per_location=30,detection_alpha=200)
    mod.train(max_epochs=3000,batch_size=None,train_size=1,accelerator="cpu")  # epoch more than 3000
    sp_adata = mod.export_posterior(sp_adata,sample_kwargs={'num_samples': 1000,'batch_size': mod.adata.n_obs,})
    sp_adata.var.drop(columns='gene_ids', inplace=True)
    sp_adata.obs[sp_adata.uns['mod']['factor_names']] = sp_adata.obsm['q05_cell_abundance_w_sf']
    sp_adata.obs[sp_adata.uns['mod']['factor_names']].to_csv(OutputFile)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-FunctionFlag', dest='FunctionFlag', default="", type=str, help='Function aim')
    parser.add_argument('-sTcell2location_scRNAFold', dest='sTcell2location_scRNAFold', default="", type=str, help='scRNA folder 10 X format')
    parser.add_argument('-sTcell2location_sTRNAFold', dest='sTcell2location_sTRNAFold', default="", type=str, help='spatial folder')
    parser.add_argument('-sTcell2location_OutputFile', dest='sTcell2location_OutputFile', default="", type=str, help='output file')
    args = parser.parse_args()

    if args.FunctionFlag == "sTcell2location":
        sTcell2location(args.sTcell2location_scRNAFold,args.sTcell2location_sTRNAFold,args.sTcell2location_OutputFile)

if __name__ == '__main__':
    main()