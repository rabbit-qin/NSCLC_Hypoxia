# The importance of hypoxia on alternative polyadenylation and tumor microenvironmental alterations in non-small cell lung cancer
### Data:
### APAEventsDiff`:
- The differentially alternative polyadenylation (APA) events calculated in The Cancer Genome Atlas (TCGA) and the Cancer Cell Line Encyclopedia (CCLE) datasets.

### NCScore:
- `Hypoxia_scores.gmt`: The hypoxia gene set used in the NC method;
- `Hypoxia_scoresTCGA` and `Hypoxia_scoresCCLE`: The hypoxia scores in the TCGA and CCLE datasets calculated by the NC method;
- `NC2NMF_cluster3TCGA.csv` and `NC2NMF_cluster3CCLE.csv`: By integrating the two methods, classify the samples in the TCGA and CCLE datasets.

### NMFScore:
- `NMF_geneset.csv`: The hypoxia gene set used in the NMF method;
- `cluster2Tcga.csv` and `cluster3Ccle.csv`: The hypoxia scores in the TCGA and CCLE datasets calculated by the NMF method;

### rawData:
- Gene expression data, PDUI value data, survival information data, and drug information data in the TCGA and CCLE datasets; PDUI value data in the single cell datasets.

### Codes:
### bulk&ccle:
- The code used for analyzing the TCGA and CCLE datasets.

### singleCell:
- The code used for analyzing the single cell datasets.

### spatialTranscriptome:
- The code used for analyzing the spatial transcriptome datasets.
