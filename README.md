# Mus-musculus-Complement-Expression
> Generating UMAPs to explore the expression of the complement cascade in Mus musculus

## Methods

**BioMart Annotation**

Gene identifier mapping for _Mus musculus_ was performed using Scanpy’s built-in BioMart interface (scanpy.queries.biomart_annotations) in a Python 3.10 environment (Wolf _et al._, 2018). We queried the Ensembl database (host: [www.ensembl.org](http://www.ensembl.org/)) for the external_gene_name and ensembl_gene_id attributes across the entire mouse genome. The returned records were converted into a pandas DataFrame and saved as a CSV file for subsequent analyses.

**Inflammatory Pathway Gene Mapping**

Mouse inflammation pathway gene sets were retrieved via the KEGG REST API using the **bioservices** Python package (Cokelaer _et al._, 2013). We initialized a KEGG client for _Mus musculus_ (organism = "mmu") and specified the complement cascade pathway (mmu04610). For each pathway, raw entries were parsed to extract gene symbols, which were collated into (1) a comprehensive pathway–gene mapping file, (2) a deduplicated gene list, and (3) individual pathway-specific gene lists. Outputs were written as CSV files into a timestamped directory for record keeping.

**Dataset acquisition and preprocessing**  
We used the **gget** Python package (Luebbert et al., 2023) to download and aggregate three curated mouse brain single-cell RNA-seq datasets—_Thyroid hormone remodels cortex to coordinate body-wide metabolism and exploration_, _Tabula Muris Senis_, and _Humoral immunity at the brain borders in homeostasis_—from the CellxGene Census (census_version = "2025-01-30") using their dataset UUIDs. Queries were restricted to healthy (“normal”) mouse brain tissue and returned Ensembl-annotated counts. Filter thresholds of ≥ 3 counts per gene and ≥ 500 genes detected per cell were applied, and the raw counts were preserved in adata.raw.

**Dimensionality reduction and integration**  
The gget .h5ad object was loaded with Scanpy (Wolf _et al._, 2018), normalized data were subjected to principal component analysis (PCA) with zero-centering disabled (scanpy.tl.pca(zero_center=False)) and harmonized across datasets using the Harmony algorithm (Korsunsky _et al._, 2019) via sc.external.pp.harmony_integrate(…, key="dataset_id"). A k-nearest neighbors graph was constructed on the Harmony-corrected PC space (scanpy.pp.neighbors(use_rep="X_pca_harmony")), followed by Uniform Manifold Approximation and Projection (UMAP; McInnes _et al._, 2018) embeddings computed with scanpy.tl.umap(). Processed data were saved as an H5AD file for downstream visualization.

**References**

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. <https://doi.org/10.1093/bioinformatics/btac836>
- Korsunsky I, Millard N, Fan J, et al. _Fast, sensitive and accurate integration of single-cell data with Harmony_. Nat Methods. 2019;16(12):1289–1296.
- Wolf FA, Angerer P, Theis FJ. _SCANPY: large-scale single-cell gene expression data analysis_. Genome Biol. 2018;19(1):15.
- McInnes L, Healy J, Melville J. _UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction_. arXiv:1802.03426. 2018.
- Cokelaer T, Pultz D, Harder L, Serra-Musach J, Vicente R, Godzik A. _Bioservices: a common Python package to access biological Web Services programmatically_. Bioinformatics. 2013;29(24):3241–3242.
- Kanehisa M, Goto S. _KEGG: Kyoto Encyclopedia of Genes and Genomes_. Nucleic Acids Res. 2000;28(1):27–30.
