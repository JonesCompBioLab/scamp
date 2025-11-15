###############################################
import pandas as pd
import anndata as ad
import numpy as np

def read_adata(anndata_file) :
    return ad.read_h5ad(anndata_file)

# Get all the correct settings for visualization
# Returns: None
#          writes anndata to file [temp_folder]/annotated_anndata.h5ad
#          writes gene set data to file [temp_folder]/ecDNA_gene_set.csv
def setup_anndata(adata, scamp_tsv, temp_folder, cn_threshold, cn_percentile_threshold, umap_name, expression_adata) :
    
    # Making gene set
    scamp = pd.read_csv(scamp_tsv, sep = "\t")
    gene_set_data = {
        "gene_set_name": [],
        "gene_set_description": [],
        "gene_symbol": [],
        "gene_description": []
    }   
    gene_set_df = pd.DataFrame(gene_set_data)

    # get only ecDNA rows
    for index, row in scamp.iterrows():
        gene = row['gene']
        if row['pred'] == True:
            # When we have gene clusters from scamp, we can use this
            # new_row = {"gene_set_name" : row["cluster"], "gene_set_description" : "Predicted ecDNA by scAmp", "gene_symbol" : gene, "gene_description" : ""}

            new_row = {"gene_set_name" : 'ecDNA', "gene_set_description" : "Predicted ecDNA by scAmp", "gene_symbol" : gene, "gene_description" : ""}
            gene_set_df.loc[len(gene_set_df)] = new_row
    
    gene_set_df.to_csv(f"{temp_folder}/ecDNA_gene_set.csv", index = False)

    # Add cell sets
    add_cell_sets(adata, gene_set_df, cn_threshold, cn_percentile_threshold)

    # Convert visualization to expression data
    adata.X = expression_adata.X.copy()
    adata.X = np.log1p(adata.X)
    adata.obsm = expression_adata.obsm.copy()

    # cellxgene needs an embedding, make sure we have one
    get_umap(umap_name, adata, temp_folder)

    adata.write(f"{temp_folder}/annotated_anndata.h5ad")

# Adds cell sets to anndata
def add_cell_sets(adata, gene_set_df, cn_threshold, cn_percentile_threshold):
    # Get list of genes
    gene_list = gene_set_df['gene_symbol'].tolist()

    # Only get ecDNA genes
    X_ecDNA = adata[:, gene_list].X
    X_ecDNA = X_ecDNA.toarray() if not isinstance(X_ecDNA, np.ndarray) else X_ecDNA

    # Get number of ecDNA positive genes
    percentile_thresholds = np.percentile(X_ecDNA, cn_percentile_threshold, axis=0)
    ecDNA_bool = (X_ecDNA > cn_threshold) | (X_ecDNA > percentile_thresholds)
    adata.obs["Number of ecDNA Positive Genes"] = ecDNA_bool.sum(axis=1)

    # Convert to log for the obs
    logX = np.log1p(X_ecDNA)

    for i, gene in enumerate(gene_list):
        adata.obs[f"{gene}_ecDNA"] = logX[:, i]
        # Ensure it goes to top of visualization
        adata.obs = adata.obs[[f"{gene}_ecDNA"] + [c for c in adata.obs.columns if c != f"{gene}_ecDNA"]]
    adata.obs = adata.obs[["Number of ecDNA Positive Genes"] + [c for c in adata.obs.columns if c != "Number of ecDNA Positive Genes"]]



    # Get only gene list as obs
    if len(gene_list) > 0 :
        print("Only leaving ecDNA postiive genes in var")
        to_keep = adata.var_names.intersection(gene_list)
        adata = adata[:, to_keep].copy()

# Converts expression tsv data into anndata
# Returns: expression anndata
def setup_expression(expression_data, cn_adata) :
    # Detect file extension
    expression_file_ext = expression_data.split('.')[-1]
    if expression_file_ext == "csv" :
        sep = ','
    elif expression_file_ext == "tsv" :
        sep = '\t'
    else :
        print(f"Unknown file type for {expression_data}... treating as csv")
        sep = ','
    
    # Read dataframe
    exression_df = pd.read_csv(expression_data, sep=sep)
    exression_df_subset = exression_df.loc[:, exression_df.columns.isin(cn_adata.var_names)]
    exp_adata = ad.AnnData(exression_df_subset)
    return exp_adata

# Converts copy number tsv to anndata
# Returns: copy number anndata
def setup_copynumber(copy_number_file) :
    # Detect file extension
    copy_number_file_ext = copy_number_file.split('.')[-1]
    if copy_number_file_ext == "csv" :
        sep = ','
    elif copy_number_file_ext == "tsv" :
        sep = '\t'
    else :
        print(f"Unknown file type for {copy_number_file}... treating as csv")
        sep = ','

    # Create anndata from tsv
    counts_df = pd.read_csv(copy_number_file, sep=sep)
    var = pd.DataFrame({
        'idx': range(1, counts_df.shape[1] + 1)
    }, index=counts_df.columns)
    obs = pd.DataFrame(index=counts_df.index)

    adata = ad.AnnData(X=counts_df.values, obs=obs, var=var)
    adata.uns['X_name'] = 'GeneScoreMatrix'

    return adata

# Create a umap if one doesn't exist
# Returns: None
#          Edits input adata
#          If new umap created, new anndata saved in {temp_folder}/umap_anndata.h5ad
def get_umap(umap_name, adata, temp_folder) :
    # If we found existing umap
    if umap_name in list(adata.obsm.keys()):
        print(f"Found UMAP in {umap_name}...")
        if (umap_name != "X_umap") :
            print(f"Renaming UMAP to X_umap for cellxgene")
            adata.obsm["X_umap"] = adata.obsm[umap_name]
        return


    # Start with pca then run UMAP
    print(f"No {umap_name} Obsm Found. - Creating UMAP in X_umap...")
    import scanpy as sc
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.umap(adata)
    adata.obsp.clear()
    adata.varm.clear()

    # Clean up for cellxgene
    for key in list(adata.obsm.keys()):
        if key not in ["X_umap"]:
            del adata.obsm[key]

    for key in list(adata.uns.keys()):
        if key not in ["X_name"]:
            del adata.uns[key]


    print(f"Saving anndata with umap in {temp_folder}/umap_anndata.h5ad")
    adata.write(f'{temp_folder}/umap_anndata.h5ad')
