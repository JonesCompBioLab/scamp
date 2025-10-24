###############################################
import pandas as pd
import anndata as ad
import numpy as np

def read_adata(anndata_file) :
    return ad.read_h5ad(anndata_file)

# Get all the correct settings for visualization
def setup_anndata(adata, scamp_tsv, temp_folder, cn_threshold, cn_percentile_threshold) :
    
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
            new_row = {"gene_set_name" : 'ecDNA', "gene_set_description" : "Predicted ecDNA by scAmp", "gene_symbol" : gene, "gene_description" : ""}
            gene_set_df.loc[len(gene_set_df)] = new_row
    
    gene_set_df.to_csv(f"{temp_folder}/ecDNA_gene_set.csv", index = False)

    

    # cellxgene needs an embedding, make sure we have one
    get_umap(adata, temp_folder)

    # Add cell sets
    for index, row in gene_set_df.iterrows():
        gene = row['gene_symbol']
        
        counts = adata[:, gene].X
        counts = counts.toarray().flatten()
        percentile_thresh = np.percentile(counts, cn_percentile_threshold)
        
        adata.obs[f"~{gene}_ecDNA"] = (counts > cn_threshold) | (counts > percentile_thresh)


    ecDNA_cols = [col for col in adata.obs.columns if col.endswith("_ecDNA")]
    adata.obs["Number of ecDNA Positive Genes"] = adata.obs[ecDNA_cols].astype(int).sum(axis=1)
    cols = ["Number of ecDNA Positive Genes"] + [c for c in adata.obs.columns if c != "Number of ecDNA Positive Genes"]
    adata.obs = adata.obs[cols]

    adata.write(f"{temp_folder}/annotated_anndata.h5ad")


def setup_copynumber(copy_numbers_file, scamp_tsv, temp_folder, cn_threshold, cn_percentile_threshold) :
    # Create anndata from tsv
    counts_df = pd.read_csv(copy_numbers_file, sep='\t')
    var = pd.DataFrame({
        'idx': range(1, counts_df.shape[1] + 1)
    }, index=counts_df.columns)
    obs = pd.DataFrame(index=counts_df.index)

    adata = ad.AnnData(X=counts_df.values, obs=obs, var=var)
    adata.uns['X_name'] = 'GeneScoreMatrix'

    # Call rest of the pipeline
    setup_anndata(adata, scamp_tsv, temp_folder, cn_threshold, cn_percentile_threshold)


# Create a umap if one doesn't exist
def get_umap(adata, temp_folder) :
    if "X_umap" not in adata.obsm :
        # Check for other umap
        for key in list(adata.obsm.keys()):
            if "umap" in key.lower():
                print(f"Using {key} as umap...")
                adata.obsm["X_umap"] = adata.obsm[key]
                return

        # Otherwise generate one
        print("No X_umap Obsm Found - Creating X_umap...")
        import scanpy as sc
        sc.pp.neighbors(adata)
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
