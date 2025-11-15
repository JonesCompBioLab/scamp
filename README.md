# scAmp
scAmp (single-cell Amplicon) is a python-based workflow for detecting and
analyzing focal amplifications from single-cell data.

It consists of three main modules (in progress):

* copy-number inference: a set of modules for inferring copy-numbers from single-cell data
* ecDNA detection: detection of extrachromsomal or chromosomal DNA amplifications
* clonal analysis: Analysis of clonal history with respect to amplifications

## Installation

You can install scAmp by cloning this directory and running `pip install .` This should install scAmp and all python dependencies. scAmp should be used with python version >= 3.9.

In addition, if you'd like to use the ATAC CNV module, you should install the following R packages:

* [ArchR](https://github.com/GreenleafLab/ArchR)
* SummarizedExperiment
* dplyr

scAmp's ATAC CNV module was tested with R version 4.3.2 (and should work with later versions).

To use the visualization module, you should run `pip install cellxgene`

## Running scAmp from the command line

You can invoke `scamp` modules from the command line by running

`scamp [module] [arguments]`

Currently there are three command line modules you can run:

* `scamp atac-cnv`: Computes copy-numbers across genes from a scATAC fragments files.
* `scamp predict`: Predicts ecDNA status from copy-number data.
* `scamp visualize`: Visualizes results on cellxgene after running `predict`

You can look at usage instructions for these modules by running `scamp [module] --help`.

## Visualization

After running visualization once on a dataset, it will generate two files in your specified temp folder:

* `annotated_anndata.h5ad`: Anndata with ecDNA cell set annotations
* `ecDNA_gene_set.csv`: Gene sets for ecDNA

If you want to visualize the same dataset again, you can then just run `cellxgene launch [temp_folder]/annotated_anndata.h5ad --gene-sets-file [temp_folder]/ecDNA_gene_set.csv --open`

Necessary inputs are scamp predict results, expression data, and copy number data. 

Supported data types for expression and copy number data are cellranger MEX, anndata, csv, and tsv. To use MEX, please give a folder as an argument. To use csv or tsv, please give a file with the extension `.tsv` or `.csv`. To use anndata, ensure the file has the extension `.h5ad`.

If your umap is not named "X_umap", provide the name to `--umap-name`. Otherwise, a umap will be created automatically and an anndata will be saved with the new umap in `[temp_folder]/umap_anndata.h5ad`

## Pretrained models

Though you can train new models using scAmp, we also provide pre-trained models in the `./pretrained_models` directory.

You can pass a path to the pretrained model directly to `scamp predict`. Otherwise, if you are using scAmp interactively, you can load this in as so:

```
from scamp import models

saved_model_path = "./pretrained_models/scamp_model_1.0"
pretrained_model = models.SCAMP.load(pretrained_model_path)
```
