# scAmp
[![lint](https://github.com/JonesCompBioLab/scamp/actions/workflows/lint.yml/badge.svg)](https://github.com/JonesCompBioLab/scamp/actions/workflows/lint.yml)

scAmp (single-cell Amplicon) is a python-based workflow for detecting and
analyzing focal amplifications from single-cell data.

It consists of three main modules (in progress):

* copy-number inference: a set of modules for inferring copy-numbers from single-cell data
* ecDNA detection: detection of extrachromsomal or chromosomal DNA amplifications
* ecDNA deconvolution: deconvolution of ecDNA+ genes into ecDNA species

The latter two processes have been improved with scRASCAL (single-cell Rare Amplicon SubClone Analysis through Latent factorization) which builds on the original scAmp predictive framework
to detect rare, subclonal ecDNA.  

## Installation

You can install scAmp by cloning this directory and running `pip install .` This should install scAmp and all python dependencies. scAmp should be used with python version >= 3.9.

We support CNV quantification in both R and python. If you choose to use the R version (`atac-cnv-legacy`), you'll need to install the following R packages:

* Rutils
* SummarizedExperiment
* dplyr

scAmp's ATAC CNV module was tested with R version 4.3.2 (and should work with later versions).

To use the visualization module, you should run `pip install cellxgene`

## Running scAmp from the command line

For a basic tutorial of the `scamp` classification pipeline, you can find an example [here](notebooks/intro_to_scamp.ipynb).

You can also invoke `scamp` modules from the command line by running

`scamp [module] [arguments]`

Currently there are three command line modules you can run:

* `scamp atac-cnv`: Computes copy-numbers across genes from a scATAC fragments files. Supports genome builds hg38, hg19, mm10, and mm39. 
* `scamp atac-cnv-legacy`: Computes copy-numbers across genes using a legacy implementation in R. Supports genome builds hg38, hg19, and mm10. 
* `scamp predict`: Predicts ecDNA status from copy-number data.
* `scamp visualize`: Visualizes results on cellxgene after running `predict`

You can look at usage instructions for these modules by running `scamp [module] --help`.

### Running scRASCAL with `scamp predict`

`scamp predict` calls each gene ecDNA+ or ecDNA− from a cell-by-gene table of copy numbers, then groups the
ecDNA+ genes into ecDNA species (deconvolution). The usage is:

```
scamp predict [OPTIONS] OUTPUT_DIR COPY_NUMBERS [MODEL_FILE]
```

For example, to run scRASCAL on the COLO320DM example data in this repository:

```
scamp predict colo320dm_results data/colo320dm_cna_subset.tsv pretrained_models/scamp_model_1.0 \
    --predict-method GMM \
    --deconvolution-method auto \
    --seed 10
```

This writes three files to `colo320dm_results/`, named after the input file (here `colo320dm_cna_subset`):

| File | Contents |
|---|---|
| `ecDNA_preds_<input>.tsv` | One row per gene that passed `--filter-copy-number`: `gene`, `mean` and `var` (the gene's copy-number mean and variance across cells), `prediction` (`True` = ecDNA+), `Species` (comma-separated ecDNA species, or `None`), and columns for the prediction method: `Probability` (SCAMP), `Probability` and `Weight` (GMM), or `NumberOfPositiveCells` (KNN). |
| `cell_by_ecDNA_<input>.tsv` | One row per cell (cell barcodes in the first column), one column per ecDNA species: the estimated amount of each species in each cell. |
| `log.log` | The settings used and progress messages. |

For a walkthrough of how the settings change the calls, and how to tune them for rare ecDNA subclones, see the
[scRASCAL tutorial](notebooks/scRASCAL_tutorial.ipynb).

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `OUTPUT_DIR` | required | Directory for the output files. Created if it does not exist. |
| `COPY_NUMBERS` | required | Cell-by-gene copy-number matrix: a `.tsv` or `.csv` file (cells as rows, genes as columns), an AnnData `.h5ad` file, or a cellranger MEX folder. |
| `MODEL_FILE` | `scamp_model_1.0` (bundled) | Path to a saved scAmp model directory. |

#### Input and output options

| Option | Default | Description |
|---|---|---|
| `--whitelist-file` | none | File with one cell barcode per line and no header. Only these cells are used. |
| `--cNMF-log-dir` | `./temp` | Directory for cNMF's intermediate files. |
| `--verbose` | off | Write more detail about the species deconvolution to `log.log`. |

#### Prediction options

| Option | Default | Description |
|---|---|---|
| `--predict-method` | `GMM` | How genes are called: `SCAMP` (the original scAmp model, scoring all of a gene's cells at once), `GMM` (split the cells into distributions with a Gaussian mixture model and score each distribution), or `KNN` (split the cells the same way, then score each cell's neighbourhood). |
| `--decision-rule` | `0.5` | Probability cutoff for calling a gene, distribution (GMM) or neighbourhood (KNN) ecDNA+. |
| `--filter-copy-number` | `2.5` | Genes whose mean copy number is at or below this are not called. |
| `--min-copy-number` | `2.0` | Cells with copy number at or below this are left out. SCAMP: when the model scores the gene. GMM and KNN: before the GMM fit. |
| `--max-percentile` | `99.9` | Cells above this percentile of the gene's copy numbers are left out, in the same places as `--min-copy-number`. |
| `--no-trim-before-gmm` | off | GMM and KNN: do not apply `--min-copy-number` and `--max-percentile` before the GMM fit. |

#### GMM and KNN options

| Option | Default | Description |
|---|---|---|
| `--one-thresh` | `1.15` | How readily the GMM settles on a single distribution. Higher values give one distribution more often. |
| `--kneedle-coeff` | `3` | Higher values make the GMM choose fewer distributions. |
| `--min-weight` | `0.005` | Distributions holding this fraction of the cells or less are ignored. |
| `--var-scale` | `0.5` | Each distribution's variance is scaled by its weight to the power of this. Higher values give fewer false positives. |
| `--no-single-dist-scamp` | off | By default, when only one distribution holds more than `--min-weight` of the cells, the gene is scored on all of its cells as published scAmp would score it. With this flag, that distribution is scored like any other. |
| `--n-draws` | `1` | GMM: number of draws of 1,000 cells sampled from each distribution. The call uses the mean probability over the draws. |
| `--k` | `100` | KNN: number of nearest neighbours in each cell's neighbourhood. |
| `--knn-min-mean` | `2.5` | KNN: a neighbourhood is scored only if its mean copy number is above this (`-inf` turns this off). |
| `--knn-min-var` | `10` | KNN: a neighbourhood is scored only if its copy-number variance is above this (`-inf` turns this off). |

#### Deconvolution options

| Option | Default | Description |
|---|---|---|
| `--deconvolution-method` | `auto` | How ecDNA+ genes are grouped into species: `hier` (hierarchical clustering), `cNMF` (consensus non-negative matrix factorization, which lets a gene belong to more than one species), or `auto` (cNMF if genes appear to be shared between species, otherwise hierarchical). |
| `--num-species` | chosen automatically | Fix the number of ecDNA species. |
| `--max-species` | number of genes − 1 | Largest number of ecDNA species to consider. |
| `--hier-ddist` | `1.3` | Distance given to null nodes in hierarchical clustering, used to choose the number of species. |
| `--overlap-thresh` | `0.55` | `auto`: threshold for detecting genes shared between species. Lower values choose cNMF less often. |
| `--error-w` | `0.25` | cNMF: weight of the reconstruction error when choosing the number of species. |
| `--deconvolution-score-cutoff` | `3` | cNMF: a gene is assigned to a species if its score for that species is above this. |
| `--ecDNA-remove-cutoff` | `0.01` | cNMF: species whose highest spectra value is below this are removed. |
| `--cNMF-num-iter` | `50` | cNMF: number of factorization runs for each number of species tried. |
| `--cNMF-density-threshold` | `0.1` | cNMF: density threshold for the consensus of the runs. |
| `--seed` | random | Random seed for cNMF, for reproducible species. |

## Visualization

Necessary inputs are scamp predict results, expression data, and copy number data. 

Supported data types for expression and copy number data are cellranger MEX, anndata, csv, and tsv. To use MEX, please give a folder as an argument. To use csv or tsv, please give a file with the extension `.tsv` or `.csv`. To use anndata, ensure the file has the extension `.h5ad`.

If your umap is not named "X_umap", provide the name to `--umap-name`. Otherwise, a umap will be created automatically and an anndata will be saved with the new umap in `[temp_folder]/umap_anndata.h5ad`

After running visualization once on a dataset, it will generate two files in your specified temp folder:

* `annotated_anndata.h5ad`: Anndata with ecDNA cell set annotations
* `ecDNA_gene_set.csv`: Gene sets for ecDNA

If you want to visualize the same dataset again, you can then just run `cellxgene launch [temp_folder]/annotated_anndata.h5ad --gene-sets-file [temp_folder]/ecDNA_gene_set.csv --open`

## Pretrained models

Though you can train new models using scAmp, we also provide pre-trained models in the `./pretrained_models` directory.

You can pass a path to the pretrained model directly to `scamp predict`. Otherwise, if you are using scAmp interactively, you can load this in as so:

```
from scamp import models

saved_model_path = "./pretrained_models/scamp_model_1.0"
pretrained_model = models.SCAMP.load(pretrained_model_path)
```

## Citing scAmp

You can cite scAmp with the following paper:

Jones, M.G., Weiser, N.E., Hung, K.L. et al. scAmp enables focal gene amplification analysis from single-cell data. Nat Commun (2026). https://doi.org/10.1038/s41467-026-76417-3
