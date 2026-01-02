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

* Rutils
* SummarizedExperiment
* dplyr

scAmp's ATAC CNV module was tested with R version 4.3.2 (and should work with later versions).

## Running scAmp from the command line

For a basic tutorial of the `scamp` classification pipeline, you can find an example [here](notebooks/intro_to_scamp.ipynb).

You can also invoke `scamp` modules from the command line by running

`scamp [module] [arguments]`

Currently there are two command line modules you can run:

* `scamp atac-cnv`: Computes copy-numbers across genes from a scATAC fragments files.
* `scamp predict`: Predicts ecDNA status from copy-number data.

You can look at usage instructions for these modules by running `scamp [module] --help`.

## Pretrained models

Though you can train new models using scAmp, we also provide pre-trained models in the `./pretrained_models` directory.

You can pass a path to the pretrained model directly to `scamp predict`. Otherwise, if you are using scAmp interactively, you can load this in as so:

```
from scamp import models

saved_model_path = "./pretrained_models/scamp_model_1.0"
pretrained_model = models.SCAMP.load(pretrained_model_path)
```
