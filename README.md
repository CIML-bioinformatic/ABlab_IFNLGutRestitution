# ABlab IFNLGutRestitution

This repository contains source code for bioinformatics analyses presented in corresponding article.

## Reference

## Authors

## Abstract

## Repository content

Two distinct dataset from different RNA-sequencing technologies ('ion-torrent bulk' and '10x genomics single-cell') were analysed separately.
This repository has two independent project folders, each containing source code used for corresponding analyses described in the article.

Both projects share a common organisation of the folder structure, and execution strategy to facilitate reproductibility.

In brief, the folder structure is hierarchically organised by:
* Project ('ion-torrent bulk' or '10x genomics single-cell')
* Experiment (Analysis of individual replicates, or merged)
* Analysis step (prefixed by a number ordering the sequential processings of current experiment)

Variables used in R scripts are generally defined outside the script, in files suffixed with `*...Params.R`. 
`GlobalParams.R` define common variables for all scripts in an experiment folder, while `AnalysisParams.R` defines variables specific for each analysis step.
See readme in each project for more details on how to load these variables and compile reports automatically (helper script).

## External ressources

### RAW data

RAW sequencing data (fastq files) required as starting point for this analysis were uploaded to Gene Expression Omnibus:
* Bulk (ion-torrent):
* Single-cell (10x genomics): [GSE246333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246333)

### Processed data

Resulting reports and files for all analysis steps were uploaded to [recherche.data.gouv.fr](recherche.data.gouv.fr):
* Bulk (ion-torrent):
* Single-cell (10x genomics): [doi/10.57745/KPIAYZ](https://doi.org/10.57745/KPIAYZ)

## Project naming

While the full project name is `IFNL Gut Restitution`, the normalized name used for analyses is `IFNL_Recovery`:
* `001_IFNL_Recovery` for ion-torrent bulk sequencing
* `002_IFNL_Recovery_scRNAseq` for 10x genomics single-cell

Project name appears as such in `globalParams.R` file. 
This name is used extensively in the scripts (mostly for input path and output filenames) so it should be preserved to facilitate reproduction of results.

## Reproducing results

To execute analyses, one needs to download Docker images from corresponding "Processed data" repository (source for Dockerfile available in github subfolder). 
It is recommended to clone the repository, and modify `globalParams.R` files to match the path where the repository has been cloned. 

Then, each analysis step should provide a script that can be executed within the appropriate container. 
For R analyses, executing `launch_reports_compilation.R` from R (or using `Rscript`) takes care of loading variables from `*...Params.R` files, and rendering associated Rmd report (and all other files) in output folder.

