# IFNL_Recovery_scRNAseq

Source code and scripts for bioinformatics analyses of 10x Genomics single-cell dataset.

## What are 'globalParams.R' and 'analysisParams.R' file ?

Both files are R scripts that aim to be sourced by other (analysis) scripts to provide pre-configured variables.

The file 'globalParams.R' defines general project-specific variables (mostly paths to input data and result folders).
The file 'analysisParams.R' is located in the folder of an analysis step and defines parameters specific to this analytic step (e.g. 01_QC).

## What is 'launch_reports_compilation.R' file ?

This file is a helper that starts the rendering of associated Rmd files. It is found next to Rmd files.
It is in charge of loading variables from 'globalParams.R' and 'analysisParams.R' before starting to render the associated Rmd file.

## Example to execute one step

### Note

Most scripts reload data from results of previous steps using RDS files (and sometimes csv files). 
One must download input files expected by the specific analysis step he is planning to execute (see `analysisParams.R` to get a list of expected inputs). 

Results for all individual steps were uploaded to "processed data" repository (links in parent folder).
Therefore intermediate results are already available for direct download, and it is possible to execute individual analysis step without having to reprocess all upstream processing.

### Preparation

After cloning repository, one must update path to the project in files `globalParams.R`. 
These variables are used in all other scripts to define relative path to data. 

Execution environment must be imported from exported docker images (tar.gz) available on "processed data" repository (links in parent folder).
Import appropriate docker image with: docker load < fileName.tar.gz

For executing first analysis step (cellranger), fastq files must be donloaded from "raw data" repository (links in parent folder).

### Execution

Change current directory to the cloned repository, and execute:
```
docker run --rm -v /mnt:/mnt \
           -e PASSWORD=Pass -e USER=$(whoami) \
           -e USERID=$(id -u) -e GROUPID=$(id -g) \
           rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
           Rscript ./01_GlobalHeterogeneity_priorClustering/launch_reports_compilation.R"
```

Replace `./01_GlobalHeterogeneity_priorClustering/launch_reports_compilation.R` by the appropriate script to be executed, and `rfenouil/r411_tidyverse_seurat4` by the imported docker image.


