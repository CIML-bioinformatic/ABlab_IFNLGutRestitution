# Version 1.0
# 02/18/2020

--------------------------------------------------------------------------------
CIBERSORTx High resolution
--------------------------------------------------------------------------------
CIBERSORTxHiRes imputes sample-level gene expression variation of distinct cell
types from a collection of bulk tissue transcriptomes (Newman et al., submitted).
Unlike CIBERSORTxGEP, the output is an expression matrix for each cell type 
rather than a single transcriptome profile. CIBERSORTxHiRes is useful for
exploring cell type expression variation without prior knowledge of biological
or functional groupings (e.g. relating cell type specific gene expression to 
survival). 

--------------------------------------------------------------------------------
INSTALLATION
--------------------------------------------------------------------------------
Install Docker Desktop (https://www.docker.com/products/docker-desktop)
Open Docker Desktop on your local computer and log in.
Then, you open a terminal window, and you type the following command:
> docker pull cibersortx/hires

The next thing you need is a token that you will provide every time you run the
CIBERSORTx executables. You can obtain the token from the CIBERSORTx website
(https://cibersortx.stanford.edu/getoken.php).

Please note that each token is uniquely tied to your user account, and tokens are
good for a specific time interval from date of request, so you will need to
request a new token when an existing one has expired.

Once you have pulled the CIBERSORTx executable from Docker, and you have obtained
a token from the CIBERSORTx website, you now have access to the CIBERSORTx HiRes
executable and can run it following the instructions below. 

--------------------------------------------------------------------------------
INSTRUCTIONS & EXAMPLES
--------------------------------------------------------------------------------
These instructions assume that you already have completed the tutorials on the 
the CIBERSORTx website (https://cibersortx.stanford.edu/tutorial.php), and have
learned how to set up a CIBERSORTx job. A typical CIBERSORTxHiRes run may take
from a few seconds with a gene subset file to potentially over an hour depending
on the number of genes and the number of samples.

Download the file packages for the examples from Download in the Menu tab on the
CIBERSORTx website.

When running the following commands, you should be able to reproduce the results 
of the examples on the website. For further information about these examples, refer
to the tutorial on the website: https://cibersortx.stanford.edu/tutorial.php.

To run CIBERSORTx HiRes locally on your computer using docker, navigate to the
directory containing the files you wish to analyze. 

Next, you have to "bind mount" directories so that they can can be accessed within
the docker container. This is done using the following command when setting up the
CIBERSORTx job:

> docker run -v {dir_path}:/src/data -v {dir_path}:/src/outdir cibersortx/hires [Options]

Please note that you have to provide absolute paths. To better understand how to
bind mount directories and setting up the CIBERSORTx job using docker, please follow
the following examples: 

# Group Level GEPs - FL (Fig. 3b-f):
# This examples imputes cell type specific gene expression profiles from bulk
# follicular lymphoma samples profiled on microarray, using the signature matrix
# LM22 collapsed to 4 major cell types. In addition the results are compared
# to ground truth reference profiles obtained from FACS-sorted cell subsets.

docker run -v absolute/path/to/input/dir:/src/data -v absolute/path/to/output/dir:/src/outdir cibersortx/hires --username email_address_registered_on_CIBERSORTx_website --token token_obtained_from_CIBERSORTx_website --mixture Fig4a-arrays-SimulatedMixtures.MAS5.txt --sigmatrix Fig4a-LM4.txt --classes Fig4a-LM4-mergedclasses.txt --window 20 --QN TRUE

# Group Level GEPs - NSCLC (Fig. 3g):
# This examples imputes cell type specific gene expression profiles from bulk
# NSCLC samples profiled by RNA-Seq and compares the results to ground truth
# reference profiles obtained from FACS-sorted cell subsets.

docker run -v absolute/path/to/input/dir:/src/data -v absolute/path/to/output/dir:/src/outdir cibersortx/hires --username email_address_registered_on_CIBERSORTx_website --token token_obtained_from_CIBERSORTx_website --mixture SuppFig11-DLBCL_CHOP_Lenz-arrays-bulktumors.MAS5.txt --sigmatrix  LM22.txt --classes SuppFig11-LM22_10_merged_classes.txt --subsetgenes SuppFig11-DLBCL-GCBABC-genes --QN TRUE


--------------------------------------------------------------------------------
Options for CIBERSORTxHiRes.R:
--------------------------------------------------------------------------------
For further details about the different options, refer to the tutorials on the 
website: https://cibersortx.stanford.edu/tutorial.php
For further details about the input file format, refer to the upload file page on
the website: https://cibersortx.stanford.edu/upload.php

--mixture       <file_name>  Mixture matrix [required]
--sigmatrix     <file_name>  Signature matrix [required]
--classes       <file_name>  Cell type groupings [optional, default: signature matrix labels]
--cibresults    <file_name>  Previous CIBERSORTx cell fractions 
                        [default: run CIBERSORT]
--label         <char>  Sample label [default: none]
--rmbatchBmode       <bool>  Run B-mode batch correction [default: FALSE]
--rmbatchSmode       <bool>  Run S-mode batch correction [default: FALSE]
--sourceGEPs    <file_name>  Signature matrix GEPs for batch correction
                        [required if rmbatch=T; default: sigmatrix]
--groundtruth   <file_name>  Ground truth GEPs [same labels as classes] [optional; 
                        default: none]
--threads       <int>   Number of parallel processes [default: No. cores - 1]
--QN            <bool>  Run quantile normalization [default: FALSE]
--nsampling     <int>   Number of subsamples for NNLS [default: 30]
--degclasses    <file_name>  Run on two classes, specified by 1, 2, 0=skip 
                        [default: none]
--redo          <bool>  Redo transcriptome imputation [default: FALSE]
--redocibersort <bool>  Redo CIBERSORT [default: FALSE]
--useadjustedmixtures <bool>  If doing batch correction, use adjusted mixtures for GEP imputation [default: TRUE]
--variableonly  <bool>  Write only genes that show variation in expression to final output [default: FALSE]

--------------------------------------------------------------------------------
CIBERSORTx High Resolution Mode - OUTPUT FILES
--------------------------------------------------------------------------------

MAIN OUTPUT FILES
--------------------------------------------------------------------------------
The main result of CIBERSORTx High Resolution is a set of expression matrix .txt
files and heatmaps for each individual cell type, showing the cell-type specific
expression of individual genes at the sample level.

The name of each file contains the name of the cell type, followed by the
window size used for the job (e.g. CIBERSORTxHiRes_job1_Mastcells_Window40.txt).

The "1" values in the expression matrix txt files are genes with insufficient
evidence of expression (these genes are either not expressed or have inadequate
statistical power to be imputed).

The NA values are genes that have inadequate statistical power to be imputed.

If you have provided a gene subset list, these files will have all the genes in
the list that are found in the original mixture file given as input. If some
genes are still missing, this could be due to different annotations or gene
symbols between the gene subset list and the mixture file.

ADDITIONAL OUTPUT FILES
--------------------------------------------------------------------------------
CIBERSORTx High Resolution Mode runs both CIBERSORTxFractions and CIBERSORTxGEP,
and a set of output files is generated from each analysis.

Output from CIBERSORTxFractions:
--------------------------------------------------------------------------------
- *_Fractions.txt: file enumerating the fractions of the different cell types in
bulks samples. 

If batch correction is performed (B-mode or S-mode):
- *_Adjusted.txt: file enumerating the fractions of the different cell types in 
mixture samples after batch correction.
- *_Mixtures_Adjusted.txt: the input mixture file after batch correction

If S-mode batch correction is performed:
- [signature matrix filename]_Adjusted: the signature matrix after batch correction 

Output from CIBERSORTx Group Mode:
--------------------------------------------------------------------------------
CIBERSORTxGEP imputes a representative transcriptome profile for each cell type
in the signature matrix from a group of bulk tissue transcriptomes  (Newman et 
al., submitted). This mode is useful for learning context-dependent changes in 
a cell's expression profile (e.g., responders vs. non-responders to a therapy). 

The main result of CIBERSORTx Group Mode is a file for cell-type specific gene
expression profiles where genes have been filtered out using a threshold to
eliminate unreliably estimated genes for each cell type. This file is called
*_GEPs_Filtered.txt.

*_GEPs.txt is the file for cell-type specific gene expression profiles where no
filtering was done.

The different statistics used for the filtering are saved in the files listed
below. Refer to Supplementary Note in Newman et al. (submitted) for further
details:
- *_GEPs_StdErrs.txt: analytically derived standard errors of the regression
coefficients.
- *GEPs_Pvals.txt: p-values used to determine the significance of the
regression coefficients.
- *GEPs_Qvals.txt: adjusted p-values (q-values) after multiple hypothesis
testing using the Benjamini-Hochberg method.
- *_GEPs_CV.txt: To further reduce confounding noise, genes were filtered based
on their geometric coefficient of variation (geometric c.v.), which are the
values listed in this file, calculated using the natural logarithm of
subsampled regression coefficients. The geometric c.v. were used to
determined the adaptive cell-type specific noise threshold.
- *_GEPs_ThresholdPlots.pdf: plot illustrating the adaptive noise threshold 
used for filtering.
- *_Weights.txt: the fractions of the different cell types after
merging them into major classes according to the merged classes file. 

If ground truth was given as input:
- *_CrossCorrelationMatrix.pdf: plots showing the correlation between estimated
genes and ground truth for each cell types. for all genes (GEP), and for the
signature matrix genes (SM). The corresponding *_CrossCorrelationMatrix.txt file
is also given as output.
- *_ScatterPlots.pdf: scatterplots showing the estimated expression values (non-
zero only) versus the observed expression values for the whole gene expression 
profile (GEP) and for the signature matrix genes (SM) after noise filtering.
- *_SM_GEPS_HeatMap.png: Heatmap illustrating the CiBERSORTx imputed gene
expression values for the signature matrix genes (y axis), compared to ground
truth.
- *GEPs_Stats.txt: set of benchmark statistics used to compare with ground truth.

--------------------------------------------------------------------------------
CIBERSORTx HiRes Source Code
--------------------------------------------------------------------------------
Source code is available for academic non-profit use upon reasonable request (mailto:cibersortx@gmail.com).
