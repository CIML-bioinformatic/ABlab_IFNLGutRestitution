#!/bin/bash

# CIBERSORTx - enumerate cell composition in bulk genomic profiles
# 
# For instructions and terms of use, see cibersortx.stanford.edu
# 
# Usage:
# docker run <bind_mounts> cibersortxfractions [Options] 
# 
# Bind Mounting:
# > 2 directories must be bind mounted to be accessed within the docker container: 
#     1) Input file dir 
#         > Format: -v {dir_path}:/src/data 
#     2) Outdir 
#         > Format: -v {dir_path}:/src/outdir 
# > Note: Absolute paths required
# 
# Authorization Parameters:
# --username      <string>  Email used for login to cibersortx.stanford.edu
# --token         <string>  Token associated with current IP address (generated on website)
# 
# Primary Options:
# --mixture       <file_name>  Mixture matrix [required for running CIBERSORTx, optional for creating a custom signature matrix only]
# --sigmatrix     <file_name>  Signature matrix [required: use preexisting matrix or create one]
# --perm          <int>   No. of permutations for p-value calculation [default: 0]
# --label         <char>  Sample label [default: none]
# --rmbatchBmode  <bool>  Run B-mode batch correction [default: FALSE]
# --rmbatchSmode  <bool>  Run S-mode batch correction [default: FALSE]
# --sourceGEPs    <file_name>  Signature matrix GEPs for batch correction [default: sigmatrix]
# --QN            <bool>  Run quantile normalization [default: FALSE]
# --absolute      <bool>  Run absolute mode [default: FALSE]
# --abs_method    <char>  Pick absolute method ['sig.score' (default) or 'no.sumto1']
# --verbose       <bool>  Print verbose output to terminal [default: FALSE]
# 
# Options for creating a custom signature matrix:
# --refsample     <file_name>  Reference profiles (w/replicates) or labeled scRNA-Seq data [required]
# --phenoclasses  <file_name>  Cell type classes [required, if single_cell = FALSE]
# --single_cell   <bool>  Create matrix from scRNA-Seq data [default: FALSE]
# --G.min         <int>   Minimum number of genes per cell type in sig. matrix [default: 50, if single_cell = TRUE: 300]
# --G.max         <int>   Maximum number of genes per cell type in sig. matrix [default: 150, if single_cell = TRUE: 500] 
# --q.value       <int>   Q-value threshold for differential expression [default: 0.3, if single_cell = TRUE: 0.01] 
# --filter        <bool>  Remove non-hematopoietic genes [default: FALSE] 
# --k.max         <int>   Maximum condition number [default: 999] 
# --remake        <bool>  Remake signature gene matrix [default: False] 
# --replicates    <int>   Number of replicates to use for building scRNAseq reference file [default: 5] 
# --sampling      <float> Fraction of available single cell GEPs selected using random sampling [default: 0.5] 
# --fraction      <float> Fraction of cells of same identity showing evidence of expression [default: 0.75] 

# --------------------------------------

# DEFINE THE INPUT DIR and OUTPUT DIR

INPUT_DIR="/mnt/DOSI/ABLAB/BIOINFO/Project/001_IFNL_Recovery/001_Iontorrent_RNAseq_211006_Deconvolution/05_Output/01_PrepareData"
OUTPUT_DIR="/mnt/DOSI/ABLAB/BIOINFO/Project/001_IFNL_Recovery/001_Iontorrent_RNAseq_211006_Deconvolution/05_Output/02_runCibersortxFractions"

# DEFINE THE EMAIL AND TOKEN ALLOWED TO RUN CIBERSORT

EMAIL="lionel.spinelli@univ-amu.fr"
TOKEN="340a6d2a39e309b2b0acafbc95466382"

# DEFINE THE DEFAULT OPTIONS

FRACTION="0"
RMBATCH_MODE="TRUE"
SINGLE_CELL_MODE="TRUE"

TISSUE_LIST=("AllTissues" "Ileum" "Duodenum" "Jejunum")

for TISSUE in "${TISSUE_LIST[@]}"
do
  
  # DEFINE THE REFERENCE SINGLE-CELL FILE TO EXTRACT THE SIGNATURE FROM
  
  REFERENCE_FILE="regional_cell_sampling_Log2TPM_fixed_nolog2_nocellid_${TISSUE}.txt"
  
  # RUN THE ANALYSIS FOR THE FIRST MIXTURE
  # --------------------------------------
  
  MIXTURE_FILE="Auto_Zanoni_LP1_093021_Trx_MOUSE_s5-torrent-server-vm_247.rpm.bcmatrix.txt"
  OUTPUT_SUB="${TISSUE}/Chip1"
  
  mkdir -p ${OUTPUT_DIR}/${OUTPUT_SUB}
  #docker run --rm -v ${INPUT_DIR}:/src/data -v ${OUTPUT_DIR}/${OUTPUT_SUB}:/src/outdir ablab_iflnrecovery_cibersort --username $EMAIL --token $TOKEN --single_cell $SINGLE_CELL_MODE --refsample $REFERENCE_FILE --mixture $MIXTURE_FILE --verbose TRUE > ${OUTPUT_DIR}/${OUTPUT_SUB}/cibersort_run.log
  docker run --rm -v ${INPUT_DIR}:/src/data -v ${OUTPUT_DIR}/${OUTPUT_SUB}:/src/outdir ablab_iflnrecovery_cibersort --perm 1000 --username $EMAIL --token $TOKEN --single_cell $SINGLE_CELL_MODE --refsample $REFERENCE_FILE --mixture $MIXTURE_FILE --fraction $FRACTION --rmbatchSmode $RMBATCH_MODE --verbose TRUE > ${OUTPUT_DIR}/${OUTPUT_SUB}/cibersort_run.log
  
  # RUN THE ANALYSIS FOR THE SECOND MIXTURE
  # ---------------------------------------
  
  MIXTURE_FILE="Auto_Zanoni_LP2_093021_Trx_MOUSE_s5-torrent-server-vm_249.rpm.bcmatrix.txt"
  OUTPUT_SUB="${TISSUE}/Chip2"
  
  mkdir -p ${OUTPUT_DIR}/${OUTPUT_SUB}
#  docker run --rm -v ${INPUT_DIR}:/src/data -v ${OUTPUT_DIR}/${OUTPUT_SUB}:/src/outdir ablab_iflnrecovery_cibersort --perm 1000 --username $EMAIL --token $TOKEN --single_cell $SINGLE_CELL_MODE --refsample $REFERENCE_FILE --mixture $MIXTURE_FILE --verbose TRUE > ${OUTPUT_DIR}/${OUTPUT_SUB}/cibersort_run.log
  docker run --rm -v ${INPUT_DIR}:/src/data -v ${OUTPUT_DIR}/${OUTPUT_SUB}:/src/outdir ablab_iflnrecovery_cibersort --perm 1000 --username $EMAIL --token $TOKEN --single_cell $SINGLE_CELL_MODE --refsample $REFERENCE_FILE --mixture $MIXTURE_FILE --fraction $FRACTION --rmbatchSmode $RMBATCH_MODE --verbose TRUE > ${OUTPUT_DIR}/${OUTPUT_SUB}/cibersort_run.log

done
