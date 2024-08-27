###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "02_runCibersortxFractions"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Run CibersortX Fractions and analyse results"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the RNA-seq metadata file
PATH_RNASEQ_METADATA = file.path(  PATH_EXPERIMENT_RAWDATA, 
                                   "Gene_expression_matrixes",
                                   "001_Iontorrent_RNAseq_211006_Metadata.csv")

# Path to the RNA-seq gene expression matices from project experiment
PATH_CIBERSORT_ANALYSIS_RESULT_CHIP1_REGEX = file.path( PATH_ANALYSIS_OUTPUT,
                                                  "{TISSUE}",
                                                 "Chip1", 
                                                 "CIBERSORTx_Adjusted.txt")

PATH_CIBERSORT_ANALYSIS_RESULT_CHIP2_REGEX = file.path( PATH_ANALYSIS_OUTPUT, 
                                                  "{TISSUE}",
                                                 "Chip2", 
                                                 "CIBERSORTx_Adjusted.txt")

# Path to the computed cell type signature
PATH_CIBERSORT_SIGNATURE_RESULT_REGEX = file.path( PATH_ANALYSIS_OUTPUT, 
                                             "{TISSUE}",
                                             "Chip1",
					     "CIBERSORTx_sigmatrix_Adjusted.txt")


