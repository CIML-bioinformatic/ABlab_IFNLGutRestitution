###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "01_prepareData"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Prepare the data from public repository for CibersortX analysis"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the public dataset that will be used as reference in the CibersortX analysis
PATH_PUBLIC_REFERENCE_DATASET = file.path( PATH_EXPERIMENT_REFERENCE, "public_reference_dataset", "regional_cell_sampling_Log2TPM_fixed.txt")

# Path to the RNA-seq gene expression matices from project experiment
PATH_RNASEQ_EXPRESSION_MATRIX_CHIP1 = file.path( PATH_EXPERIMENT_RAWDATA, 
                                                 "Gene_expression_matrixes", 
                                                 "Chip1", 
                                                 "Auto_Zanoni_LP1_093021_Trx_MOUSE_s5-torrent-server-vm_247.rpm.bcmatrix.csv")

PATH_RNASEQ_EXPRESSION_MATRIX_CHIP2 = file.path( PATH_EXPERIMENT_RAWDATA, 
                                                 "Gene_expression_matrixes", 
                                                 "Chip2", 
                                                 "Auto_Zanoni_LP2_093021_Trx_MOUSE_s5-torrent-server-vm_249.rpm.bcmatrix.csv")

# Path tot the RNA-seq metadata file
PATH_RNASEQ_METADATA = file.path(  PATH_EXPERIMENT_RAWDATA, 
                                   "Gene_expression_matrixes",
                                   "001_Iontorrent_RNAseq_211006_Metadata.csv")
