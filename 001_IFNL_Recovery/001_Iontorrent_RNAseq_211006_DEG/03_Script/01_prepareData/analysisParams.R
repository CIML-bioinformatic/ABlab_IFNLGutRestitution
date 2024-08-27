###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "01_prepareData"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Prepare the data for DEG analysis"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the RNA-seq gene expression matices from project experiment
PATH_RNASEQ_EXPRESSION_MATRIX_CHIP1 = file.path( PATH_EXPERIMENT_RAWDATA, 
                                                 "Gene_expression_matrixes", 
                                                 "Chip1", 
                                                 "Auto_Zanoni_LP1_093021_Trx_MOUSE_s5-torrent-server-vm_247.bcmatrix.tsv")

PATH_RNASEQ_EXPRESSION_MATRIX_CHIP2 = file.path( PATH_EXPERIMENT_RAWDATA, 
                                                 "Gene_expression_matrixes", 
                                                 "Chip2", 
                                                 "Auto_Zanoni_LP2_093021_Trx_MOUSE_s5-torrent-server-vm_249.bcmatrix.tsv")

# Path tot the RNA-seq metadata file
PATH_RNASEQ_METADATA = file.path(  PATH_EXPERIMENT_RAWDATA, 
                                   "Gene_expression_matrixes",
                                   "001_Iontorrent_RNAseq_211006_Metadata.csv")
