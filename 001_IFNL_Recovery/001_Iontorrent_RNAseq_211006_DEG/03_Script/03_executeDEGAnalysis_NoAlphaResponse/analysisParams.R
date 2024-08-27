###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "03_executeDEGAnalysis_NoAlphaResponse"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Produce the DEG analysis over several pairs of conditions, removing genes from the HALLMARK Interferon alpha response"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the RNA-seq gene expression matices from project experiment
PATH_RNASEQ_EXPRESSION_MATRIX_CHIP1 = file.path( PATH_EXPERIMENT_OUTPUT,
                                                 "01_prepareData",
                                                 "Auto_Zanoni_LP1_093021_Trx_MOUSE_s5-torrent-server-vm_247.bcmatrix.txt")

PATH_RNASEQ_EXPRESSION_MATRIX_CHIP2 = file.path( PATH_EXPERIMENT_OUTPUT, 
                                                 "01_prepareData",
                                                 "Auto_Zanoni_LP2_093021_Trx_MOUSE_s5-torrent-server-vm_249.bcmatrix.txt")

# Path tot the RNA-seq metadata file
PATH_RNASEQ_METADATA = file.path(  PATH_EXPERIMENT_OUTPUT, 
                                   "01_prepareData",
                                   "001_Iontorrent_RNAseq_211006_Metadata.txt")

# Path to Hallmark definition file
PATH_HALLMARKS_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "GSEA_hallmarkGeneSets", "h.all.v7.4.symbols.gmt")

# Path to Hallmark definition file
PATH_HALLMARKS_MOUSE_FILE = file.path( PATH_ANALYSIS_OUTPUT, "GSEA_hallmarkGeneSets", "h.all.v7.4.symbols_mouse.gmt")

# Path to the file with GOBP_CELL_POPULATION_PROLIFERATION
PATH_GOBP_CELL_POPULATION_PROLIFERATION_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "GeneOntology", "GOBP_CELL_POPULATION_PROLIFERATION.txt")

# Path to STAPP Pathway definition file
PATH_STAPP_PATHWAYS_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "STAPP_pathways.gmt")

# List of HALLMARK to remove from the expression matrix
HALLMARKS_TO_REMOVE = c( "HALLMARK_INTERFERON_ALPHA_RESPONSE")

# Minimal number of read counts across samples for genes to be kept for analysis
FILTER_MINIMUM_ROWSUM_COUNT = 10

# First species error to select the DEG
DEG_PADJ_THRESHOLD = 0.1
