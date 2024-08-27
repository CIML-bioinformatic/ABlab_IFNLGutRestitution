###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "01_IFNL_vs_CDUC"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analyse relationship between IFNL and CD/UC disease states"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the RNA-seq gene expression matrices from the CD/UC project experiment of PNlab
PATH_EXPRESSION_MATRIX_FILE = file.path( PATH_EXPERIMENT_RAWDATA, "expression_data_ibd_log2_vst_AB.xlsx")

# Class the various disease conditions:
CONDITION_ALL_CONTROL = c( "CTRL")
CONDITION_ALL_QUIESCENT = c( "QNT", "QTNF", "Q5ASA")
CONDITION_ALL_ACTIVE = c("AMNT", "AMTNF", "AM5ASA", "ASCort")

# Definition of status
STATUS_CONTROL = "Control"
STATUS_QUIESCENT = "Quiescent"
STATUS_ACTIVE = "Active"

# Path to the file with HALLMARK definitions
PATH_HALLMARKS_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "GSEA_hallmarkGeneSets", "h.all.v7.4.symbols.gmt")

# Define the Hallmark to analyse
HALLMARKS_TO_ANALYSE = c( "HALLMARK_INTERFERON_ALPHA_RESPONSE")
