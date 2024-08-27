###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#
# Parameters for individual samples (e.g. SAMPLE_NAME, PATH_DATA) must be set
# in the file 'sampleParams.R' (stored in the sample folder). Any value defined 
# there will override values defined in this file.
#


    

#### General

PROJECT_PARENT_NAME = "IFNLGutRestitution"
PROJECT_NAME = "ABlab - IFNL Recovery - Data of 10/06/21";
PROJECT_SHORT_NAME = "IFNL_Recovery"

PATH_PROJECT = "/mnt/DOSI/ABLAB/BIOINFO/Project/IFNLGutRestitution/001_IFNL_Recovery"

EXPERIMENT_NAME = "001_Iontorrent_RNAseq_211006_DEG"

#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_EXPERIMENT_RAWDATA = file.path( PATH_PROJECT, EXPERIMENT_NAME, "00_RawData")
PATH_EXPERIMENT_REFERENCE = file.path( PATH_PROJECT, EXPERIMENT_NAME, "01_Reference")
PATH_EXPERIMENT_OUTPUT = file.path( PATH_PROJECT, EXPERIMENT_NAME, "05_Output")

#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



