###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "05_GenesList_ExpressionFigures"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Gene Expression"


# !!! THIS ANALYSIS STEP NEEDS TO BE STARTED TWICE !!!
# On first time, it does all computation and generates figures as external files
# (using png/pdf and dev.off) which causes the rmarkdown layout to fail (figures
# in wrong tabset).
# On second execution, the RDATA result of previous run is loaded and previously
# generated external figures are not rendered again, but just integrated to
# rmarkdown using computed file path. It also skips time consuming executions
# (compute DEG analyses and enrichments).
# See chunk 'rmd_loadData' in RMD


# Path to Seurat Object previously saved as RMD file
PATH_RDS_SEURAT_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, 
                                    "04_GlobalHeterogeneity_groupedClusters", 
                                    "002_IFNL_Recovery_scRNAseq_20240418_WT_KO_MergeAllCells_seuratObject_final.RDS");

# List of genes to be analyzed / plotted
GENES_LIST = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
                                            "03_MonitoredGenes",
                                            "230928_genesForHeatmapsAndViolin.tsv"),
                                  sep = "\t",
                                  header = TRUE,
                                  stringsAsFactors = FALSE,
                                  row.names = NULL, 
                                  fill = TRUE));

# Remove empty elements 
GENES_LIST = Map('[', GENES_LIST, lapply(GENES_LIST, function(x){ which( nchar( x)>0)})); # Remove empty strings

# Factor levels used for ordering (specially in heatmaps)
#HTO_FACTOR_LEVELS = c("KO1", "KO2", "KO3", "WT1", "WT2", "WT3")#, "KONeg", "WTNeg"); # Order of levels to use for HTOs (for plots, must match with actual values).
HTO_FACTOR_LEVELS = c("WT", "KO") # Need to specify like this as they get truncated on loading to group individuals together (circumventing the problem with suboptimal HTOs)
HEATMAPS_CLUSTERS_ORDERING = c( "Enterocyte (distal)", "Enterocyte (proximal)", "Enterocyte (distal damaged)", "Goblet", "Stem/Progenitor", "Enterocyte (proximal damaged)", "Enterocyte progenitors II", "Dead", "Enterocyte progenitor cycling", "Paneth", "Enteroendocrine", "Tuft")

MODULES_CONTROL_SIZE = 100 # Number of genes to sample as control for each group of 'module score'

#### Colors

clustersColor = NULL # Named vector or NULL to set ggplot defaults
HTOsColor = NULL # Named vector or NULL to set ggplot defaults



#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (using 'future' for Seurat3, mclapply for
# other loops)
NBCORES = 4;

# Set the max global amount of "shared" memory for paralellization (future)
options(future.globals.maxSize= 891289600)



