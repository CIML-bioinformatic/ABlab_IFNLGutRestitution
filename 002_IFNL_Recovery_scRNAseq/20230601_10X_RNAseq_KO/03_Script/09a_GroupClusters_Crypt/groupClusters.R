# Merge clusters from TSV file (taken from 01_QC) to be used in downstream analyses

library( funr)
library(forcats);


# Get path of current script
WORKING_DIR = dirname( sys.script())

# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original clustering output and read as data.frame
inputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                      "05a_GlobalHeterogeneity_Crypt",
                      "ClusteringRes_1",
                      "Extra",
                      "002_IFNL_Recovery_scRNAseq_20230601_10X_RNAseq_KO_cellsClusterIdentity.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "Goblet" = "0", 
                                   "Goblet" = "1", 
                                   "Goblet" = "6", 
                                   "Enterocytes" = "3",
                                   "Transient Ampli Progenitor" = "5",
                                   "Cycling CBC (stem)" = "2",
                                   "Revival Stem" = "7",
                                   "Paneth" = "8",
                                   "Enteroendocrine" = "9",
                                   "Tuft" = "11",
                                   "Dead" = "4",
                                   "Paneth/Stem" = "10");

# Write results
outputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                       "09a_GroupClusters_Crypt",
                       "002_IFNL_Recovery_scRNAseq_20230601_10X_RNAseq_KO_cellsClusterIdentity.tsv")

# Make sure output directory exists
dir.create( dirname( outputTSV))

# Save cells cluster identity as determined with 'FindClusters'
write.table( idents, 
             file = outputTSV, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




