# Merge clusters from TSV file (taken from 01_QC) to be used in downstream analyses

library( funr)
library(forcats);


# Get path of current script
WORKING_DIR = dirname( sys.script())

# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original clustering output and read as data.frame
inputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                      "02_GlobalHeterogeneity",
                      "ClusteringRes_1",
                      "Extra",
                      "002_IFNL_Recovery_scRNAseq_20230601_10X_RNAseq_KO_cellsClusterIdentity.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "Dying enterocytes mito Hi" = "8", 
                                   "Dying enterocytes mito Hi" = "3", 
                                   "Dying enterocytes mito Hi" = "4", 
                                   "Mature distal enterocytes" = "0", 
                                   "Mature distal enterocytes" = "1", 
                                   "Mature distal enterocytes" = "2", 
                                   "Mature distal enterocytes" = "5", 
                                   "Proximal enterocytes" = "6",
                                   "Goblet" = "10",
                                   "Goblet" = "7",
                                   "Mix Progenitors/Stem/Tuft" = "11",
                                   "Cycling+Stem" = "12",
                                   "Enterocytes progenitors" = "9",
                                   "Paneth+Stem" = "13",
                                   "Enteroendocrine" = "14");

# Write results
outputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                       "14a_GroupClusters_AllCells_From02",
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




