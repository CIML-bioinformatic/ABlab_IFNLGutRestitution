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
                      "ClusteringRes_0.8",
                      "Extra",
                      "002_IFNL_Recovery_scRNAseq_20230601_10X_RNAseq_WT_cellsClusterIdentity.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "Enterocytes" = "0", 
                                   "Enterocytes" = "1", 
                                   "Enterocytes" = "9", 
                                   "Dying enterocytes mito hi" = "2", 
                                   "Dying enterocytes mito hi" = "3", 
                                   "Dying enterocytes mito hi" = "4", 
                                   "Dying enterocytes mito hi" = "7", 
                                   "Cycling+Stem" = "11",
                                   "Goblet" = "8",
                                   "Goblet" = "6",                                   
                                   "Paneth" = "13",
                                   "Enterocytes progenitors" = "10",
                                   "Secretory cycling progenitors" = "12",
                                   "Cycling(G2) ribo hi progenitors" = "5");

# Write results
outputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                       "14a_GroupClusters_AllCells_From02",
                       "002_IFNL_Recovery_scRNAseq_20230601_10X_RNAseq_WT_cellsClusterIdentity.tsv")

# Make sure output directory exists
dir.create( dirname( outputTSV))

# Save cells cluster identity as determined with 'FindClusters'
write.table( idents, 
             file = outputTSV, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




