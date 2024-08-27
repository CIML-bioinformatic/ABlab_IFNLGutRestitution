# Merge clusters from TSV file (taken from 01_QC) to be used in downstream analyses

library( funr)
library(forcats);


# Get path of current script
WORKING_DIR = dirname( sys.script())

# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original clustering output and read as data.frame
inputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                      "02_GlobalHeterogeneity_reClustering",
                      "ClusteringRes_1",
                      "Extra",
                      "002_IFNL_Recovery_scRNAseq_20240418_WT_KO_MergeAllCells_cellsClusterIdentity.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "Goblet" = "8", 
                                   "Goblet" = "10", 
                                   "Goblet" = "11", 
                                   "Enterocyte progenitor cycling" = "13",
                                   "Enterocyte progenitors II" = "7",
                                   "Enterocyte (proximal)" = "2",
                                   "Enterocyte (proximal)" = "3",
                                   "Enterocyte (proximal)" = "5", # "Enterocyte (proximal damaged)" = "5",
                                   "Enterocyte (distal)" = "1",
                                   "Enterocyte (distal)" = "4",
                                   "Enterocyte (distal)" = "0", # "Enterocyte (distal damaged)" = "0",
                                   "Stem/Progenitor" = "6",
                                   "Stem/Progenitor" = "12",
                                   "Stem/Progenitor" = "15",
                                   "Stem/Progenitor" = "17",
                                   "Paneth" = "14",
                                   "Enteroendocrine" = "16",
                                   "Tuft" = "18",
                                   "Dead" = "9");

# Write results
outputTSV = file.path( PATH_EXPERIMENT_OUTPUT,
                       "03_GroupClusters",
                       "002_IFNL_Recovery_scRNAseq_20240418_WT_KO_MergeAllCells_cellsClusterIdentity.tsv")

# Make sure output directory exists
dir.create( dirname( outputTSV))

# Save cells cluster identity as determined with 'FindClusters'
write.table( idents, 
             file = outputTSV, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




