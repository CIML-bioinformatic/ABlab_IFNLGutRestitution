# Load seurat object and subset it based on cell identity as specified in tsv file
# Used for downstream analysis when clustering is not taken as resulting from 
# subclustering script (which gives a subset Seurat object too).


library( funr)
library( forcats);
library( Seurat);


# Get path of current script
WORKING_DIR = dirname( sys.script())
# debug: WORKING_DIR = getwd()


# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));
source( file.path( WORKING_DIR, "..", "sampleParams.R"));

# Get path to original seurat object (take it from subclustering results)
inputSeurat = dir( file.path( PATH_EXPERIMENT_OUTPUT,
                              "02_GlobalHeterogeneity",
                              "ClusteringRes_0.4"),
                   pattern = "_seuratObject_final.RDS",
                   full.names = TRUE)

# Load Seurat object
sc10x = readRDS(inputSeurat)


# Get path to TSV files identifying cells of each population (manual selection with interactive umap)
inputFilesTSV = dir( file.path( PATH_EXPERIMENT_REFERENCE,
                                "06_CellsSelections",
                                "01_Crypt_VS_Enterocytes"),
                      pattern = SAMPLE_NAME,
                      full.names = TRUE)

names(inputFilesTSV) = gsub(".csv", "", basename(inputFilesTSV))

cellsGroups = lapply(inputFilesTSV, 
                     read.csv, 
                     row.names = 1,
                     header = TRUE);




###
# Subset seurat object and save result

# Make sure output directory exists
outputFolder = file.path( PATH_EXPERIMENT_OUTPUT,
                          "04_subsetSeurat_Crypt_Enterocytes")
dir.create( outputFolder)

for(currentSetName in names(cellsGroups))
{
  message(currentSetName)
  # Subset the original object using cell names
  resultSubset = sc10x[, rownames(cellsGroups[[currentSetName]])]
  # Update identity (in case downstream don't want to load identity csv file separately)
  Idents(resultSubset) = cellsGroups[[currentSetName]][["Identity"]]
  
  # Save as binary file for downstream analyses
  saveRDS( object = resultSubset,
           file = file.path(outputFolder, paste0("seuratObject_subset_", currentSetName, ".RDS")))
  
}



