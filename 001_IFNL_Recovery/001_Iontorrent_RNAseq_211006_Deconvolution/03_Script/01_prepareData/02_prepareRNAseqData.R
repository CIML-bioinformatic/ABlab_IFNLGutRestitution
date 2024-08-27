# ########################################################################
# This script aims to load the data from public reference dataset
# and convert it to become the reference dataset for the CibersortX
# analysis
# ########################################################################

## @knitr loadRNAseqData

# Read the RNA-seq gene expression data from the project experiment
original_rnaseq_chip1 <- read.table( PATH_RNASEQ_EXPRESSION_MATRIX_CHIP1, header = TRUE, sep = ",")
original_rnaseq_chip2 <- read.table( PATH_RNASEQ_EXPRESSION_MATRIX_CHIP2, header = TRUE, sep = ",")

cat("<BR>Gene expression data loaded from", PATH_RNASEQ_EXPRESSION_MATRIX_CHIP1)
cat("<BR>Number of RNA-seq=", ncol( original_rnaseq_chip1))
cat("<BR>Number of genes=", nrow( original_rnaseq_chip1))

cat("<BR><BR>Gene expression data loaded from", PATH_RNASEQ_EXPRESSION_MATRIX_CHIP2)
cat("<BR>Number of RNA-seq=", ncol( original_rnaseq_chip2))
cat("<BR>Number of genes=", nrow( original_rnaseq_chip2))

# Remove the Target column from the data
if( "Target" %in% names( original_rnaseq_chip1)){
  original_rnaseq_chip1$Target = NULL
}
if( "Target" %in% names( original_rnaseq_chip2)){
  original_rnaseq_chip2$Target = NULL
}

# Load the meta-data of the RNA-seq experiment
cat("<BR>Experiment meta-data loaded from", PATH_RNASEQ_METADATA)
original_rnaseq_metadata = read.table( PATH_RNASEQ_METADATA, header = TRUE, sep = ",", skip = 1)

# Retrieve the Ioncode from the run name
original_rnaseq_metadata$IonCode = stringr::str_match( original_rnaseq_metadata$run_name, "\\IonCode_\\d{4}")

# Define a string (ExperimentID) to summarize the experiment conditions
original_rnaseq_metadata$ExperimentID = paste0( original_rnaseq_metadata$genotype, "_",
                                              original_rnaseq_metadata$treatment, "_", 
                                              original_rnaseq_metadata$sample_ID)

# Plot the datatable of meta-data
DT::datatable( original_rnaseq_metadata)

# Export meta-data to file
export_metadata_file_path = file.path( PATH_ANALYSIS_OUTPUT, 
                                       paste0( 
                                          tools::file_path_sans_ext( basename( PATH_RNASEQ_METADATA)), 
                                          ".txt"))
cat("<BR>Export RNA-seq meta-data to file:", export_metadata_file_path)

write.table( original_rnaseq_metadata,
             file = export_metadata_file_path,
             col.names = TRUE,
             row.names = FALSE,
             sep="\t",
             quote=FALSE)



# Assign the ExperimentID to the column names of the gene expression data according to the matching IonCode
names( original_rnaseq_chip1) = sapply( names( original_rnaseq_chip1), function( name){
  ion_code_index = which( original_rnaseq_metadata$IonCode == as.character( name))
  if( length( ion_code_index) > 0 && ion_code_index > 0){
    return( original_rnaseq_metadata[ ion_code_index, "ExperimentID"])
  }else{
    return( name)
  }
})

names( original_rnaseq_chip2) = sapply( names( original_rnaseq_chip2), function( name){
  ion_code_index = which( original_rnaseq_metadata$IonCode == as.character( name))
  if( length( ion_code_index) > 0 && ion_code_index > 0){
    return( original_rnaseq_metadata[ ion_code_index, "ExperimentID"])
  }else{
    return( name)
  }
})

# Export the modified gene expression data to file
export_rnaseq_chip1_file_path = file.path( PATH_ANALYSIS_OUTPUT, 
                                           paste0( 
                                             tools::file_path_sans_ext( basename( PATH_RNASEQ_EXPRESSION_MATRIX_CHIP1)), 
                                             ".txt"))
cat("<BR>Export RNA-seq chip1 data to file:", export_rnaseq_chip1_file_path)

write.table( original_rnaseq_chip1,
             file = export_rnaseq_chip1_file_path,
             col.names = TRUE,
             row.names = FALSE,
             sep="\t",
             quote=FALSE)


export_rnaseq_chip2_file_path = file.path( PATH_ANALYSIS_OUTPUT, 
                                          paste0( 
                                            tools::file_path_sans_ext( basename( PATH_RNASEQ_EXPRESSION_MATRIX_CHIP2)), 
                                            ".txt"))
cat("<BR>Export RNA-seq chip2 data to file:", export_rnaseq_chip2_file_path)

write.table( original_rnaseq_chip2,
             file = export_rnaseq_chip2_file_path,
             col.names = TRUE,
             row.names = FALSE,
             sep="\t",
             quote=FALSE)
