# ########################################################################
# This script aims to load the data from public reference dataset
# and convert it to become the reference dataset for the CibersortX
# analysis
# ########################################################################

## @knitr loadReferenceData

# Load the reference data from file
original_reference_df = read.table( PATH_PUBLIC_REFERENCE_DATASET, header= TRUE, sep="\t")

cat("<BR>Reference data loaded from", PATH_PUBLIC_REFERENCE_DATASET)
cat("<BR>Number of cells=", ncol( original_reference_df))
cat("<BR>Number of genes=", nrow( original_reference_df))

# Push the column GENE as row names and remove the column for easy math computations
row.names( original_reference_df) = original_reference_df$GENE
original_reference_df$GENE = NULL

# Convert log2(TPM+1) to linear TPM
cat("<BR><BR>Converting log2TPM data to TPM")
original_reference_df = 2^original_reference_df - 1

# Remove from cell name the cell ID and mouse ID to keep only cell type and tissue
cat("<BR>Changing names of cells to recover only cell type")
bc_cell_type = gsub( "_m\\d{1}_", "_", gsub( "[ATCG]{14}_", "", names( original_reference_df)))
names( original_reference_df) = bc_cell_type

# Display a table summarizing the Tissue/cell type frequency
DT::datatable( as.data.frame.table( table( names( original_reference_df))), colnames = c( "Tissue/Cell type", "Frequency"))


for( REFERENCE_FILTER in REFERENCE_FILTER_LIST){

  cat("<H4>Select data for reference tissue", REFERENCE_FILTER, "</H4>")
  
  # Restrict the data to a specific tissue or take all the tissues
  if( REFERENCE_FILTER != TISSUE_ALL){
    # Keep only cell type/tissue  of a specific tissue
    selected_samples = grepl( REFERENCE_FILTER, names( original_reference_df))
    selected_samples_names = names( original_reference_df)[ selected_samples]
    filtered_reference_df = original_reference_df[ ,  selected_samples]
  }else{
    # Remove the tissue of origin from the cell type/tissue name so to get merging of cell type for all tissues
    selected_samples_names = names( original_reference_df)
    for( filter in REFERENCE_FILTER_LIST){
      selected_samples_names = gsub( paste0( filter, "_"), "", selected_samples_names)
    }
    filtered_reference_df = original_reference_df
  }
  
  # Add the Gene symbol in first column
  filtered_reference_df = cbind( GeneSymbol = row.names( original_reference_df), filtered_reference_df)
  names( filtered_reference_df) = c( "GeneSymbol", selected_samples_names)
  
  # Export the table to file
  export_reference_file_path = file.path( PATH_ANALYSIS_OUTPUT, 
                                          paste0( 
                                            tools::file_path_sans_ext( basename( PATH_PUBLIC_REFERENCE_DATASET)), 
                                            "_nolog2_nocellid_", REFERENCE_FILTER, ".txt"))
  
  cat("<BR><BR>Export reference data to file:", export_reference_file_path)
  write.table( filtered_reference_df,
               file = export_reference_file_path,
               col.names = TRUE,
               row.names = FALSE,
               sep="\t",
               quote=FALSE)

}
