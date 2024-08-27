# ##################################################
# Load the data and prepare the dataframe
# ##################################################

## @knitr prepare_data

# Load the data on gene expression
# ________________________________
expression_matrix_df = as.data.frame( read_excel( PATH_EXPRESSION_MATRIX_FILE, sheet = "countData"))
gene_data_df = read_excel( PATH_EXPRESSION_MATRIX_FILE, sheet = "rowData")
row.names( expression_matrix_df) = gene_data_df$Gene_name

# Load the metadata on samples
# ______________________________
sample_metadata_df = as.data.frame( read_excel( PATH_EXPRESSION_MATRIX_FILE, sheet = "colData"))
row.names( sample_metadata_df) = sample_metadata_df$sample.id

# Select the samples by group of interest
CONTROL_SAMPLES = sample_metadata_df$sample.id[ which( sample_metadata_df$condition.name %in% CONDITION_ALL_CONTROL)]
QUIESCENT_SAMPLES = sample_metadata_df$sample.id[ which( sample_metadata_df$condition.name %in% CONDITION_ALL_QUIESCENT)]
ACTIVE_SAMPLES = sample_metadata_df$sample.id[ which( sample_metadata_df$condition.name %in% CONDITION_ALL_ACTIVE)]

# Add a global group information called "status"
sample_metadata_df$status = NA
sample_metadata_df[ CONTROL_SAMPLES, "status"] = STATUS_CONTROL
sample_metadata_df[ QUIESCENT_SAMPLES, "status"] = STATUS_QUIESCENT
sample_metadata_df[ ACTIVE_SAMPLES, "status"] = STATUS_ACTIVE
sample_metadata_df$status = factor( sample_metadata_df$status, levels = c( STATUS_CONTROL, STATUS_QUIESCENT, STATUS_ACTIVE))

# Show the sample metatdata table on report
sample_metadata_df[ , c("sample.name", "condition.name", "status")] %>%
  kbl() %>%
  kable_styling()


# Load the hallmarks definitions
# ______________________________

HALLMARKS_LIST = list()
hallmark_con = file( PATH_HALLMARKS_FILE, "r")
while ( TRUE ) {
  line = readLines( hallmark_con, n = 1)
  if ( length( line) == 0 ) {
    break
  }
  tokens = strsplit( line, "\t")[[1]]
  HALLMARKS_LIST[[ tokens[ 1]]] = tokens[3:length( tokens)]
}

close( hallmark_con)

