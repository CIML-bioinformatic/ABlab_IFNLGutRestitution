# ##################################################
# Load the expression data and experiments meta-data
# ##################################################

## @knitr prepareData

# Load the gene expression matrices
expression_data_chip1_df = read.table( PATH_RNASEQ_EXPRESSION_MATRIX_CHIP1, header= TRUE, sep = "\t", row.names = 1)
expression_data_chip2_df = read.table( PATH_RNASEQ_EXPRESSION_MATRIX_CHIP2, header= TRUE, sep = "\t", row.names = 1)

# Load the meta-data of the experiments
metadata_df = read.table( PATH_RNASEQ_METADATA, header = TRUE, sep ="\t")

row.names( metadata_df) = metadata_df$ExperimentID
metadata_df$Condition = as.factor( metadata_df$Condition)

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
  tokens = tokens[ which( nchar( tokens) > 0)]
  HALLMARKS_LIST[[ tokens[ 1]]] = tokens[3:length( tokens)]
}

close( hallmark_con)

# If file containing the list of hallmarks with MOUSE genes does not exists, create it by converting the HUMAN lists
if( !file.exists( PATH_HALLMARKS_MOUSE_FILE)){
  mouse_hallmarks_list = convertHumanGeneAllList( HALLMARKS_LIST)
  dir.create( dirname( PATH_HALLMARKS_MOUSE_FILE), recursive = TRUE, showWarnings = FALSE)
  sink( PATH_HALLMARKS_MOUSE_FILE)
  for( hallmark in names( mouse_hallmarks_list)){
    cat( hallmark, "\t", paste( mouse_hallmarks_list[[ hallmark]], collapse="\t"), "\n", sep="")
  }
  sink()
}

# Load the list of hallmarks with MOUSE genes 
HALLMARKS_LIST_MOUSE = list()
hallmark_mouse_con = file( PATH_HALLMARKS_MOUSE_FILE, "r")
while ( TRUE ) {
  line = readLines( hallmark_mouse_con, n = 1)
  if ( length( line) == 0 ) {
    break
  }
  tokens = strsplit( line, "\t")[[1]]
  tokens = tokens[ which( nchar( tokens) > 0)]
  HALLMARKS_LIST_MOUSE[[ tokens[ 1]]] = tokens[3:length( tokens)]
}

close( hallmark_mouse_con)

# Load the list of STAPP pathways
STAPP_PATHWAYS_LIST = list()
stapp_pathways_con = file( PATH_STAPP_PATHWAYS_FILE, "r")
while ( TRUE ) {
  line = readLines( stapp_pathways_con, n = 1)
  if ( length( line) == 0 ) {
    break
  }
  tokens = strsplit( line, "\t")[[1]]
  tokens = tokens[ which( nchar( tokens) > 0)]
  STAPP_PATHWAYS_LIST[[ tokens[ 1]]] = tokens[3:length( tokens)]
}

close( stapp_pathways_con)


# Load the Gene Ontology data
# ______________________________

GOBP_CELL_POPULATION_PROLIFERATION_geneset = unique( read.table( PATH_GOBP_CELL_POPULATION_PROLIFERATION_FILE)$V1)



