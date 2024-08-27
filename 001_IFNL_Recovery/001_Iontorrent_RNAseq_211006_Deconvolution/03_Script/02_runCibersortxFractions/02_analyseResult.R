# ##############################################################################
# This script aims to analyse the results of the Cibersortx Franction analysis
# ##############################################################################

## @knitr analyseResult

# Load meta-data
# rnaseq_metadata = read.table( PATH_RNASEQ_METADATA, sep="\t", header=TRUE)

# __________________________________________________
# Load Cybersort results
# __________________________________________________

cat("<BR><b>Loading Cibersort results from:</b>")
cat("<BR>* ", PATH_CIBERSORT_ANALYSIS_RESULT_CHIP1)
cat("<BR>* ", PATH_CIBERSORT_ANALYSIS_RESULT_CHIP2)
cibersort_chip1 = read.table( PATH_CIBERSORT_ANALYSIS_RESULT_CHIP1, sep="\t", header=TRUE)
cibersort_chip2 = read.table( PATH_CIBERSORT_ANALYSIS_RESULT_CHIP2, sep="\t", header=TRUE)

# __________________________________________________
# Look at the signature of cell type
# __________________________________________________
cat("<H4>Analysis of the cell type signature</H4>")
cat("<BR>Here we present the cell type signature infered from the single-cell data. We first show the heatmap (normalized by row)")
cat("of the signature. Then for each cell type, we look in a table at the most caracteristic genes.")
cat("<BR><BR>")

# Load the signature data
cibersort_signature = read.table( PATH_CIBERSORT_SIGNATURE_RESULT, sep="\t", header=TRUE)
row.names( cibersort_signature) = cibersort_signature$GeneSymbol
cibersort_signature$GeneSymbol = NULL

# Remove the lines with only zeros
row_sums = apply( cibersort_signature, 1, sum)
cibersort_signature = cibersort_signature[ !(row_sums==0), ]

# Normalize the gene expression by cell type
cibersort_signature_rownorm = apply( cibersort_signature, 1, function(x){ return( (x-mean(x))/(sd(x)))})

# Show the heatmap
pheatmap::pheatmap( cibersort_signature_rownorm)

# For each cell type, get the most expressed genes
signature_cell_type_df = data.frame()
for( cell_type in row.names( cibersort_signature_rownorm)){
  cell_type_df = cibersort_signature[ which( cibersort_signature_rownorm[ cell_type, ] >= 0.95*max( cibersort_signature_rownorm[ cell_type, ])),]
  cell_type_df = cbind( CellType=cell_type, cell_type_df)
  signature_cell_type_df = rbind( signature_cell_type_df, cell_type_df)
}

DT::datatable( signature_cell_type_df)

# __________________________________________________
# Look at quality of results
# __________________________________________________
cat("<H4>Statistical analysis of results of cell type inference</H4>")

cat("<BR>Here we present the statistical information on th cell type inference to check for the inference quality")
cat("<BR><BR>")
DT::datatable( rbind( cibersort_chip1[, c( "Mixture", "P.value", "Correlation", "RMSE")],
                      cibersort_chip2[, c( "Mixture", "P.value", "Correlation", "RMSE")]),
               caption = "Statistical result of Cibersort run"
             )

# __________________________________________________
# Look at the cell type proportions inference
# __________________________________________________

# Convert cell type proportions to long format
cibersort_chip1_long <- cibersort_chip1[, 1:(ncol( cibersort_chip1)-3)] %>% gather(CellType, Value, -c(Mixture))
cibersort_chip2_long <- cibersort_chip2[, 1:(ncol( cibersort_chip2)-3)] %>% gather(CellType, Value, -c(Mixture))

# Pool the data of cell type proportions
cibersort_chip12_long = rbind( cibersort_chip1_long, cibersort_chip2_long)

# Plot the results of cell type proportions
cat("<H4>Result of the cell type inference</H4>")
cat("<BR>Here we present the result of the cell type inference, in a cumulative barplot and in a numerical datatable.")
cat("<BR><BR>")
ggplot( cibersort_chip12_long) + 
  geom_bar( aes( x=Mixture, y=Value, fill=CellType), stat="identity") +
  scale_fill_manual( values=RColorBrewer::brewer.pal(n = length( unique( cibersort_chip12_long$CellType)), name = "Set3")) +
  theme_classic() + theme( axis.text.x = element_text( angle=45, h=1)) +
  ggtitle( "Cell type inference result")

ggplot( cibersort_chip12_long) + 
  geom_bar( aes( x=Mixture, y=Value, fill=CellType), stat="identity") +
  facet_wrap( .~CellType) +
  scale_fill_manual( values=RColorBrewer::brewer.pal(n = length( unique( cibersort_chip12_long$CellType)), name = "Set3")) +
  theme_classic() + theme( axis.text.x = element_text( angle=45, h=1)) +
  ggtitle( "Cell type inference result")

DT::datatable( rbind( cibersort_chip1[, 1:(ncol( cibersort_chip1)-3)],
                      cibersort_chip2[, 1:(ncol( cibersort_chip1)-3)]),
               caption = "Cell type inference result"
)


