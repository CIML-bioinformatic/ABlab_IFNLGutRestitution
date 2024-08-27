# ##################################################
# Load the expression data and experiments meta-data
# ##################################################

## @knitr executeDEGAnalysis


cat("<H3>", TEST_CONDITION, "versus", REFERENCE_CONDITION, "</H3>")

# select the experiment corresponding to the reference and test conditions
reference_experiment_set = metadata_df[ which( metadata_df$Condition == REFERENCE_CONDITION), "ExperimentID"]
test_experiment_set = metadata_df[ which( metadata_df$Condition == TEST_CONDITION), "ExperimentID"]

# group the test and reference experiments
all_experiment_set = c( reference_experiment_set, test_experiment_set)

# Choose the right chip of gene expression according to the selected experiments
if( sum( sapply( all_experiment_set, function( x){ return( x %in% names( expression_data_chip1_df))})) == length( all_experiment_set)){
  cat("<BR>Using Chip1 data")
  gene_expression_df = expression_data_chip1_df[ , all_experiment_set]
}else if( sum( sapply( all_experiment_set, function( x){ return( x %in% names( expression_data_chip2_df))})) == length( all_experiment_set)){
  cat("<BR>Using Chip2 data")
  gene_expression_df = expression_data_chip2_df[ , all_experiment_set]
}else{
  stop( "Error: analysis across chip asked.")
}

# Select the meta-data of the chosen experiments
coldata = metadata_df[ all_experiment_set, ]
coldata$Condition = factor( coldata$Condition, levels = c( REFERENCE_CONDITION, TEST_CONDITION))

# Build the Deseq2 object
dds <- DESeqDataSetFromMatrix(countData = gene_expression_df,
                              colData = coldata,
                              design = ~ Condition)

# Filter genes with low count across conditions
keep <- rowSums( counts( dds)) >= FILTER_MINIMUM_ROWSUM_COUNT
dds <- dds[ keep, ]

# Compute the Differential expressed genes
dds = DESeq(dds)

# Extract the result for the DEG analysis using log shrinkage and reordering by adjusted pvalue
# res = results(dds, contrast=c( "Condition", TEST_CONDITION, REFERENCE_CONDITION))
resLFC <- lfcShrink(dds, coef=paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION), type="apeglm")
resLFC <- resLFC[ order( resLFC$padj),]

# Export complete result to file
write.csv( as.data.frame( resLFC), 
           file=file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, ".csv")))

# Export top200 most DE result to file
write.csv( as.data.frame( resLFC)[ 1:200, ], 
           file=file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_Top200.csv")))


# Study regularized expression
# _________________________________________

cat("<H4>Study of regularized expressions</H4>")

# Compute log2(n + 1) expression
ntd <- normTransform(dds)
# Compute variance stabilized transformed (VST) expressions
vsd <- vst( dds, blind=FALSE)
# Compute log transformed (rld) expression
rld <- rlog( dds, blind=FALSE)

pairs(log10( 1+ gene_expression_df), cex = 1, pch = 1, cex.labels = 1, font.labels = 1, panel = panel.smooth,
          bg = "light blue", horOdd=TRUE)

# Plot the data in the various regularization
meanSdPlot(assay(ntd), plot = FALSE)$gg + ggtitle( "Mean versus SD plot of log transformed expression")
meanSdPlot(assay(vsd), plot = FALSE)$gg + ggtitle( "Mean versus SD plot of variance stabilized transformed (VST) expression")
meanSdPlot(assay(rld), plot = FALSE)$gg + ggtitle( "Mean versus SD plot of regularized log transformed (rld) expression")

# Plot the PCA of the various regularization
plotPCA( ntd, intgroup=c("Condition")) + ggtitle( "PCA of log transformed expression")
plotPCA( vsd, intgroup=c("Condition")) + ggtitle( "PCA of variance stabilized transformed (VST) expression")
plotPCA( rld, intgroup=c("Condition")) + ggtitle( "PCA of regularized log transformed (rld) expression")

# Heatmap of the sample-to-sample distances
# _________________________________________

cat("<H4>Heatmap of the sample-to-sample distances (VST expression)</H4>")

sampleDists <- dist( t( assay( vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Study the result of DEG
# _________________________________________

cat("<H4>Study of DEG result</H4>")

# Look at the MA plot
plotMA( resLFC, ylim=c( min(resLFC$log2FoldChange, na.rm = TRUE), max( resLFC$log2FoldChange, na.rm = TRUE)))

# Look at the expression levels of the 9 most DEG
# best_genes = row.names( resLFC)[1:9]
# best_genes_plot = list()
# 
# for( gene_index in 1:length( best_genes)){
#   gene_data = plotCounts(dds, gene=which( row.names( resLFC) == best_genes[ gene_index]), intgroup="Condition", returnData=TRUE)
#   best_genes_plot[[ gene_index]] = ggplot( gene_data, aes(x=Condition, y= count)) + 
#                                   geom_point(position=position_jitter(w=0.1,h=0)) + 
#                                   scale_y_log10() + ggtitle( best_genes[ gene_index]) +
#                                   theme_minimal() 
# }
# 
# grid.arrange( best_genes_plot[[1]], best_genes_plot[[2]], best_genes_plot[[3]],
#               best_genes_plot[[4]], best_genes_plot[[5]], best_genes_plot[[6]],
#               best_genes_plot[[7]], best_genes_plot[[8]], best_genes_plot[[9]], ncol =3)


EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0( "Condition ", TEST_CONDITION, " vs ", REFERENCE_CONDITION),
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 2.0,
                labSize = 6.0,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 2.0,
                drawConnectors = TRUE)

DT::datatable( as.data.frame( resLFC[ 1:200,]), caption = "Differential expression data of the 200 top DEG")
