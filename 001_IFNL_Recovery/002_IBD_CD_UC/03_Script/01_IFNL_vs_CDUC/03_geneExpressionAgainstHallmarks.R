# ##################################################
# Analyse the behavior of gene IFNL1 agains
# immunological Hallmarks
# ##################################################

## @knitr gene_expression_against_hallmarks

for( current_hallmark in HALLMARKS_TO_ANALYSE){
  
  cat("<H5>Analysing", current_hallmark, "</H5>")
  
  hallmark_genes = HALLMARKS_LIST[[ current_hallmark]]
  cat("<BR>", length( hallmark_genes), "genes in", current_hallmark)
  hallmark_genes = intersect( hallmark_genes, row.names( expression_matrix_df))
  cat("<BR>", length( hallmark_genes), "genes of the hallmark found with expression")
  
  # -- Create dataframe with meta-data and gene expression
  hallmark_expression_df = sample_metadata_df[ c( CONTROL_SAMPLES, QUIESCENT_SAMPLES, ACTIVE_SAMPLES), c("sample.name", "condition.name", "status", "group")]
  hallmark_expression_df = cbind( hallmark_expression_df, t( expression_matrix_df[ c( "GSDMC", hallmark_genes) , c( CONTROL_SAMPLES, QUIESCENT_SAMPLES, ACTIVE_SAMPLES)]))
  hallmark_expression_df$hallmark_mean = unlist( apply( hallmark_expression_df, 1, function( row){
    current_mean = mean( as.numeric( row[ hallmark_genes]), na.rm=TRUE)
    return( current_mean)
  }))
  
  # -- Write data to file
  write.table( hallmark_expression_df[ , c("sample.name", "condition.name", "status", "group", "hallmark_mean")], file = file.path( PATH_ANALYSIS_OUTPUT, paste0( current_hallmark, "_expression.csv")),
               col.names = NA, row.names = TRUE, quote = FALSE, sep=",")
  
  
  print(
    ggplot( hallmark_expression_df) +
    geom_boxplot( aes( x=status, y=hallmark_mean, fill = status)) +
    geom_jitter( aes( x=status, y=hallmark_mean, color = group), width = 0.2, size = 3) +
    ggtitle( paste( "Mean Expression of ", current_hallmark, "across status (boxplot)")) + 
    theme_minimal()
  )
  
  print( 
  ggplot( hallmark_expression_df) +
    geom_violin( aes( x=status, y=hallmark_mean, fill = status)) +
    geom_jitter( aes( x=status, y=hallmark_mean, color = group), width = 0.2, size = 3) +
    ggtitle( paste( "Mean Expression of ", current_hallmark, "across status (violin)")) + 
    theme_minimal()
  )
  
  print(
    ggplot( hallmark_expression_df) +
      geom_point( aes( x= GSDMC, y=hallmark_mean, col = status, pch=group), size = 5) +
      ggtitle( paste( "Mean Expression of ", current_hallmark, "\nagainst expression of GSDMC")) + 
      theme_minimal()
  )
  
  print(
    ggplot( hallmark_expression_df) +
      geom_point( aes( x= GSDMC, y=hallmark_mean, col = group, pch=status), size = 5) +
      ggtitle( paste( "Mean Expression of ", current_hallmark, "\nagainst expression of GSDMC")) + 
      theme_minimal()
  )
  
  regression_plot = ggplot( hallmark_expression_df) +
      geom_point( aes( x= GSDMC, y=hallmark_mean, col = status), size = 3) +
      geom_smooth( aes( x= GSDMC, y=hallmark_mean, col = status), method = "lm") +
      facet_wrap( .~status, ncol = 3) +
      # ggtitle( paste( "Mean Expression of ", current_hallmark, "\nagainst expression of GSDMC with linear regression by status")) + 
      theme_minimal() +
      ylab( paste( current_hallmark)) +
      theme( axis.text = element_text( size = 8),
             axis.title = element_text( size = 8),
             legend.text = element_text( size = 8))
  
  print( regression_plot)
  ggsave( plot = regression_plot, 
          filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( current_hallmark, "_vs_GSDMC_plot.pdf")),
          height = unit( 3, "cm"), width = unit( 6, "cm"))
  
  
  control_cor = cor.test( hallmark_expression_df[ CONTROL_SAMPLES, "hallmark_mean"],
            hallmark_expression_df[ CONTROL_SAMPLES, "GSDMC"],
            method = "spearman")
  
  cat("<BR><BR> On CONTROL samples, the correlation of hallmark against GSDMC has a spearman rho of",
      control_cor$estimate, "with a p-value of", control_cor$p.value)
  
  quiescent_cor = cor.test( hallmark_expression_df[ QUIESCENT_SAMPLES, "hallmark_mean"],
                          hallmark_expression_df[ QUIESCENT_SAMPLES, "GSDMC"],
                          method = "spearman")
  
  cat("<BR> On QUIESCENT samples, the correlation of hallmark against GSDMC has a spearman rho of",
      quiescent_cor$estimate, "with a p-value of", quiescent_cor$p.value)
  
 active_cor = cor.test( hallmark_expression_df[ ACTIVE_SAMPLES, "hallmark_mean"],
                            hallmark_expression_df[ ACTIVE_SAMPLES, "GSDMC"],
                            method = "spearman")
  
  cat("<BR> On ACTIVE samples, the correlation of hallmark against GSDMC has a spearman rho of",
      active_cor$estimate, "with a p-value of", active_cor$p.value)
  
}

