# ##################################################
# Load the expression data and experiments meta-data
# ##################################################

## @knitr functionalAnalysis


# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character( rownames( resLFC))

# Extract significant results
signif_res = list()
signif_res[[ "UPREGULATED" ]] <- resLFC[ resLFC$padj < DEG_PADJ_THRESHOLD & !is.na( resLFC$padj) & resLFC$log2FoldChange >= 1.5, ]
signif_res[[ "DOWNREGULATED" ]] <- resLFC[ resLFC$padj < DEG_PADJ_THRESHOLD & !is.na( resLFC$padj) & resLFC$log2FoldChange <= -1.5, ]

# Enrichment analysis using enrichGO
# __________________________________

for( signif_genes_category in names( signif_res)){

  cat("<H4>GO Biological Processes enrichment using EnrichGO for", signif_genes_category, "genes</H4>")
  
  signif_genes = row.names( signif_res[[ signif_genes_category ]])
  ego <- enrichGO(gene = signif_genes,
                  universe = all_genes,
                  keyType = "SYMBOL",
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = FALSE)
  
  # Output results from GO analysis to a table
  cluster_summary <- data.frame(ego)
  print( DT::datatable( cluster_summary[ , c( "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue")], caption = "EnrichGO result on Biological Processes"))
  write.csv( cluster_summary, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_", signif_genes_category, "_EnrichGO.csv")))
  
  # Show the enrichgo result in dot p lot
  enrichgo_plot = dotplot(ego, showCategory=15, font.size = 8) + 
    theme( legend.text = element_text( size=6),
           legend.title = element_text( size=6),
           axis.text.x = element_text( size = 8),
           axis.text.y = element_text( size = 8),
           axis.title.x = element_text( size = 8),
           axis.title.y = element_text( size = 8))
  print( enrichgo_plot)
  
  ggsave( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_", signif_genes_category, "_EnrichGO.pdf")),
          plot = enrichgo_plot,
          height = 12, width = 15, units = "cm")
  
  # To color genes by log2 fold changes
  signif_res_lFC <- signif_res$log2FoldChange
  print( cnetplot(ego,
           categorySize="pvalue",
           showCategory = 5,
           foldChange= signif_res_lFC,
           vertex.label.font=6)
        )
  
  
  # Use Revigo to reduce the GO annotations
  simMatrix <- calculateSimMatrix(ego$ID,
                                  orgdb="org.Mm.eg.db",
                                  ont="BP",
                                  method="Rel")
  scores <- setNames( -log10( ego$qvalue), ego$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Mm.eg.db")
  
  # Plot the reduced terms in different ways
  print( scatterPlot(simMatrix, reducedTerms))
  treemapPlot(reducedTerms)
  
}


# Enrichment analysis using GSEA on Hallmarks
# ______________________________________________________________

cat("<H4>Hallmarks enrichment using GSEA</H4>")

lFC <- sort( resLFC$log2FoldChange, decreasing = TRUE)
fgseaRes <- fgseaMultilevel( HALLMARKS_LIST_MOUSE, lFC, minSize=15, maxSize = 500)

fgseaRes_df = as.data.frame( fgseaRes)
fgseaRes_df = fgseaRes_df[ order( fgseaRes_df$padj, decreasing = FALSE), ]

DT::datatable( as.data.frame( fgseaRes_df), caption= "Hallmarks enrichment")

# Plot Enrichment graph for HALLMARK_INTERFERON_ALPHA_RESPONSE
cat("<H5>Analysis detail for HALLMARK_INTERFERON_ALPHA_RESPONSE</H5>")
enrich_hallmark_plot = plotEnrichment( HALLMARKS_LIST_MOUSE[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], lFC) + 
                          theme( axis.text.x = element_text( size = 8),
                                 axis.text.y = element_text( size = 8),
                                 axis.title.x = element_text( size = 8),
                                 axis.title.y = element_text( size = 8))
print( enrich_hallmark_plot)

ggsave( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_HALLMARK_INTERFERON_ALPHA_RESPONSE_Score.pdf")),
        plot = enrich_hallmark_plot,
        height = 3, width = 7, units = "cm")

# Plot a heatmap of the expression of the leading edge
leadingedge_zscore_df = t( scale( t( gene_expression_df[unlist( fgseaRes_df[ fgseaRes_df$pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE", "leadingEdge"]),])))
custom_height = 3.42+(0.1785714*nrow( leadingedge_zscore_df))
Heatmap( leadingedge_zscore_df, name ="Z-score",
         column_names_gp = gpar(fontsize = 8),
         row_names_gp = gpar(fontsize = 8),
         width = unit( 5, "cm"), height = unit( custom_height, "cm")) %>% 
  save_pdf( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_HALLMARK_INTERFERON_ALPHA_RESPONSE_Heatmap.pdf")),
            height = custom_height+4, width = 10, units = "cm")


# Enrichment analysis using GSEA on STAPP Pathways
# ______________________________________________________________

cat("<H4>STAPP Pathways enrichment using GSEA</H4>")

lFC <- sort( resLFC$log2FoldChange, decreasing = TRUE)
fgseaRes <- fgseaMultilevel( STAPP_PATHWAYS_LIST, lFC, minSize=15, maxSize = 1000)

fgseaRes_df = as.data.frame( fgseaRes)
fgseaRes_df = fgseaRes_df[ order( fgseaRes_df$padj, decreasing = FALSE), ]

DT::datatable( as.data.frame( fgseaRes_df), caption= "STAPP Pathway enrichment")

# Plot Enrichment graph for PATHWAY STAPP_REGENERATIVE
cat("<H5>Analysis detail for pathway STAPP_REGENERATIVE</H5>")
enrich_hallmark_plot = plotEnrichment( STAPP_PATHWAYS_LIST[["STAPP_REGENERATIVE"]], lFC) + 
  theme( axis.text.x = element_text( size = 8),
         axis.text.y = element_text( size = 8),
         axis.title.x = element_text( size = 8),
         axis.title.y = element_text( size = 8))
print( enrich_hallmark_plot)

ggsave( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_STAPP_REGENERATIVE_Score.pdf")),
        plot = enrich_hallmark_plot,
        height = 3, width = 7, units = "cm")

# Plot a heatmap of the expression of the leading edge
leadingedge_zscore_df = t( scale( t( gene_expression_df[unlist( fgseaRes_df[ fgseaRes_df$pathway == "STAPP_REGENERATIVE", "leadingEdge"]),])))
custom_height = 3.42+(0.1785714*nrow( leadingedge_zscore_df))
Heatmap( leadingedge_zscore_df, name ="Z-score", 
         column_names_gp = gpar(fontsize = 8),
         row_names_gp = gpar(fontsize = 8),
         width = unit( 5, "cm"), height = unit( custom_height, "cm")) %>% 
  save_pdf( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_STAPP_REGENERATIVE_Heatmap.pdf")),
            height = custom_height+4, width = 13, units = "cm")

# Plot Enrichment graph for PATHWAY STAPP_GOBLET
cat("<H5>Analysis detail for pathway STAPP_GOBLET</H5>")
enrich_hallmark_plot = plotEnrichment( STAPP_PATHWAYS_LIST[["STAPP_GOBLET"]], lFC) + 
  theme( axis.text.x = element_text( size = 8),
         axis.text.y = element_text( size = 8),
         axis.title.x = element_text( size = 8),
         axis.title.y = element_text( size = 8))
print( enrich_hallmark_plot)

ggsave( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_STAPP_GOBLET_Score.pdf")),
        plot = enrich_hallmark_plot,
        height = 3, width = 7, units = "cm")

# Plot a heatmap of the expression of the leading edge
leadingedge_zscore_df = t( scale( t( gene_expression_df[unlist( fgseaRes_df[ fgseaRes_df$pathway == "STAPP_GOBLET", "leadingEdge"]),])))
custom_height = 3.42+(0.1785714*nrow( leadingedge_zscore_df))
Heatmap( leadingedge_zscore_df, name ="Z-score",
         column_names_gp = gpar(fontsize = 8),
         row_names_gp = gpar(fontsize = 8),
         width = unit( 5, "cm"), height = unit( custom_height, "cm")) %>% 
  save_pdf( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_STAPP_GOBLET_Heatmap.pdf")),
            height = custom_height+4, width = 14, units = "cm")



# Enrichment analysis using GSEA on GOBP_CELL_POPULATION_PROLIFERATION
# _____________________________________________________________________

cat("<H4>GOBP_CELL_POPULATION_PROLIFERATION enrichment using GSEA</H4>")

lFC <- sort( resLFC$log2FoldChange, decreasing = TRUE)
fgseaRes <- fgseaMultilevel( list( GOBP_CELL_POPULATION_PROLIFERATION = GOBP_CELL_POPULATION_PROLIFERATION_geneset), lFC, minSize=15, maxSize = 2000)

fgseaRes_df = as.data.frame( fgseaRes)
fgseaRes_df = fgseaRes_df[ order( fgseaRes_df$padj, decreasing = FALSE), ]

DT::datatable( as.data.frame( fgseaRes_df), caption= "GOBP_CELL_POPULATION_PROLIFERATION enrichment")

# Plot Enrichment graph for GOBP_CELL_POPULATION_PROLIFERATION
cat("<H5>Analysis detail for GOBP_CELL_POPULATION_PROLIFERATION</H5>")
enrich_hallmark_plot = plotEnrichment( GOBP_CELL_POPULATION_PROLIFERATION_geneset, lFC) + 
  theme( axis.text.x = element_text( size = 8),
         axis.text.y = element_text( size = 8),
         axis.title.x = element_text( size = 8),
         axis.title.y = element_text( size = 8))
print( enrich_hallmark_plot)

ggsave( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_GOBP_CELL_POPULATION_PROLIFERATION_Score.pdf")),
        plot = enrich_hallmark_plot,
        height = 3, width = 7, units = "cm")

# Plot a heatmap of the expression of the leading edge
leadingedge_zscore_df = t( scale( t( gene_expression_df[unlist( fgseaRes_df[ fgseaRes_df$pathway == "GOBP_CELL_POPULATION_PROLIFERATION", "leadingEdge"]),])))
custom_height = 3.42+(0.1785714*nrow( leadingedge_zscore_df))
Heatmap( leadingedge_zscore_df, name ="Z-score", 
         column_names_gp = gpar(fontsize = 8),
         row_names_gp = gpar(fontsize = 8),
         width = unit( 5, "cm"), height = unit( custom_height, "cm")) %>% 
  save_pdf( filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Condition", "_", TEST_CONDITION, "_vs_", REFERENCE_CONDITION, "_GSEA_GOBP_CELL_POPULATION_PROLIFERATION_Heatmap.pdf")),
            height = custom_height+4, width = 13, units = "cm")


# Convert genes names to EntrezID and order the genes by decreasing LogFC
# ________________________________________________________________________

# Get the biomart ensembl reference
httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- biomaRt::useDataset("mmusculus_gene_ensembl", biomaRt::useMart("ensembl"))

# Get the entrezID of the genes
all_genes_entrezid <- biomaRt::getBM( filters="external_gene_name",
                         attributes=c("entrezgene_id", "external_gene_name"),
                         values= all_genes,
                         mart=mart)

# Look at NA in the entrezid and remove them
indNA = which(is.na(all_genes_entrezid$entrezgene_id))
all_genes_entrezid_noNA <- all_genes_entrezid[-indNA,]

# Look at the duplicates in entrezID and remove them
indnodup = which(duplicated( all_genes_entrezid_noNA$ entrezgene_id) == F)
all_genes_entrezid_noNA_nodup <- all_genes_entrezid_noNA[ indnodup,]

# Keep only the entreid that are no NA and not replicated in the DEseq result
lFC_entrezid <- resLFC$log2FoldChange[ -indNA]
lFC_entrezid <- lFC_entrezid[ indnodup]
names( lFC_entrezid) <- all_genes_entrezid_noNA_nodup$entrezgene_id

# Sort fold changes in decreasing order
lFC_entrezid <- sort( lFC_entrezid, decreasing = TRUE)

# Enrichment analysis using GSEA on KEGG Pathways
# ______________________________________________________________

cat("<H4>KEGG Pathway enrichment using GSEA</H4>")

# Execute GSEA KEGG analysis
gseaKEGG <- gseKEGG(geneList = lFC_entrezid,
                    organism = "mmu",
                    # nPerm = 1000, # default number permutations
                    minGSSize = 5, # minimum gene set size
                    pvalueCutoff = 0.2, # padj cutoff value
                    verbose = FALSE)
# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

# Replace the entrezID by symbols in core enrichment feature
gseaKEGG_results$core_enrichment_symbol = sapply( gseaKEGG_results$core_enrichment, function( id_string){
  token_ids = strsplit( id_string, "/", fixed= TRUE)[[1]]
  token_symbols = vector()
  for( token in token_ids){
    token_symbol = all_genes_entrezid_noNA_nodup[ which( all_genes_entrezid_noNA_nodup$entrezgene_id == token), "external_gene_name"]
    token_symbols = append( token_symbols, token_symbol)
  }
  return( paste( token_symbols, collapse = "/"))
})

DT::datatable( gseaKEGG_results, caption = "Enriched KEGG pathways")
