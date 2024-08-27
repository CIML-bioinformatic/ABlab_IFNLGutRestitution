# ########################################################
# This script applies the selected clustering algorithm to
# the dataset and plots all 
# ########################################################


###################
# IDENTIFY CLUSTERS
###################

# ..............................................................................
## @knitr heterogeneity_identifyClusters
# ..............................................................................

# Identify clusters of cells by graph approach
nbPC_findclusters=FINDCLUSTERS_USE_PCA_NBDIMS


nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(RECOMPUTE_CLUSTERING)
{
  cat("\n\n(Re-)computing clustering...\n")
  
  if(FINDCLUSTERS_USE_PCA_NBDIMS>nbPC)
  {
    warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'findclusters' (", FINDCLUSTERS_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
    nbPC_findclusters = nbPC
  }
  
  sc10x <- FindNeighbors(object    = sc10x,
                         k.param   = FINDNEIGHBORS_K, 
                         reduction = "pca",
                         dims      = 1:nbPC_findclusters,
                         verbose   = .VERBOSE);
  
  sc10x <- FindClusters(object             = sc10x,
                        resolution         = FINDCLUSTERS_RESOLUTION,
                        algorithm          = FINDCLUSTERS_ALGORITHM,
                        temp.file.location = "/tmp/",
                        verbose            = .VERBOSE);
  
  # Rename the clusters if required
  if( exists( "RENAME_CLUSTERS") && SAMPLE_NAME %in% names( RENAME_CLUSTERS))
    {
      Idents( sc10x) = "seurat_clusters"
      sc10x = AddMetaData( sc10x, metadata = Idents( sc10x), col.name = "original_seurat_clusters")
      new_names = Idents( sc10x)
      for( cluster_id in levels( Idents( sc10x)))
      {
        new_names[ which( Idents( sc10x) == cluster_id)] = RENAME_CLUSTERS[[ SAMPLE_NAME]][[ cluster_id]]
      }
    sc10x = AddMetaData( sc10x, metadata = new_names, col.name = "seurat_clusters")
    
    Idents( sc10x) = "seurat_clusters"
  }
  
  # Define a set of colors for clusters (based on ggplot default)
  clustersColor = hue_pal()( nlevels( Idents( sc10x)));
  names( clustersColor) = levels( Idents( sc10x));
  
} else
{
  cat("\n\nUsing predefined clustering...\n")
}


# Save cells cluster identity
write.table( data.frame(sc10x[["numID"]], identity = Idents(sc10x)), 
             file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsClusterIdentity.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Also save cluster color attribution for reference
# Save cells cluster identity as determined with 'FindClusters'
write.table( clustersColor, 
             file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "clustersColor.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = FALSE,
             sep="\t");


# Compute the number of cells by cluster
clustersCount = as.data.frame( table( Cluster = Idents(sc10x)), responseName = "Total");

# Compute the cluster counts sliced by sample
clustersCountBySample = table(cbind(Idents( sc10x), sc10x[["factorHTO"]]))

# # Compute the cluster counts sliced by condition
clustersCountByCondition = table(cbind(Idents( sc10x), sc10x[["orig.ident"]]))


# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( GetAssayData( sc10x)), # expression values
                                   1,                                # by rows
                                   tapply,                           # apply by group
                                   INDEX = Idents( sc10x),           # clusters IDs
                                   mean,                             # summary function
                                   simplify = FALSE));               # do not auto coerce
# Save it as 'tsv' file
write.table( geneExpByCluster, 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCluster.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");


# Compute a matrix of average expression value by sample (for each gene)
geneExpBySample = do.call( rbind, 
                           apply( as.matrix( GetAssayData( sc10x, 
                                                           slot = "data", 
                                                           assay="RNA")),  # normalized expression values
                                  1,                                       # by rows
                                  tapply,                                  # apply by group
                                  INDEX = sc10x[["factorHTO", drop = TRUE]],                  # clusters IDs
                                  mean,                                    # summary function
                                  simplify = FALSE));                      # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpBySample, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionBySample.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");


# Compute a matrix of average expression value by sample (for each gene)
geneExpByCondition = do.call( rbind, 
                           apply( as.matrix( GetAssayData( sc10x, 
                                                           slot = "data", 
                                                           assay="RNA")),  # normalized expression values
                                  1,                                       # by rows
                                  tapply,                                  # apply by group
                                  INDEX = sc10x[["orig.ident", drop = TRUE]],                  # clusters IDs
                                  mean,                                    # summary function
                                  simplify = FALSE));                      # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByCondition, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCondition.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");


# Compute a matrix of average expression value by cluster and Sample (for each gene)
geneExpByClusterAndSample = do.call( rbind, 
                                     apply( as.matrix( GetAssayData( sc10x, 
                                                                     slot = "data", 
                                                                     assay="RNA")),  # normalized expression values
                                            1,                                       # by rows
                                            tapply,                                  # apply by group
                                            INDEX = paste( Idents( sc10x),
                                                           sc10x[["factorHTO", drop = TRUE]], 
                                                           sep = "_"),               # clusters + samples IDs
                                            mean,                                    # summary function
                                            simplify = FALSE));                      # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByClusterAndSample, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByClusterAndSample.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");


# Compute a matrix of average expression value by cluster and Sample (for each gene)
geneExpByClusterAndCondition = do.call( rbind, 
                                     apply( as.matrix( GetAssayData( sc10x, 
                                                                     slot = "data", 
                                                                     assay="RNA")),  # normalized expression values
                                            1,                                       # by rows
                                            tapply,                                  # apply by group
                                            INDEX = paste( Idents( sc10x),
                                                           sc10x[["orig.ident", drop = TRUE]], 
                                                           sep = "_"),               # clusters + samples IDs
                                            mean,                                    # summary function
                                            simplify = FALSE));                      # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByClusterAndCondition, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByClusterAndCondition.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");






