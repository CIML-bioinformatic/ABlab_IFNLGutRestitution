# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# MODULES ANALYSIS
##################

## @knitr heterogeneity_modules

# Just remind the warning for genes names not in object, or modules that were transfered to individual monitoring of genes
if(any( is.na( matchModulesGenes)))
{
    warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:",
                    paste( modulesGenesNotFound, collapse=" - ")));
}

if(any(modulesToTransfer))
{
    warning( paste0( "Following Module(s) contained very few genes (<",
                     MONITORED_GENES_SMALL_MODULE,
                     "). These genes were transfered to 'Monitored genes' to be analyzed individually: ",
                     paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
}




## @knitr heterogeneity_modules_scoring

# Compute the score of the cells according to group of monitored genes
for( moduleName in names( MODULES_GENES))
{
  message(moduleName);
  if( length( MODULES_GENES[[moduleName]]) == 0)
  {
    warning( paste0( "List of genes in module '", moduleName, "' is empty, ignoring..."));
  } else
  {
    sc10x <- AddModuleScore( object = sc10x,
                             features = MODULES_GENES[ moduleName], # Must be a list
                             ctrl = MODULES_CONTROL_SIZE,           #length(MODULES_GENES[[ moduleName]]),
                             name = moduleName,
                             seed = SEED);
  }
}




## @knitr heterogeneity_modules_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = topMarkersDF[["gene"]])));

# Get the matrix of module scores (not expression) for each cell and associated clusters from Seurat object
modulesScoreMat = t( as.matrix( sc10x[[ paste0(names( MODULES_GENES),1)]]));
clusterID = Idents( sc10x);

# Remove the extra numeric character added to modules names by Seurat
rownames( modulesScoreMat) = substr( rownames( modulesScoreMat), 1, nchar( rownames( modulesScoreMat))-1);

# Select genes in modules and reorder cells to group clusters together
clusterOrdering = order( clusterID);

modulesScoreMat = modulesScoreMat[, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Prepare rows and columns annotation bars (module and cluster respectively)
#rowsAnnot = data.frame( Module = names( MODULES_GENES));
colsAnnot = data.frame( Cluster = clusterID);


# Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
originalRowNames = rownames( modulesScoreMat);
originalColNames = colnames( modulesScoreMat);
rownames( modulesScoreMat) = make.unique( originalRowNames);
colnames( modulesScoreMat) = make.unique( originalColNames);
#rownames( rowsAnnot) = rownames( modulesScoreMat);
rownames( colsAnnot) = colnames( modulesScoreMat);

# Plot the 'non-interactive' heatmap
pheatmap( modulesScoreMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
#          annotation_row = rowsAnnot,
          annotation_col = colsAnnot,
          labels_row = originalRowNames,
          annotation_colors = list( Cluster = clustersColor),
          show_colnames = FALSE);




## @knitr heterogeneity_modules_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot scores of modules on dimensionality reduction figures
invisible( lapply( names(MODULES_GENES), function(featureName)
{
  print( FeaturePlot( sc10x, features = paste0(featureName, "1"), reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
           ggtitle( label = featureName) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none"));
}));

cat(" \n \n"); # Required for '.tabset'




## @knitr heterogeneity_modules_expression_violin

# Violinplot for module values by cluster
invisible( lapply( paste0(names(MODULES_GENES), 1), violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor, yLabel = "Score", addStats = FALSE, trimTitle = 1));

cat(" \n \n"); # Required for '.tabset'




## @knitr heterogeneity_modules_expression_projection_pngFile
# Plot expression values of individual module genes as png files (not in report)
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console

  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_", ifelse(exists("useReduction"), useReduction, "umap")), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);

  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);

    print(FeaturePlot( sc10x, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
      theme( axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none"))

    dev.off(); # Close file descriptor
  }));

}));




## @knitr heterogeneity_modules_expression_violin_pngFile
# Plot expression values of monitored genes as violinplot in a png file for each cluster (TODO: message if list empty)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console

  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_violin"), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);

  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);

    violinFeatureByCluster(featureName, seuratObject = sc10x, clustersColor = clustersColor);

    dev.off(); # Close file descriptor
  }))

}))


