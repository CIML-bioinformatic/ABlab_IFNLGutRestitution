# ##############################################################################
# Plots for provided genes lists
# ##############################################################################

#### MODULE SCORE

## @knitr heterogeneity_modules_scoring

MODULES_GENES = GENES_LIST

# Compute the score of the cells according to lists of genes
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


# # Plot the 'non-interactive' heatmap
# pheatmap( modulesScoreMat,
#           cluster_rows = FALSE,
#           cluster_cols = FALSE,
#           #          annotation_row = rowsAnnot,
#           annotation_col = colsAnnot,
#           labels_row = originalRowNames,
#           annotation_colors = list( Cluster = clustersColor),
#           show_colnames = FALSE);

scaleRowsHeatmap = TRUE

Heatmap( if(scaleRowsHeatmap) t(scale(t(modulesScoreMat))) else modulesScoreMat, 
         name = "Module Score", 
         col = RColorBrewer::brewer.pal(name = "Reds", n = 4),
         rect_gp = gpar(col = "white", lwd = 0), # Border of cells
         column_names_rot = 0, #75
         show_column_names = FALSE,
         cluster_rows = FALSE,
         cluster_row_slices = FALSE, # Do not cluster rows split categories (based on logFC) to keep order of slices under control
         cluster_columns = TRUE,
         left_annotation = NULL,
         top_annotation = HeatmapAnnotation( df = colsAnnot,
                                             col = list("Cluster" = clustersColor)),
         # Make groups of columns based on categories in dataframe
         column_split = colsAnnot[["Cluster"]],
         column_title = character(0), # Default: character(0), Disable: NULL
         column_title_rot = 0,
         # Make groups of rows
         row_split = NULL,
         row_title = character(0), # Default: character(0), Disable: NULL
         row_title_rot = 0)



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


### In case a reordering is required for figures ###
# Idents(sc10x) = factor(Idents(sc10x), levels = sort(levels(Idents(sc10x)))) 
####################################################

## @knitr heterogeneity_modules_expression_violin

# Violinplot for module values by cluster
invisible( lapply( paste0(names(MODULES_GENES), 1), violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor, yLabel = "Score", addStats = FALSE, trimTitle = 1));

cat(" \n \n"); # Required for '.tabset'





## @knitr heterogeneity_modules_expression_violin_byHTO

# Violinplot for module values by cluster
invisible( lapply( paste0(names(MODULES_GENES), 1), violinFeatureByMetadata, seuratObject = sc10x, yLabel = "Score", addStats = FALSE, trimTitle = 1));

cat(" \n \n"); # Required for '.tabset'




## @knitr heterogeneity_modules_expression_violin_HTO_Cluster_TWAnova
require(xtable)
library(Hmisc)

for(currentModuleName in  paste0(names(MODULES_GENES), 1))
{
  cat(paste0("\n\n##### ",  currentModuleName, "\n\n"))
  
  featureValues = FetchData( sc10x, vars = currentModuleName);
  data = data.frame(moduleScore = sc10x[[currentModuleName, drop = TRUE]], HTO = factor(sc10x[["factorHTO", drop = TRUE]]), Cluster = factor(Idents(sc10x)))

  # Count number of cells in each group
  countTable = tapply(data[["moduleScore"]], as.list(data[c("HTO", "Cluster")]), length)
  cat(" \n \n");
  formatedTable = 
  #names(formatedTable) = colnames(data)
  # print( kable( lapply( countTable, function(x)
  # {
  #   kableExtra::cell_spec(x, color = ifelse(x<30, "red", "black"))
  # }), caption = "Groups counts"))
  emphasize.strong.cells(which(countTable<30, arr.ind = TRUE))
  cat(pander(countTable, split.table=500, split.cell = 20, caption="Cell count"))
  cat(" \n \n");
  
  # Create model
  mod <- aov(moduleScore ~ HTO * Cluster, data = data)
  #cat(" \n \n");
  #cat(pander(mod))
  #cat(" \n \n");
  
  
  # Flag to check whether we can do a parametric test or not (depends on following tests)
  paramTestOK = TRUE
  
  # Check whether model satisfies residuals normality (if less than ~30 observation in each group) 
  cat("###### Testing Residuals normality")
  library(car)
  cat(" \n \n");
  qqPlot(mod$residuals, id = FALSE)
  hist(mod$residuals)
  cat(" \n \n");
  # Isolate residual (and eventually trim) for shapiro test (cannot be done on more thant 5000, tends to always reject null hypothesis on large samples)
  residuals = mod$residuals
  if(length(residuals)>=5000)
  {
    cat("\n\nWARNING: Too many samples for shapiro test, trimming to 5000 first residuals only.\n")
    residuals = residuals[1:5000]
  }
  shapiTest = shapiro.test(residuals)
  cat(pander(broom::tidy(shapiTest), split.table=500, split.cell = 20)) # If Pval <0.05 the null hypothesis of normality is rejected (not good)
  cat(" \n \n");
  if(shapiTest[["p.value"]]<0.05)
  {
    cat(" \n \n");
    cat("\nResiduals normality null hypothesis is rejected by Shapiro test...\n")
    cat(" \n \n");
    paramTestOK = FALSE
  }
  
  
  # Test equality of variance between groups (homoscedasticity)
  cat("###### Testing equality of variance between groups")
  cat(" \n \n");
  plot(mod, which = 3)
  cat(" \n \n");
  levTest = leveneTest(mod) 
  cat(" \n \n");
  cat(pander(broom::tidy(levTest), split.table=500, split.cell = 20)) # If Pval <0.05 the null hypothesis of variance equality is rejected (not good)
  cat(" \n \n");
  if(levTest[["Pr(>F)"]][1]<0.05)
  {
    cat(" \n \n");
    cat("\nResiduals heteroscedasticity null hypothesis is rejected by Levene test...\n")
    cat(" \n \n");
    paramTestOK = FALSE
  }
  
  cat("###### 2-Way Anova test")
  cat(" \n \n");
  cat(pander(anova(mod), split.table=500, split.cell = 20))
  cat(" \n \n");
  
  cat("###### TukeyHSD post-hoc tests")
  
  oldpar = par(mar = c(5.1, 20, 4.1, 2.1))
  
  tukeyTest=TukeyHSD(mod, which = "HTO")
  cat(" \n \n");
  cat(pander(broom::tidy(tukeyTest), split.table=500, split.cell = 20))
  cat(" \n \n");
  plot(TukeyHSD(mod, which = "HTO"), las = 2) # rotate x-axis ticks
  cat(" \n \n");
  
  tukeyTest=TukeyHSD(mod, which = "Cluster")
  cat(" \n \n");
  cat(pander(broom::tidy(tukeyTest), split.table=500, split.cell = 20))
  cat(" \n \n");
  plot(TukeyHSD(mod, which = "Cluster"), las = 2) # rotate x-axis ticks
  cat(" \n \n");
  
#  tukeyTest=TukeyHSD(mod, which = "HTO:Cluster")
#  cat(" \n \n");
#  cat(pander(broom::tidy(tukeyTest), split.table=500, split.cell = 20))
#  cat(" \n \n");
#  plot(TukeyHSD(mod, which = "HTO:Cluster"), las = 2) # rotate x-axis ticks
#  cat(" \n \n");
  
  par(oldpar)
}

cat(" \n \n"); # Required for '.tabset'




## @knitr heterogeneity_modules_expression_KWallis_HTO_byCluster

for(currentModuleName in  paste0(names(MODULES_GENES), 1))
{
  cat(paste0("\n\n##### ",  currentModuleName, " {.tabset .tabset-fade .tabset-pills}\n\n"))
  
  featureValues = FetchData( sc10x, vars = currentModuleName);
  dataAll = data.frame(moduleScore = sc10x[[currentModuleName, drop = TRUE]], HTO = factor(sc10x[["factorHTO", drop = TRUE]]), Cluster = factor(Idents(sc10x)))
  
  for(currentCluster in levels(factor(Idents(sc10x))))
  {
    cat(paste0("\n\n###### ",  currentCluster, "\n\n"))
    
    # Isolate data for current subgroup
    data = dataAll[ dataAll[["Cluster"]] %in% currentCluster,]
    
    # Perform KT test for global differences
    kt = kruskal.test(moduleScore ~ HTO, data = data)
    
    cat(" \n \n");
    cat(pander(kt, split.table=500, split.cell = 20))
    cat(" \n \n");
    
    
    # Perform post-hoc test for differences group by group
    dt = dunn_test(data = data, moduleScore ~ HTO)
    
    cat(" \n \n");
    cat(pander(dt, split.table=500, split.cell = 20))
    cat(" \n \n");
    
  }
}

cat(" \n \n"); # Required for '.tabset'






## @knitr heterogeneity_modules_expression_KWallis_Cluster_byHTO

for(currentModuleName in  paste0(names(MODULES_GENES), 1))
{
  cat(paste0("\n\n##### ",  currentModuleName, " {.tabset .tabset-fade .tabset-pills}\n\n"))
  
  featureValues = FetchData( sc10x, vars = currentModuleName);
  dataAll = data.frame(moduleScore = sc10x[[currentModuleName, drop = TRUE]], HTO = factor(sc10x[["factorHTO", drop = TRUE]]), Cluster = factor(Idents(sc10x)))
  
  for(currentHTO in levels(factor(dataAll[["HTO"]])))
  {
    cat(paste0("\n\n###### ",  currentHTO, "\n\n"))
    
    # Isolate data for current subgroup
    data = dataAll[ dataAll[["HTO"]] %in% currentHTO,]
    
    # Perform KT test for global differences
    kt = kruskal.test(moduleScore ~ Cluster, data = data)
    
    cat(" \n \n");
    cat(pander(kt, split.table=500, split.cell = 20))
    cat(" \n \n");
    
    
    # Perform post-hoc test for differences group by group
    dt = dunn_test(data = data, moduleScore ~ Cluster)
    
    cat(" \n \n");
    cat(pander(dt, split.table=500, split.cell = 20))
    cat(" \n \n");
    
  }
}

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
  pathCurrentModule = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_", ifelse(exists("useReduction"), useReduction, "umap")), moduleName);
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
  pathCurrentModule = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_violin"), moduleName);
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
  
}));





## @knitr heterogeneity_aucell

BiocManager::install("AUCell")
library("AUCell")

MODULES_GENES = GENES_LIST

expMat = as.matrix(GetAssayData(sc10x, slot = "count"))

cellsByHTO = split( colnames(expMat), sc10x[["factorHTO"]])
expMatByHTO = lapply( cellsByHTO, function(x){ expMat[, x] })


# Calculate enrichment scores
auc_rankings_HTOlist = lapply(expMatByHTO, AUCell_buildRankings, GENES_LIST )#, aucMaxRank=nrow(cells_rankings)*0.05)
auc_HTOlist = lapply(auc_rankings_HTOlist,function(x) {AUCell_calcAUC(GENES_LIST, x)})


AUCell_calcAUC(GENES_LIST, auc_rankings_list[[1]])

# Optional: Set the assignment thresholds
par(mfrow=c(3,3))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)


