# ########################################################
# This script wraps up analysis (clean, save, render, ...)
# ########################################################




## @knitr final_saveSessionImage

# Save an image of the current R session containing all objects
sessionImagePath = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sessionImage_final.RDATA"));
save(list = ls( all.names = TRUE), 
     file = sessionImagePath, 
     envir = environment());
#save.image( sessionImagePath); # Saves everything in .Globalenv which is not what we want depending on rendering context

# Convert character metadata to factor and save binary file of Seurat object
characterMetadata = colnames(sc10x[[]])[sapply(sc10x[[]], is.character)]
sc10x[[characterMetadata]]=lapply(sc10x[[characterMetadata]], as.factor)

seuratObjectPath = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "seuratObject_final.RDS"));
saveRDS( object = sc10x, file = seuratObjectPath);



