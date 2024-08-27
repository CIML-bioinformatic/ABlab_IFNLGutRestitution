# ##################################################
# Global declarations and libraries for the analysis
# ##################################################

## @knitr load_general_deps

######## Loading of libraries

library( digest)
library( fs)
library( pander)

library( DESeq2)
library( apeglm)
library( vsn)
library( EnhancedVolcano)

library( org.Mm.eg.db)
library( DOSE)
library( pathview)
library( clusterProfiler)
library( AnnotationHub)
library( ensembldb)
library( biomaRt)
library( fgsea)

library( ggplot2)
library( gridExtra)
library( RColorBrewer)
library( pheatmap)
library( rrvgo)
library( ComplexHeatmap)
library( tidyHeatmap)

######## DEFINITION OF FUNCTIONS #########################

# ########################################################
# Convert a list of Human genes to their homolog in Mouse
# ########################################################

convertHumanGeneList <- function(x){
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  mousex <- unique( genesV2[, 2])
  return( mousex)
}

# ########################################################
# Convert a list of Human genes to their homolog in Mouse
# ########################################################

convertHumanGeneAllList <- function(x){
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  all_list = list()
  for( list_name in names( x)){
    current_list = x[[ list_name]]
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = current_list , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    mousex <- unique( genesV2[, 2])  
    all_list[[ list_name]] = mousex
  }
  
  return( all_list)
}


