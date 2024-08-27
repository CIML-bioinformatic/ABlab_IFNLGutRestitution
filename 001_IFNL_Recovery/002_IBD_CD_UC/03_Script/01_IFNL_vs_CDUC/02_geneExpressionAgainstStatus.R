# ##################################################
# LOkk at the gene expression of some gene agains
# the status of smaples
# ##################################################

## @knitr gene_expression_against_status

# Look at the expression of GSDMC in the different status
# ________________________________________________________

cat("<H4>Expression of GSDMC</H4>")

# -- Create dataframe with meta-data and gene expression
GSDMC_expression_df = sample_metadata_df[ c( CONTROL_SAMPLES, QUIESCENT_SAMPLES, ACTIVE_SAMPLES), c("sample.name", "condition.name", "status", "group")]
GSDMC_expression_df$GSDMC = t( expression_matrix_df[ "GSDMC", c( CONTROL_SAMPLES, QUIESCENT_SAMPLES, ACTIVE_SAMPLES)])

# -- Write data to file
write.table( GSDMC_expression_df, file = file.path( PATH_ANALYSIS_OUTPUT, "GSDMC_expression.csv"),
           col.names = NA, row.names = TRUE, quote = FALSE, sep=",")

# -- Display the gene expression agains the statuse
ggplot( GSDMC_expression_df) +
  geom_boxplot( aes( x=status, y=GSDMC, fill = status)) +
  geom_jitter( aes( x=status, y=GSDMC, color = group), width = 0.2, size = 3) +
  ggtitle( "Expression of GSDMC across status (boxplot)") + 
  theme_minimal()

ggplot( GSDMC_expression_df) +
  geom_violin( aes( x=status, y=GSDMC, fill = status)) +
  geom_jitter( aes( x=status, y=GSDMC, color = group), width = 0.2, size = 3) +
  ggtitle( "Expression of GSDMC across status (violin plot)") + 
  theme_minimal()

#-- Compute the statistics on gene expression againt status using non-parametric test
cat("<H5>Kruskall-Wallis test and Dunn post-hoc tests</H5>")
GSDMC_kruskaltest = kruskal_test( GSDMC ~ status, data = GSDMC_expression_df)
kable( GSDMC_kruskaltest, format = "html", caption = "Krukal-Wallis test of GSDMC gene expression against status ")
GSDMC_kruskaltest %>%
  kbl() %>%
  kable_styling()

GSDMC_dunntest = dunn_test( GSDMC ~ status, data = GSDMC_expression_df, p.adjust.method = "bonferroni") 
GSDMC_dunntest %>%
  kbl() %>%
  kable_styling()

# Look at the expression of ZBP1 in the different status
# ________________________________________________________

cat("<H4>Expression of ZBP1</H4>")

# -- Create dataframe with meta-data and gene expression
ZBP1_expression_df = sample_metadata_df[ c( CONTROL_SAMPLES, QUIESCENT_SAMPLES, ACTIVE_SAMPLES), c("sample.name", "condition.name", "status", "group")]
ZBP1_expression_df$ZBP1 = t( expression_matrix_df[ "ZBP1", c( CONTROL_SAMPLES, QUIESCENT_SAMPLES, ACTIVE_SAMPLES)])

# -- Write data to file
write.table( ZBP1_expression_df, file = file.path( PATH_ANALYSIS_OUTPUT, "ZBP1_expression.csv"),
           col.names = NA, row.names = TRUE, quote = FALSE, sep=",")

# -- Display the gene expression agains the statuse
ggplot( ZBP1_expression_df) +
  geom_boxplot( aes( x=status, y=ZBP1, fill = status)) +
  geom_jitter( aes( x=status, y=ZBP1, color = group), width = 0.2, size = 3) +
  ggtitle( "Expression of ZBP1 across status (boxplot)") + 
  theme_minimal()

ggplot( ZBP1_expression_df) +
  geom_violin( aes( x=status, y=ZBP1, fill = status)) +
  geom_jitter( aes( x=status, y=ZBP1, color = group), width = 0.2, size = 3) +
  ggtitle( "Expression of ZBP1 across status (violin plot)") + 
  theme_minimal()

#-- Compute the statistics on gene expression againt status using non-parametric test
cat("<H5>Kruskall-Wallis test and Dunn post-hoc tests</H5>")
ZBP1_kruskaltest = kruskal_test( ZBP1 ~ status, data = ZBP1_expression_df)
ZBP1_kruskaltest %>%
  kbl() %>%
  kable_styling()

ZBP1_dunntest = dunn_test( ZBP1 ~ status, data = ZBP1_expression_df, p.adjust.method = "bonferroni") 
ZBP1_dunntest %>%
  kbl() %>%
  kable_styling()

