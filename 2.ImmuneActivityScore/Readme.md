# Score of anti-cancer immune activity

`[Data]`  
The underlying data that used for computation.

**annotation.RData**   contains the signature genes involved in the seven steps of cancer-immunity cycle.  
**human_gene2ensembl2symbol_list.RData**   is used for gene identifier transformation.  
**signature.GeneSymbol.list.RData**   contains IDs of signature genes from different databases.  
**genename_length3.0.Rdata**  contains gene length information.  


`[core scripts]`  
This collection of code is the core source code used for the calculation of activity levels of anti-cancer immunity signatures based on single sample Gene Set Enrichment Analysis (ssGSEA) algorithm.

**ssGSEAPermutation.R**  is to build the random matrix 'permutation.exp'.  
**ssgsea.core.R**  contains several subfunctions needed to calculate ES score of one geneset.  
**Count2TPM3.0.R**  for expression transformation.  


`[scripts for visualization]`   
This collection of code is used to change the format of results for visualization on the TIP platform.  

**changeStepNames.R**  is to change the names of different steps for clear visualization.  
**filterOutliers.R**  is to filter outliers in a matrix.  
**makeHeatmapData.R**  change the data format for heatmap visualization.  
**PCAScatterPlot.R**  change the data format for principal component analysis visualization.  
**BarPlotFormat.R**  change the data format for line chart of individual sample.  
**score_boxplot.R**  change the data format for boxplot of individual sample.  


`[Main scripts]`  
**ProcessMultipleSample.R**   This is the main function for analyzing expression profile with multiple samples.  
**ProcessSingleSample.R**   This is the main function for analyzing expression profile with a single sample.  
**immunityScore.server.R**  This is the main function for calculation of immune activity levels based on ssGSEA algorithm.  
