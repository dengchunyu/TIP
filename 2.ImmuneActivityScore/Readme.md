# Calculating score of anti-cancer immune activity

`[Data]`  
The underlying data that used for computation.

**annotation.RData**   A matrix contains the signature genes involved in the seven steps of cancer-immunity cycle, with gene IDs from different databases, number of steps, direction of gene function and corresponding cell type.  
**human_gene2ensembl2symbol_list.RData**   A list used for gene identifier transformation.  
**signature.GeneSymbol.list.RData**   A list contains gene symbol of signature genes.  
**genename_length3.0.Rdata**  A data.frame consist of Entrze geneID, gene name and gene length required for normalizing each gene count to TPM.  


`[core scripts]`  
This collection of code is the core source code used for the calculation of activity levels of anti-cancer immunity signatures based on single sample Gene Set Enrichment Analysis (ssGSEA) algorithm.

**ssGSEAPermutation.R**  is to build the random expression matrix 'permutation.exp' that required for normalizing ES score of real samples.  
**ssgsea.core.R**  contains several subfunctions required for calculating ES score of one geneset.  
**Count2TPM3.0.R**  for normalization from gene count to TPM.  


`[scripts for visualization]`   
This collection of code is used to change the format of results for visualization on the TIP platform.  

**changeStepNames.R**  is to change the names of every steps for clear visualization.  
**filterOutliers.R**  is to filter outliers in a matrix.  
**makeHeatmapData.R**  change the data format for heatmap visualization.  
**PCAScatterPlot.R**  change the data format for principal component analysis visualization.  
**BarPlotFormat.R**  change the data format for line chart visualization of individual sample.  
**score_boxplot.R**  change the data format for boxplot visualization of individual sample.  


`[Main scripts]`  
**ProcessMultipleSample.R**   This is the main function for analysising expression profile with multiple samples, required to complete 'immunityScore.server.R'.  
**ProcessSingleSample.R**  This is the main function for analysising expression profile with a single sample, required to complete 'immunityScore.server.R'.  
**immunityScore.server.R**  This is the main function for calculation of immune activity levels based on ssGSEA algorithm.  
