# clusterMLD: Clustering Multivariate Longitudinal Data
Support to efficiently cluster multivariate longitudinal data with sparse (very little observed points for each subject) and irregular (the observation occassions are not aligned across subjects) observations. 

Longitudinal data with multiple outcomes is supported as well, where some of them could be pure noise (non-distinguishable). 

Support the case with potentially unbalanced cluster size, i.e., the number of subjects in some clusters are way outnumbered by the others.

In sum, the package is capable of clustering **sparse, irregular, unbalanced, and multivariate continuous** longitudinal data  

## Installation
To install the package: 
```r
devtools::install_github("junyzhou10/clusterMLD")
```

## Usage
Use main function `LongDataCluster(x, Y, id, ...)` to cluster longitudinal data in long format. Parallel computing is supported by specifying `parallel = TRUE`. 

`DendroPlot(Cluster.object)` yields corresponding dendrogram, where Cluster.object is the output from `LongDataCluster()`

`MeanPlot(Cluster.object)` yields corresponding mean curves of each detected cluster.

## Exploration
The author created an R shiny app for a better illustration/visualization of the package. Please refer to [clusterMLD_ShinyApp](https://github.com/junyzhou10/clusterMLD_ShinyApp) for more details.

The App is published at https://junyzhou.shinyapps.io/clusterMLD_ShinyApp/, with two toy examples uploaded already. Please note that it will take more time to run parallel online than local.

## Limitations
The current version mainly support for the continuous outcomes. Algorithm for other outcome types, such as binary and time-to-event, are now under development. 
