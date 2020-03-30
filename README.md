# NestedFCV.SNF

Optimize Similarity Network Fusion algorithm (https://www.nature.com/articles/nmeth.2810) hyper-parameters using a nested fold cross-validation (NFCV). The optimization is based on the proportion of samples assigments to the same clusters, and the proportions might be corrected by the null distribution of random assigments.

### Prerequisites

`install.packages(c("SNFtool","ggsci","ggplot2"))`

## Getting Started
 
```
install.packages("devtools")
devtools::install_github("bartg01/NestedFCV.SNF")
library(NestedFCV.SNF)
```

## Prepare Input Datasets

Multiple datasets must be used in the analysis. Each layer of information must be formated in a matrix with features in rows and samples in columns. Every layer must present the same samples and in the the same order (column names must be unique and identify each sample).

```
input.data <- list()
input.data$Transcriptome <- transcriptome.values[,order(colnames(transcriptome.values))]
input.data$Methylome <- methylome.values[,order(colnames(methylome.values))]
```

## Running Example

In this example NFCV.SNF is used in a discovery cohort to optimize the hyper-parameters, and then the trained model is applied to cluster the samples of a validation cohort.<br/><br/>

###Prepare discovery and validation cohorts:
```
discovery.data <- list()
validation.data <- list()
discovery.data$Transcriptome <- input.data[samples,]
discovery.data$Methylome <- input.data[samples,]
validation.data$Transcriptome <- input.data[-samples,]
validation.data$Methylome <- input.data[-samples,]
```
###Running nested fold cross-validation:



