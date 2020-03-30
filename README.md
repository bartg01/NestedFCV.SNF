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

### Prepare discovery and validation cohorts:
```
discovery.data <- list()
validation.data <- list()
discovery.data$Transcriptome <- input.data[samples,]
discovery.data$Methylome <- input.data[samples,]
validation.data$Transcriptome <- input.data[-samples,]
validation.data$Methylome <- input.data[-samples,]
```
### Run nested fold cross-validation:

NestedCrossValidation function is the main function for running the algorithm, in this example the algorithm was run 10 times and tested from 2 to 20 clusters in an 5 inner and 10 outer nested fold cross-validation and it was parallelized in 10 threads.

`Nested.FCV.results <- NestedCrossValidation(discovery.data,10,20,5,10,10)`

### Plot NFCV results:

Plot inner FCV results:
```Plot.InnerValues.Corrected <- plot.Inner.Corrected(Nested.FCV.results,20)
Plot.InnerValues.Corrected$type <- factor(Plot.InnerValues.Corrected$type, levels = c('Real','Random','Corrected'))
Plot.InnerValues.Corrected$pos <- factor(Plot.InnerValues.Corrected$pos, levels=c('Datasets',"Corrected"))
ggplot(Plot.InnerValues.Corrected, aes(x=Clusters, y=avg, linetype=type, colour=Hyperparameter,
                                        group=interaction(type, Neighbors,Hyperparameter))) +
  facet_grid(pos ~ Neighbors, scales = "free_y", space = "free_y") +
  scale_color_igv() +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.1) +
  geom_line() +
  geom_point() +
  scale_x_discrete(limits=2:20) +
  labs(linetype = "Datasets",colour="Hyperparamer") +
  ylab("Proportions of Labels Recovery") +
  xlab("Number of Clusters") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"))
```  


Plot outer FCV results:
```
Plot.OuterValues.Corrected <- plot.Outer.Corrected(Nested.FCV.results,20)
Plot.OuterValues.Corrected <- as.data.frame(Plot.OuterValues.Corrected)
Plot.OuterValues.Corrected$type <- factor(Plot.OuterValues.Corrected$type, levels = c('Real','Random','Corrected'))
Plot.OuterValues.Corrected$pos <- factor(Plot.OuterValues.Corrected$pos, levels=c('Datasets',"Corrected"))
Plot.OuterValues.Corrected$avg <- as.numeric(as.character(Plot.OuterValues.Corrected$avg))
Plot.OuterValues.Corrected$sd <- as.numeric(as.character(Plot.OuterValues.Corrected$sd))
ggplot(Plot.OuterValues.Corrected, aes(x=Clusters, y=avg, linetype=type,
                                 group=interaction(type))) +
  facet_grid(pos ~ .) +
  scale_color_igv() +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.1) +
  geom_line() +
  geom_point() +
  scale_x_discrete(limits=2:20) +
  scale_y_continuous(limits=0:1) +
  labs(linetype = "Datasets",colour="Hyperparamer") +
  ylab("Proportions of Labels Recovery") +
  xlab("Number of Clusters") +
  theme_bw()
```


