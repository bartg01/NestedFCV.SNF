library(SNFtool)
library(ggplot2)
library(ggsci)

#SNF Nested Cross-Validation main functions
NestedCrossValidation <- function(values,repetitions,numClusters,InnerFC,OuterFC,threads) {
  require('NMF')
  require('SNFtool')

  KCheck <- c(10,15,20,25,30)
  names(KCheck) <- c("K10","K15","K20","K25","K30")
  alphaCheck <- c(0.3,0.4,0.5,0.6,0.7,0.8)
  names(alphaCheck) <- c("Alpha3","Alpha4","Alpha5","Alpha6","Alpha7","Alpha8")

  resultsFC <- list()
  Datasets <- values

  ####################################
  #### LOOP FOR 10FCV REPETITIONS ####
  for (i in 1:repetitions) {

    print(paste0("Starting Nested CrossValidation ",i))

    ####################################################################
    #### FCV Preparation (Sort matrices by samples randomization) ####
    Datasets[[1]]<-Datasets[[1]][,sample(ncol(Datasets[[1]]))]
    for (j in 2:length(Datasets)) {Datasets[[j]] <- Datasets[[j]][,match(colnames(Datasets[[1]]),colnames(Datasets[[j]]))]}
    ####################################################################

    #########################################################
    #### Create outer equally size folds & Prepare datasets ####
    folds <- cut(seq(1,ncol(Datasets[[1]])),breaks=OuterFC,labels=FALSE)
    Datasets.OuterFCV <- list()
    for (j in 1:OuterFC) {
      Datasets.OuterFCV[[j]] <- list()
      Datasets.OuterFCV[[j]]$train <- list()
      Datasets.OuterFCV[[j]]$test <- list()
      index <- which(folds==j,arr.ind=TRUE)
      for (k in 1:length(Datasets)) {
        Datasets.OuterFCV[[j]]$train[[k]] <- t(Datasets[[k]][,-index])
        Datasets.OuterFCV[[j]]$test[[k]] <- t(Datasets[[k]][,index])
      }
    }
    #########################################################

    ####################################
    #### Perform inner 5FCV Model Optimization ####
    print(paste0("Performing Model Optimization (InnerFCV) SNF Clustering ",i," ..."))
    resultsModelOptimization <- mclapply(c(1:OuterFC),mainModelOptimization,values=Datasets.OuterFCV,
                                         numClusters=numClusters, KCheck=KCheck,
                                         alphaCheck=alphaCheck, InnerFC = InnerFC,
                                         mc.cores = threads)
    ###################################
    #### Select Optimal Parameters ####
    SelectParameters <- list()
    for (j in 1:length(resultsModelOptimization)) { #Each FCV from OuterFCV
      for (k in 1:length(resultsModelOptimization[[j]])) { #Each FCV from InnerFCV
        for (l in 1:length(resultsModelOptimization[[j]][[k]]$Stability)) { #K Values
          for (m in 1:length(resultsModelOptimization[[j]][[k]]$Stability[[l]])) { #Alpha Values
            for (n in 1:length(resultsModelOptimization[[j]][[k]]$Stability[[l]][[m]])) { #Clusters
              val <- paste0(names(resultsModelOptimization[[j]][[k]]$Stability)[l],'_',names(resultsModelOptimization[[j]][[k]]$Stability[[l]])[m],'_',n)
              SelectParameters[[val]] <- c(SelectParameters[[val]],resultsModelOptimization[[j]][[k]]$Stability[[l]][[m]][[n]]$RealPropAVG)
            }
          }
        }
      }
    }
    ParamValues <- sapply(SelectParameters,mean)
    SelectedValues <- names(which(ParamValues==max(ParamValues)))
    ########################################
    ######### Perform outer 10 FCV #########
    print(paste0("Performing outer OuterFCV SNF Clustering ",i," ..."))
    resultsOuterFCV <- mclapply(c(1:OuterFC),mainOuterFCV,values=Datasets.OuterFCV,
                                  numClusters=numClusters, SelectedValues=SelectedValues,
                                  KCheck=KCheck,alphaCheck=alphaCheck,
                                  mc.cores = threads)
    ########################################
    resultsFC[[i]] <- list(modelOpt = resultsModelOptimization, modelOptParams = ParamValues, selectValues = SelectedValues, outValidation = resultsOuterFCV)
  }
  ##################### END REPETITIONS ########################
  return(resultsFC)
}
mainModelOptimization <- function(k,values,numClusters,KCheck,alphaCheck,InnerFC) {
  ####################################################################
  #### InnerFCV Preparation (Sort matrices by samples randomization) ####
  values[[k]]$train[[1]]<-values[[k]]$train[[1]][sample(nrow(values[[k]]$train[[1]])),]
  for (j in 2:length(values[[k]]$train)) {values[[k]]$train[[j]] <- values[[k]]$train[[j]][match(rownames(values[[k]]$train[[1]]),rownames(values[[k]]$train[[j]])),]}
  foldsModel <- cut(seq(1,nrow(values[[k]]$train[[1]])),breaks=InnerFC,labels=FALSE)
  ####################################################################

  #############################################
  #### Perform InnerFCV for model optimization ####
  modelOptimization <- list()
  for (j in 1:InnerFC) {
    modelOptimization[[j]] <- list()
    train <- list()
    train$norm.real <- list()
    train$dist.real <- list()
    train$norm.rand <- list()
    train$dist.rand <- list()
    test <- list()
    test$norm.real <- list()
    test$dist.real <- list()
    test$norm.rand <- list()
    test$dist.rand <- list()
    #Prepare Datasets
    model.index <- which(foldsModel==j,arr.ind=TRUE)
    for (l in 1:length(values[[k]]$train)) {
      train.temp <- values[[k]]$train[[l]][-model.index,]
      test.temp <- values[[k]]$train[[l]][model.index,]
      #Real Datasets (Normalize Datasets & Calculate Distances)
      train$norm.real[[l]] <- standardNormalization(train.temp)
      train$dist.real[[l]] <- as.matrix(dist(train$norm.real[[l]]))
      test$norm.real[[l]] <- standardNormalization(test.temp)
      test$dist.real[[l]] <- as.matrix(dist(test$norm.real[[l]]))
      #Random dataset (Normalize Datasets & Calculate Distances)
      train$norm.rand[[l]] <- randomize(train.temp)
      rownames(train$norm.rand[[l]]) <- rownames(train.temp)
      colnames(train$norm.rand[[l]]) <- colnames(train.temp)
      train$norm.rand[[l]] <- standardNormalization(train$norm.rand[[l]])
      train$dist.rand[[l]] <- as.matrix(dist(train$norm.rand[[l]]))
      test$norm.rand[[l]] <- randomize(test.temp)
      rownames(test$norm.rand[[l]]) <- rownames(test.temp)
      colnames(test$norm.rand[[l]]) <- colnames(test.temp)
      test$norm.rand[[l]] <- standardNormalization(test$norm.rand[[l]])
      test$dist.rand[[l]] <- as.matrix(dist(test$norm.rand[[l]]))
    }
    modelOptimization[[j]] <- modelOpt(train,test,KCheck,alphaCheck,numClusters)
  }
  #############################################
  return(modelOptimization)
}
mainOuterFCV <- function(k,values,numClusters,SelectedValues,KCheck,alphaCheck) {
  ##################################################
  #################### Outer FCV ###################
  train <- list()
  train$norm.real <- list()
  train$dist.real <- list()
  train$norm.rand <- list()
  train$dist.rand <- list()
  test <- list()
  test$norm.real <- list()
  test$dist.real <- list()
  test$norm.rand <- list()
  test$dist.rand <- list()
  for (l in 1:length(values[[k]]$train)) {
    train.temp <- values[[k]]$train[[l]]
    test.temp <- values[[k]]$test[[l]]
    #Real Datasets (Normalize Datasets & Calculate Distances)
    train$norm.real[[l]] <- standardNormalization(train.temp)
    train$dist.real[[l]] <- as.matrix(dist(train$norm.real[[l]]))
    test$norm.real[[l]] <- standardNormalization(test.temp)
    test$dist.real[[l]] <- as.matrix(dist(test$norm.real[[l]]))
    #Random dataset (Normalize Datasets & Calculate Distances)
    train$norm.rand[[l]] <- randomize(train.temp)
    rownames(train$norm.rand[[l]]) <- rownames(train.temp)
    colnames(train$norm.rand[[l]]) <- colnames(train.temp)
    train$norm.rand[[l]] <- standardNormalization(train$norm.rand[[l]])
    train$dist.rand[[l]] <- as.matrix(dist(train$norm.rand[[l]]))
    test$norm.rand[[l]] <- randomize(test.temp)
    rownames(test$norm.rand[[l]]) <- rownames(test.temp)
    colnames(test$norm.rand[[l]]) <- colnames(test.temp)
    test$norm.rand[[l]] <- standardNormalization(test$norm.rand[[l]])
    test$dist.rand[[l]] <- as.matrix(dist(test$norm.rand[[l]]))
  }
  vals <- unlist(strsplit(SelectedValues,'_'))
  OutValidation <- modelOpt(train,test,KCheck[vals[1]],alphaCheck[vals[2]],numClusters)
  #################################################
  return(OutValidation)
}

#SNF Nested Cross-Validation functions
modelOpt <- function(train,test,KCheck,alphaCheck,numClusters) {
  #Run Clustering algorithm with model.train.dataset/model.test.dataset
  Stability <- list()
  Distances <- list()
  for (K in 1:length(KCheck)) {
    Stability[[names(KCheck)[K]]] <- list()
    Distances[[names(KCheck)[K]]] <- list()
    for (alpha in 1:length(alphaCheck)) {
      Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]] <- list()
      Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]] <- list()
      SNF.Train <- SNF.Analysis(train$dist.real,train$dist.rand,KCheck[K],alphaCheck[alpha],numClusters)
      SNF.Test <- SNF.Analysis(test$dist.real,test$dist.rand,KCheck[K],alphaCheck[alpha],numClusters)
      for (l in 2:numClusters) {
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]] <- list()
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]] <- list()
        #################################
        #### Real clusters stability ####
        predict.real <-groupPredict(train$norm.real,test$norm.real,
                                    SNF.Train$Real$Clusters[[l]], K=KCheck[K], alpha=alphaCheck[alpha], t=T, method=1)
        names(predict.real) <- c(rownames(train$norm.real[[1]]),rownames(test$norm.real[[1]]))
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RealPredict <- predict.real
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RealClusters <- SNF.Test$Real$Clusters[[l]]
        #Calculate proportions of recovered labels
        stab.real <- Prop.Stability(SNF.Test$Real$Clusters[[l]],predict.real,t(train$norm.real[[1]]),t(test$norm.real[[1]]),l)
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Real <- stab.real
        max.vals.real <- c()
        for (m in 1:nrow(stab.real)) {max.vals.real <- c(max.vals.real,max(stab.real[m,],na.rm = TRUE))}
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RealPropAVG <- mean(max.vals.real,na.rm = TRUE)
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RealPropMED <- median(max.vals.real,na.rm = TRUE)
        #### Real Distances ####
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Train.Real.EigenVals <- SNF.Train$Real$Metrics$t1[l-1]
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Train.Real.RotationCost <- SNF.Train$Real$Metrics$t2[l-1]
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Test.Real.EigenVals <- SNF.Test$Real$Metrics$t1[l-1]
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Test.Real.RotationCost <- SNF.Test$Real$Metrics$t2[l-1]
        #################################

        #################################
        #### Rand clusters stability ####
        predict.random <-groupPredict(train$norm.rand,test$norm.rand,
                                      SNF.Train$Random$Clusters[[l]], K=KCheck[K], alpha=alphaCheck[alpha], t=T, method=1)
        names(predict.random) <- c(rownames(train$norm.rand[[1]]),rownames(test$norm.rand[[1]]))
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RandomPredict <- predict.random
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RandomClusters <- SNF.Test$Random$Clusters[[l]]
        stab.random <- Prop.Stability(SNF.Test$Random$Clusters[[l]],predict.random,t(train$norm.rand[[1]]),t(test$norm.rand[[1]]),l)
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Random <- stab.random
        max.vals.random <- c()
        for (m in 1:nrow(stab.random)) {max.vals.random <- c(max.vals.random,max(stab.random[m,],na.rm = TRUE))}
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RandomPropAVG <- mean(max.vals.random,na.rm = TRUE)
        Stability[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$RandomPropMED <- median(max.vals.random,na.rm = TRUE)
        #### Random Distances ####
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Train.Rand.EigenVals <- SNF.Train$Random$Metrics$t1[l-1]
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Train.Rand.RotationCost <- SNF.Train$Random$Metrics$t2[l-1]
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Test.Rand.EigenVals <- SNF.Test$Random$Metrics$t1[l-1]
        Distances[[names(KCheck)[K]]][[names(alphaCheck)[alpha]]][[l]]$Test.Rand.RotationCost <- SNF.Test$Random$Metrics$t2[l-1]
        #################################
      }
    }
  }
  return(list(Stability=Stability,Distances=Distances))
}
Prop.Stability <- function(SNFSub2,Stab,GenesSub1,GenesSub2,Clusters) {
  CC.original <- SNFSub2
  names(CC.original) <- colnames(GenesSub2)
  CC.new <- Stab[(ncol(GenesSub1)+1):(ncol(GenesSub1)+ncol(GenesSub2))]
  CC.index.stability<- list(original=list(),new=list())
  for (i in 1:Clusters){
    CC.index.stability$original[[i]] <- names(CC.original)[CC.original==i]
    CC.index.stability$new[[i]] <- names(CC.new)[CC.new==i]
  }
  CC.prop.stability <- data.frame()
  cNames <- c()
  rNames <- c()
  for (i in 1:length(CC.index.stability$original)){
    cNames <- append(cNames,paste0("SNF.Cluster.",i))
    for (j in 1:length(CC.index.stability$new)){
      if (i==1) {rNames <- append(rNames,paste0("Predicted.Cluster.",j))}
      sub1 <- CC.index.stability$original[[i]]
      sub2 <- CC.index.stability$new[[j]]
      comun <- intersect(sub1, sub2)
      total <- unique(c(sub1,sub2))
      val = length(comun) / length(total)
      #val <- cluster_similarity(sub1,sub2,similarity="jaccard")
      CC.prop.stability[i, j] <- val
    }
  }
  colnames(CC.prop.stability) <- cNames
  rownames(CC.prop.stability) <- rNames
  return(CC.prop.stability)
}
SNF.Analysis <- function(values,values.rand,K,alpha,numClusters) {
  require(SNFtool)
  require(NMF)
  ###################
  ### REAL VALUES ###
  #Combine Matrix
  values.aff <- list()
  for (i in 1:length(values)) {values.aff[[i]] <- affinityMatrix(values[[i]],K,alpha)}
  W <- SNF(values.aff,K,20)
  #Get Metrics
  est <- estClustersPlot(W,2:numClusters)
  #Get Assignments
  Clusters <- list()
  for (i in 2:numClusters) {Clusters[[i]] <- spectralClustering(W,i)}
  ###################
  ###################

  #####################
  ### RANDOM VALUES ###
  #Combine Matrix randoms
  values.rand.aff <- list()
  for (i in 1:length(values.rand)) {values.rand.aff[[i]] <- affinityMatrix(values.rand[[i]],K,alpha)}
  W.rand <- SNF(values.rand.aff,K,20)
  #Get Metrics Randoms
  est.rand <- estClustersPlot(W.rand,2:numClusters)
  #Get Assignments randoms
  Clusters.rand <- list()
  for (i in 2:numClusters) {Clusters.rand[[i]] <- spectralClustering(W.rand,i)}
  ######################
  ######################

  return(list(Real = list(Metrics = est,Clusters = Clusters),Random = list(Metrics = est.rand, Clusters = Clusters.rand)))
}
estClustersPlot <- function(W,NUMC) {
  W = (W + t(W))/2
  diag(W) = 0
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    # compute unnormalized Laplacian
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1 / sqrt(degs))
    L = Di %*% L %*% Di
    # compute the eigenvectors corresponding to the k smallest
    # eigs$valuess
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return=T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = eigengap * (1 - eigs$values[1:length(eigs$values) - 1] ) / (1 - eigs$values[2:length(eigs$values)])
    quality = list()
    for (c_index in 1:length(NUMC)) {
      ck = NUMC[c_index]
      UU = eigs$vectors[, 1:ck]
      EigenvectorsDiscrete <- .discretisation(UU)[[1]]
      EigenVectors = EigenvectorsDiscrete^2

      # MATLAB: sort(EigenVectors,2, 'descend');
      temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), function(i) EigenVectors[, i])), ]
      temp1 <- t(apply(temp1, 1, sort, TRUE))

      quality[[c_index]] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) *
        sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*% temp1[, 1:max(2, ck-1)] ))
    }
    #t1 <- sort(eigengap[NUMC], decreasing=TRUE, index.return=T)$ix
    t1 <- eigengap[NUMC]
    # K1 = NUMC[t1[1]]
    # K12 = NUMC[t1[2]]
    #t2 <- sort(unlist(quality), index.return=TRUE)$ix
    t2 <- unlist(quality)
    # K2 <- NUMC[t2[1]]
    # K22 <- NUMC[t2[2]]
  }
  return(list(t1=t1,t2=t2))
}
.discretisation <- function(eigenVectors) {

  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))

  n = nrow(eigenVectors)
  k = ncol(eigenVectors)

  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])

  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }

  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }

  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)

    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]

    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
      break

    lastObjectiveValue = NcutValue
    R = V %*% t(U)

  }

  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}
.discretisationEigenVectorData <- function(eigenVector) {
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1

  return(Y)
}

#SNF nested cross-validation representation
plot.Inner.Corrected <- function(obj,numClusters) {
  require(reshape2)
  require(clusteval)
  require(plotrix)
  require(stringr)

  modelOpt <- list()
  modelOpt.Rand <- list()
  modelOpt.Corrected <- list()
  #ModelOptimization Analysis
  #Get Values
  for (i in 2:numClusters) {modelOpt[[i]] <- list()}
  for (i in 2:numClusters) {modelOpt.Rand[[i]] <- list()}
  for (i in 2:numClusters) {modelOpt.Corrected[[i]] <- list()}

  for (reps in 1:length(obj)) {
    for (i in 1:length(obj[[reps]]$modelOpt)) {
      for (j in 1:length(obj[[reps]]$modelOpt[[i]])) {
        for (k in names(obj[[reps]]$modelOpt[[i]][[j]]$Stability)) {
          for (l in names(obj[[reps]]$modelOpt[[i]][[j]]$Stability[[k]])) {
            for (m in 2:numClusters) {
              params <- paste0(k,'_',l)
              #Props calculate
              modelOpt[[m]][[params]] <- c(modelOpt[[m]][[params]],obj[[reps]]$modelOpt[[i]][[j]]$Stability[[k]][[l]][[m]]$RealPropAVG)
              modelOpt.Rand[[m]][[params]] <- c(modelOpt.Rand[[m]][[params]],obj[[reps]]$modelOpt[[i]][[j]]$Stability[[k]][[l]][[m]]$RandomPropAVG)
              modelOpt.Corrected[[m]][[params]] <- modelOpt[[m]][[params]] - modelOpt.Rand[[m]][[params]]
            }
          }
        }
      }
    }
  }

  Real <- list()
  Rand <- list()
  Corrected <- list()

  #Calculate Summary
  Real$AVG <- matrix(0,numClusters-1,length(modelOpt[[2]]))
  rownames(Real$AVG) <- c(2:numClusters)
  colnames(Real$AVG) <- names(modelOpt[[2]])
  Real$AVG <- as.data.frame(Real$AVG)
  Real$SD <- matrix(0,numClusters-1,length(modelOpt[[2]]))
  rownames(Real$SD) <- c(2:numClusters)
  colnames(Real$SD) <- names(modelOpt[[2]])
  Real$SD <- as.data.frame(Real$SD)

  Rand$AVG <- matrix(0,numClusters-1,length(modelOpt.Rand[[2]]))
  rownames(Rand$AVG) <- c(2:numClusters)
  colnames(Rand$AVG) <- names(modelOpt.Rand[[2]])
  Rand$AVG <- as.data.frame(Rand$AVG)
  Rand$SD <- matrix(0,numClusters-1,length(modelOpt.Rand[[2]]))
  rownames(Rand$SD) <- c(2:numClusters)
  colnames(Rand$SD) <- names(modelOpt.Rand[[2]])
  Rand$SD <- as.data.frame(Rand$SD)

  Corrected$AVG <- matrix(0,numClusters-1,length(modelOpt.Corrected[[2]]))
  rownames(Corrected$AVG) <- c(2:numClusters)
  colnames(Corrected$AVG) <- names(modelOpt.Corrected[[2]])
  Corrected$AVG <- as.data.frame(Corrected$AVG)
  Corrected$SD <- matrix(0,numClusters-1,length(modelOpt.Corrected[[2]]))
  rownames(Corrected$SD) <- c(2:numClusters)
  colnames(Corrected$SD) <- names(modelOpt.Corrected[[2]])
  Corrected$SD <- as.data.frame(Corrected$SD)

  for (i in 2:numClusters) {
    for (j in names(modelOpt[[i]])) {
      Real$AVG[i-1,j] <- mean(modelOpt[[i]][[j]],na.rm = TRUE)
      Real$SD[i-1,j] <- std.error(modelOpt[[i]][[j]],na.rm = TRUE)
      Rand$AVG[i-1,j] <- mean(modelOpt.Rand[[i]][[j]],na.rm = TRUE)
      Rand$SD[i-1,j] <- std.error(modelOpt.Rand[[i]][[j]],na.rm = TRUE)
      Corrected$AVG[i-1,j] <- mean(modelOpt.Corrected[[i]][[j]],na.rm = TRUE)
      Corrected$SD[i-1,j] <- std.error(modelOpt.Corrected[[i]][[j]],na.rm = TRUE)
    }
  }

  RealVals <- cbind(melt(as.matrix(Real$AVG)),melt(as.matrix(Real$SD)))
  RealVals$type <- "Real"
  RealVals$pos <- "Datasets"
  RealVals <- cbind(RealVals,str_split_fixed(RealVals[,2], "_", 2))
  RealVals <- RealVals[,c(1,9,10,3,6,7,8)]
  colnames(RealVals) <- c("Clusters","Neighbors","Hyperparameter","avg","sd","type","pos")

  RandVals <- cbind(melt(as.matrix(Rand$AVG)),melt(as.matrix(Rand$SD)))
  RandVals$type <- "Random"
  RandVals$pos <- "Datasets"
  RandVals <- cbind(RandVals,str_split_fixed(RandVals[,2], "_", 2))
  RandVals <- RandVals[,c(1,9,10,3,6,7,8)]
  colnames(RandVals) <- c("Clusters","Neighbors","Hyperparameter","avg","sd","type","pos")

  CorrectedVals <- cbind(melt(as.matrix(Corrected$AVG)),melt(as.matrix(Corrected$SD)))
  CorrectedVals$type <- "Corrected"
  CorrectedVals$pos <- "Corrected"
  CorrectedVals <- cbind(CorrectedVals,str_split_fixed(CorrectedVals[,2], "_", 2))
  CorrectedVals <- CorrectedVals[,c(1,9,10,3,6,7,8)]
  colnames(CorrectedVals) <- c("Clusters","Neighbors","Hyperparameter","avg","sd","type","pos")

  Vals <- rbind(RealVals,RandVals,CorrectedVals)

  return(Vals)
}
plot.Outer.Corrected <- function(obj,numClusters) {
  require(reshape2)
  require(clusteval)
  require(plotrix)
  require(stringr)

  prop <- list()
  prop.Rand <- list()
  prop.Corrected <- list()
  # jaccard <- list()
  # jaccard.Rand <- list()
  # rand <- list()
  # rand.Rand <- list()
  #ModelOptimization Analysis
  for (i in 2:numClusters) {prop[[i]] <- list()}
  for (i in 2:numClusters) {prop.Rand[[i]] <- list()}
  for (i in 2:numClusters) {prop.Corrected[[i]] <- list()}

  #Get Values
  params <- unlist(strsplit(obj[[1]]$selectValues,'_'))
  for (i in 1:length(obj[[1]]$outValidation)) {
    for (j in 1:length(obj[[1]]$outValidation[[i]]$Stability)) {
      for (m in 2:numClusters) {
        #Props calculate
        prop[[m]] <- c(prop[[m]],obj[[1]]$outValidation[[i]]$Stability[[params[[1]]]][[params[[2]]]][[m]]$RealPropAVG)
        prop.Rand[[m]] <- c(prop.Rand[[m]],obj[[1]]$outValidation[[i]]$Stability[[params[[1]]]][[params[[2]]]][[m]]$RandomPropAVG)
        prop.Corrected[[m]] <- as.numeric(prop[[m]]) - as.numeric(prop.Rand[[m]])
      }
    }
  }
  Results <- matrix(NA,ncol = 6)
  for (i in 2:numClusters) {
    Results <- rbind(Results,cbind(mean(unlist(prop[[i]])),std.error(unlist(prop[[i]])),paste0(params[1],"_",params[2]),"Real",i,"Datasets"))
    Results <- rbind(Results,cbind(mean(unlist(prop.Rand[[i]])),std.error(unlist(prop.Rand[[i]])),paste0(params[1],"_",params[2]),"Random",i,"Datasets"))
    Results <- rbind(Results,cbind(mean(unlist(prop.Corrected[[i]])),std.error(unlist(prop.Corrected[[i]])),paste0(params[1],"_",params[2]),"Corrected",i,"Corrected"))
  }
  Results <- Results[-1,]
  colnames(Results) <- c("avg","sd","Params","type","Clusters","pos")
  return(as.data.frame(Results))
}
