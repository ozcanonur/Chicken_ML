#----Imports----
library(caret)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(randomForest)
library(rms)
library(e1071)
library(caretEnsemble)
library(doParallel)
library(factoextra)
library(FactoMineR)
library(tidyverse)
library("viridis")
library(corrplot)
source("Util.R")
library(prospectr)
library(soil.spec)
library(psych)
library(hyperSpec)
options(scipen=999)
#----Read and pre-process----
rawdata <- read.table("FTIR_R1-R3_RemovedOutliers.csv",sep = ",", header = TRUE)

rawdata <- rawdata[1:98,] # Batch 1
rawdata <- rawdata[99:nrow(rawdata),] # Batch 2

rownames(rawdata) <- rawdata[,1]
rawdata <- rawdata[,-1]
TVC <- rawdata$TVC

data <- rawdata

data <- keepwave(data, 750:ncol(data)) # Keep wanted wave lengths

data <- applypreprocess(data, "baseline") # Baseline correction
data <- applypreprocess(data, "snv") # Standard normal derivate
data <- applypreprocess(data, "derivation", deriv = 2) # Standard normal derivate

data$TVC <- TVC # Attach TVC back

# More pre-processing options
# data <- prep(rawdata, c("YeoJohnson")) # RMSE: 1.7
# data <- prep(rawdata, c("nzv")) # RMSE: 1.19
# data <- prep(rawdata, c("corr")) # RMSE: 1.63, leaves only 2 features ??
# data <- prep(rawdata, c("center", "scale", "nzv")) # RMSE: TERRIBLE 3.35
# data <- prep(rawdata, c("expoTrans")) # RMSE: TERRIBLE 13.55

# data <- findandkeepwavesRF(data, 500, 0) # Feature selection by RF, 2nd param topx, 3rd param threshold
#----Find outliers by spectra.outliers----
outlier <- spectra.outliers("C:/Users/Onur/Desktop/THESIS/DATA/forspectraoutliers_B1.csv")
outlierlist <- outlier[["all.outliers"]]
data <- data[-outlierlist, ]
#----Split to train and test----
# 0 = naive, 1 = time, 2 = temp, 3 = TVC, 4=Kennard-Stone
return <- split(data, 4) 
train <- return$train
test <- return$test

# Get good train
train <- data[which(rownames(data) %in% rownames(ACCbestnnet$trainingData)), ]
test <- data[which(!rownames(data) %in% rownames(ACCbestnnet$trainingData)), ]

#----Train models----
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
cl <- makePSOCKcluster(3)
registerDoParallel(cl)

gridlist <- gridListFTIR(data)
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats=3
                           ,savePredictions = "final"
                           #,index = createResample(data$TVC, 25)
                           ,allowParallel = T)
plsgrid <- expand.grid(ncomp=1:15)
models <-  caretEnsemble::caretList(TVC~., data=data, trControl = fitControl, metric = "RMSE"
                                    ,preProcess = c("center","scale")
                                    ,tuneList = list(
                                      #pls=caretModelSpec(method="pls", tuneGrid=plsgrid)
                                      #svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                                      # ,svmRadial=caretModelSpec(method="svmRadial", tuneGrid=gridlist$svmRadialgrid)
                                      # ,rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid, ntree=1000)
                                      # ,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                                       #,pcr=caretModelSpec(method="pcr", tuneGrid=gridlist$plsgrid)
                                      # ,ridge=caretModelSpec(method="ridge")
                                      lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid, use.Gram=F)
                                      #nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, linout = TRUE, maxit=1000)
                                    )
)

stopCluster(cl)
registerDoSEQ()

custompredict(models$nnet, test) # Predict a single model and plot

plotvariableimportance(models$pls) # Plot variable importance for a single model

results <- resamples(models) # Check CV results
summary(results, metric="RMSE")
print(bwplot(results,metric = "RMSE", main="RMSE over 25 resamples"))

testresults <- multiplemodelpredict(models) # Results on set plotted
testresults

#----Wave vs Intensity for fresh/spoiled means----
dataSub <- meanintensitygraph(data, 4, 6) # 
#----PCA stuff----
# Plot PCA
pca <- PCA(data[,-which(colnames(data) %in% "TVC")], graph=F)
plot(pca)

# Check eigenvalues
eig.val <- get_eigenvalue(pca) 

# Check contribution to components by each feature
contribution <- as.data.frame(get_pca_var(pca)$cos2)
contribution$sum <- rowSums(contribution[,1:2]) # contrib to PC1 and PC2

filterwavesbyPCAcontrib(data, contribution, 0.6) # Filter by pca contribution

# Check contribution to components by each sample
contrib <- (fviz_contrib(pca, choice = "ind", axes = 1:2))$data # See contributions to pca 1 and 2 by samples
contrib <- contrib[order(-contrib$contrib),] # Sort by highest contribution

# PCA plot with ellipses and low-med-high seperation
quantiles <- quantile(data$TVC, c(.33, .66)) # Find quantiles for TVC values
data$TVCinterval <- ifelse(data$TVC <= quantiles[1], 'low',
                           ifelse(data$TVC > quantiles[1] & data$TVC <= quantiles[2], 'med',
                                  ifelse(data$TVC > quantiles[2], 'high', '')))
ggbiplot::ggbiplot(pca, choices = c(1,2), var.axes = FALSE,
                         groups = data$TVCinterval, ellipse = T, labels = rownames(data))

data <- data[,!colnames(data) %in% "TVCinterval"]

# PCA communality (Literature suggest to drop below 0.6)
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
pca2 <- principal(data[,-which(colnames(data) %in% "TVC")], nfactors = 7)
stopCluster(cl)
registerDoSEQ()
View(pca2$communality)

# Filter waves with communality value under threshold
filterlist <- rownames(as.data.frame(pca2$communality)[pca2$communality > 0.6,,drop=F])
data <- data[, filterlist]

#----STD things----
std <- plotSTD(data) # Plot SD over wave lengths
data <- filterwavesfromSTD(data, std, 0.1) # Filter waves by SD threshold
#----Check lm regression coef----
coef <- lmcoefficient(data)
#----Compare waveIntensity in two selected samples----
par(mfrow=c(1,2))
compareIntensity(data, c("0C_86h_A.csv", "0C_86h_B.csv"))
compareIntensity(data, c("15C_14h_A.csv", "15C_14h_B.csv"))
#----Find linear dependencies and remove them----
comboInfo <- findLinearCombos(data[,-1])
data <- data[,-comboInfo$remove]
#----Correlation check----
corr <- cor(data[,-which(colnames(data) %in% "TVC")])
highlyCorDescr <- findCorrelation(corr, cutoff = .99) # Find correlation above threshold
corr <- as.data.frame(corr[,-highlyCorDescr]) # Filter out correlated waves

data <- data[colnames(data) %in% colnames(corr)]
data$TVC <- rawdata$TVC # Add TVC back
#----Multiple iterations----
# Run iterations with different train and tests with 8 models
iters <- 25
results <- multipleiterations(data, iters) 

rmseResults <- results$rmseResults
predAcc <- results$predAcc
RMSEbestlist <- results$RMSEbestlist
ACCbestlist <- results$ACCbestlist

# RESULTS
rmseResults$min <- apply(rmseResults,1,min) # Min as extra col
rmseResults$mean <- apply(rmseResults,1,FUN=mean) # Median as extra col
rmseResults$max <- apply(rmseResults,1,max) # Max as extra col
rmseResults <- rmseResults[order(rmseResults$mean),] # Order by median

predAcc$min <- apply(predAcc,1,min) # Same as above
predAcc$mean <- apply(predAcc,1,FUN=mean)
predAcc$max <- apply(predAcc,1,max)
predAcc <- predAcc[order(-predAcc$mean),]

# Box plot for all rmse 100 iter
boxplot(t(rmseResults[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(11), 
        main="RMSE over 100 iterations", cex.axis=0.7)

# Box plot for all accuracy 100 iter
boxplot(t(predAcc[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(11), 
        main="Accuracy over 100 iterations", cex.axis=0.7)

# Best results by RMSE and Accuracy
rmseResults <- rmseResults[order(rownames(rmseResults)),] # Order by names to get in proper format
predAcc <- predAcc[order(rownames(predAcc)),]

finalresults <- data.frame(matrix(nrow=8,ncol=2)) # Create results DF
rownames(finalresults) <- rownames(predAcc)
colnames(finalresults) <- c("RMSE", "Accuracy")
finalresults[,1] <-  round(rmseResults$min,2) # Minimum RMSE for models
finalresults[,2] <- round(predAcc$max,2) # Max accuracy for models
finalresults$name <- rownames(finalresults)

#----Plot heatmap----
finalresultsB1 <- finalresultsB1[order(rownames(finalresultsB1)),][,-3]
finalresultsB2 <- finalresultsB2[order(rownames(finalresultsB2)),][,-3]
finalresultsB1andB2 <- finalresultsB1andB2[order(rownames(finalresultsB1andB2)),][,-3]

finalresultsB1RMSE <- finalresultsB1[,1, drop=F]
finalresultsB1RMSE$model <- rownames(finalresultsB1)
finalresultsB1RMSE <- finalresultsB1RMSE[,-1,drop=F]
finalresultsB1RMSE$RMSE <- finalresultsB1[,1]
finalresultsB1RMSE$batch <- rep("Batch 1",8)

finalresultsB2RMSE <- finalresultsB2[,1, drop=F]
finalresultsB2RMSE$model <- rownames(finalresultsB2)
finalresultsB2RMSE <- finalresultsB2RMSE[,-1,drop=F]
finalresultsB2RMSE$RMSE <- finalresultsB2[,1]
finalresultsB2RMSE$batch <- rep("Batch 2",8)

finalresultsB1andB2RMSE <- finalresultsB1andB2[,1, drop=F]
finalresultsB1andB2RMSE$model <- rownames(finalresultsB1andB2)
finalresultsB1andB2RMSE <- finalresultsB1andB2RMSE[,-1,drop=F]
finalresultsB1andB2RMSE$RMSE <- finalresultsB1andB2[,1]
finalresultsB1andB2RMSE$batch <- rep("All data",8)

finalresultsB1ACC <- finalresultsB1[,2, drop=F]
finalresultsB1ACC$model <- rownames(finalresultsB1)
finalresultsB1ACC <- finalresultsB1ACC[,-1,drop=F]
finalresultsB1ACC$Acc <- finalresultsB1[,2]
finalresultsB1ACC$batch <- rep("Batch 1",8)

finalresultsB2ACC <- finalresultsB2[,2, drop=F]
finalresultsB2ACC$model <- rownames(finalresultsB2)
finalresultsB2ACC <- finalresultsB2ACC[,-1,drop=F]
finalresultsB2ACC$Acc <- finalresultsB2[,2]
finalresultsB2ACC$batch <- rep("Batch 2",8)

finalresultsB1andB2ACC <- finalresultsB1andB2[,2, drop=F]
finalresultsB1andB2ACC$model <- rownames(finalresultsB1andB2)
finalresultsB1andB2ACC <- finalresultsB1andB2ACC[,-1,drop=F]
finalresultsB1andB2ACC$Acc <- finalresultsB1andB2[,2]
finalresultsB1andB2ACC$batch <- rep("All data", 8)

forheatmapRMSE <- rbind(finalresultsB1RMSE, finalresultsB2RMSE, finalresultsB1andB2RMSE)

forheatmapAcc <- rbind(finalresultsB1ACC, finalresultsB2ACC, finalresultsB1andB2ACC)
forheatmapAcc$Acc <- forheatmapAcc$Acc * 100

forheatmapRMSE <- forheatmapRMSE[which(!forheatmapRMSE$model %in% c("knn","rf")),]
forheatmapAcc <- forheatmapAcc[which(!forheatmapAcc$model %in% c("knn","rf")),]

ggplot(forheatmapAcc, aes(batch, model)) +
  geom_tile(aes(fill=Acc)) +
  scale_fill_gradientn(colors=colorRampPalette(c("red","green"))(3)) +
  theme_bw(base_size = 20) +
  labs(x="", y="") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  geom_text(data=forheatmapAcc,aes(batch,model, label=paste("Accuracy:",Acc,"%")), vjust=-1, size=5) +
  geom_text(data=forheatmapRMSE,aes(batch,model,label=paste("RMSE:",RMSE)),vjust=2, size=5) +
  theme(legend.position = "none")

#----Time pred----
rawdata <- read.table("FTIR_R1-R3_1000to1700_Fulloutliers.csv" ,sep = ",", header = TRUE)
rawdata <- rawdata[1:116,] # Batch 1
#rawdata <- rawdata[96:nrow(rawdata),] # Batch 2

rownames <- as.character(rawdata[,1])
rawdata <- rawdata[,-1]
rawdata <- rawdata[,-1]
data <- rawdata

spc <- createspc(data)
bl <- spc.fit.poly.below(spc)
spc <- spc - bl
data <- as.data.frame(spc$spc) # Xfer spc to data
data <- as.data.frame(standardNormalVariate(data)) 

data$Sample <- rownames
data <- transform(data, Temperature=stringr::str_extract(Sample,"[0-9]*((?=C)|(?=B)|(?=A))"))
data <- data[,-which(colnames(data) %in% "Sample"),]
rownames(data) <- rownames

return <- split(data, 0) 
train <- return$train
test <- return$test

train_index = caret::createDataPartition(y=data$Temperature, p=0.7, list=FALSE, times=1)
train <- data[train_index,]
test <- data[-train_index,]

cl <- makePSOCKcluster(3)
registerDoParallel(cl)
gridlist <- gridListFTIR(train)
plsgrid <- expand.grid(.ncomp = seq(from=1,to=20,by=1))
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                           savePredictions = "final"
                           #,preProcOptions = list(thresh=0.95)
                           ,index = createResample(train$TVC, 25)
                           ,allowParallel = T)
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
models <-  caretEnsemble::caretList(Temperature~., data=train, trControl = fitControl
                                    #,preProcess = c("center","scale")
                                    ,tuneList = list(
                                      #pls=caretModelSpec(method="pls", tuneGrid=plsgrid)
                                      #,svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                                      #,svmRadial=caretModelSpec(method="svmRadial", tuneGrid=gridlist$svmRadialgrid)
                                      rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid)
                                      #,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                                      #,pcr=caretModelSpec(method="pcr", tuneGrid=gridlist$plsgrid)
                                      # lm=caretModelSpec(method="lm")
                                      #,lmStepAIC=caretModelSpec(method="lmStepAIC")
                                      #,ridge=caretModelSpec(method="ridge")
                                      #,lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid, use.Gram=F)
                                      #lasso=caretModelSpec(method="lasso")
                                      #nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, linout = TRUE)
                                    )
)
stopCluster(cl)
registerDoSEQ()

print(models$rf)
pred <- predict(models$pls,test)



#----Josh's thing----
rawdata <- read.csv("Copy of VM_Air.csv",header = T,sep = ",")
bacterial <- read.csv("micro_Air.csv", header=T,sep=",")[,c(1,4)]
merged <- merge(rawdata,bacterial,by="Sample")
colnames(merged)[ncol(merged)] <- "TVC"
merged[,"TVC"] <- log10(merged[,"TVC"])
rownames(merged) <- merged[,1]
merged <- merged[,-1]

return <- split(merged, 0) 
train <- return$train
test <- return$test

set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
cl <- makePSOCKcluster(3)
registerDoParallel(cl)

gridlist <- gridListVM(merged)
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats=3
                           ,savePredictions = "final"
                           #,preProcOptions = list(thresh=0.95)
                           ,index = createResample(train$TVC, 25)
                           ,allowParallel = T)
models <-  caretEnsemble::caretList(TVC~., data=train, trControl = fitControl, metric = "RMSE"
                                    ,preProcess = c("center","scale")
                                    ,tuneList = list(
                                      # pls=caretModelSpec(method="pls", tuneGrid=gridlist$plsgrid)
                                      # ,svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                                      # ,svmRadial=caretModelSpec(method="svmRadial", tuneGrid=gridlist$svmRadialgrid)
                                      rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid, ntree=500)
                                      # ,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                                      # ,pcr=caretModelSpec(method="pcr", tuneGrid=gridlist$plsgrid)
                                      #  ,ridge=caretModelSpec(method="ridge")
                                      # ,lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid, use.Gram=F)
                                      ,nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, linout = TRUE)
                                    )
)

stopCluster(cl)
registerDoSEQ()

results <- resamples(models) # Check CV results
summary(results, metric="RMSE")

testresults <- multiplemodelpredict(models) # Results on set plotted
testresults