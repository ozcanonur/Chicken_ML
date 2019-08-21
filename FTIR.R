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
library(extrafont)
library(extrafontdb)
library(MLmetrics)
options(scipen=999)
#----Read and pre-process----
rawdata <- readFTIR(outliersRemoved = F, batch = 0)
b1 <- readFTIR(outliersRemoved = T, batch = 1, waves = c(850:1700, 3200:3400))
b2 <- readFTIR(outliersRemoved = F, batch = 2, waves = c(850:1700, 3200:3400))

# Apply baseline correction, binning, normalisation, center and scale
preprocessed <- apply_paper_preprocess(datalist = list(b1,b2), bin = 4, normalarea = F, centerandscale = T)
b1 <- preprocessed[[1]]
b2 <- preprocessed[[2]]
#----Cook's distance test----
cooksdistancetest(b1) # read the data without outliers removed and preprocesses before
#----(NOT USED IN PAPER)----
# data <- keepwave(rawdata, filterlist)
# data <- findandkeepwavesRF(data, 500, 0) # Feature selection by RF, 2nd param topx, 3rd param threshold
#----(NOT USED IN PAPER)----
# # Find outliers by spectra.outliers
# outlier <- spectra.outliers("forspectraoutliersPREPPED_B2.csv", plots=FALSE,k=10)
# outlier <- spectra.outliers("forspectraoutliersPREPPED_B1.csv", plots=FALSE,k=10)
#----Single iteration----
# Split if going for train/test split
# 0 = TVC, 1 = time, 2 = temp, 3 = TVC, 4=Kennard-Stone
return <- split(rawdata, 0) 
train <- return$train
test <- return$test
# Train models
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
# Set CV method
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats=3
                           ,savePredictions = "final"
                           ,allowParallel = T)
gridlist <- gridList() # Get gridlist

cl <- makePSOCKcluster(3) # Input cores to use for parallel processing
registerDoParallel(cl) # Start parallel processing

# Test all the models listed with their respective grid lists in CV
# Change data parameter of caretList to change training data
models <-  caretList(TVC~., data=b1, trControl = fitControl, metric = "RMSE",
                     tuneList = list(
                      pls=caretModelSpec(method="pls",tuneGrid=gridlist$plsgrid)
                      ,ridge=caretModelSpec(method="ridge",tuneGrid=gridlist$ridgegrid)
                      ,svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                      ,svmRadial=caretModelSpec(method="svmRadial", tuneGrid=gridlist$svmRadialgrid)
                      ,rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid, ntree=500)
                      ,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                      ,pcr=caretModelSpec(method="pcr",tuneGrid=gridlist$plsgrid)
                      ,lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid, use.Gram=F)
                      ,nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, 
                                           linout = TRUE, maxit=1000, MaxNWts=7000)
                      )
)

stopCluster(cl) # Stop parallel processing
registerDoSEQ()

# Predict test data for all models
multiplemodelpredict(models = models, test = b2, plot = T) # Results on set plotted

# Predict test data for a selected model from caretList and plot predicted vs actual
plot <- custompredict(model = models$nnet, test = b2)

#----(NOT USED IN PAPER)----
# # Print CV resample results & boxplot for model comparison
# results <- resamples(models) # Check CV results
# summary(results, metric="RMSE")
# print(bwplot(results,metric = "RMSE", main="RMSE over 25 resamples"))
#----(NOT USED IN PAPER)----
# # Wave vs Intensity for fresh/spoiled means
# dataSub <- meanintensitygraph(data = data, low = 4, high = 6) # low and high are thresholds to split data
#----PCA stuff----
# Plot PCA
plot <- doPca(data = b1, ncomp = c(1,2))
print(plot)
#----(NOT USED IN PAPER)----
# pca <- PCA(data[,-which(colnames(data) %in% "TVC")], graph = F)
# # Check eigenvalues
# eig.val <- get_eigenvalue(pca) 
# # Check contribution to components by each feature
# contribution <- as.data.frame(get_pca_var(pca)$cos2)
# contribution$sum <- rowSums(contribution[,1:2]) # contrib to PC1 and PC2
# 
# filterwavesbyPCAcontrib(data, contribution, 0.6) # Filter by pca contribution
# 
# # Check contribution to components by each sample
# contrib <- (fviz_contrib(pca, choice = "ind", axes = 1:2))$data # See contributions to pca 1 and 2 by samples
# contrib <- contrib[order(-contrib$contrib),] # Sort by highest contribution
# 
# # PCA communality (Literature suggest to drop below 0.6)
# pca2 <- principal(rawdata[,-which(colnames(rawdata) %in% "TVC")], nfactors = 4)
# View(pca2$communality)
# 
# # Filter waves with communality value under threshold
# filterlist <- rownames(as.data.frame(pca2$communality)[pca2$communality > 0.6,,drop=F])
# data <- data[, filterlist]
#----(NOT USED IN PAPER)----
# # STD things
# std <- plotSTD(rawdata) # Plot SD over wave lengths
# data <- filterwavesfromSTD(data, std, 0.1) # Filter waves by SD threshold
# # Check lm regression coef
# coef <- lmcoefficient(data)
#----Compare waveIntensity in two selected samples----
# plot <- ftir_compare_two_samples(data = data, sample_one = "0C_86h_A.csv",
#                                  sample_two = "0C_86h_B.csv")
# print(plot)
#----(NOT USED IN PAPER)----
# # Find linear dependencies and remove them
# comboInfo <- findLinearCombos(data[,-which(colnames(data) %in% "TVC")])
# data <- data[,-comboInfo$remove]
#----(NOT USED IN PAPER)----
# # Correlation check
# corr <- cor(data)
# highlyCorDescr <- findCorrelation(corr, cutoff = 0.99, names = T, exact=T)
# 
# data <- data[!colnames(data) %in% highlyCorDescr]
#----Multiple iterations----
# Run iterations with different train and tests with 8 models
iters <- 100
results <- multipleiterations(rawdata, iters, parallelcores = 3) 

rmseResults <- results$rmseResults
maeResults <- results$maeResults
r2Results <- results$r2Results
accResults <- results$accResults

# Combine the results
finalresults <- getfinalresults(rmseResults, maeResults, r2Results, accResults)
View(finalresults)
# Add finalresults to resultlist (list(finalresults1, finalresults2...)) from different workspaces
# for finalheatmap() that will happen later on
#----(NOT USED IN PAPER)----
# # Box plot for all rmse 100 iter
# boxplot(t(rmseResults[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(11), 
#         main="RMSE over 100 iterations", cex.axis=0.7)
# 
# # Box plot for all accuracy 100 iter
# boxplot(t(accResults[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(11), 
#         main="Accuracy over 100 iterations", cex.axis=0.7)
# 
# # Best results by RMSE and Accuracy
# rmseResults <- rmseResults[order(rownames(rmseResults)),] # Order by names to get in proper format
# accResults <- accResults[order(rownames(accResults)),]
#----Plot heatmap----
plot <- finalheatmap(resultlist = resultlist, names = c("Batch 1", "Batch 2", "Combined"), limits = c(65, 90))
#----Variable importance by iters----
imp <- variableimportancetest(data = b1, iters = 2, model = "pls", parallelcores = 3)
iter_result<- imp$iterresult
mean_importance <- imp$meanimp
#----(NOT USED IN PAPER)----
# # impmean <- transform(nnetimp, wave=as.numeric(stringr::str_extract(rownames(nnetimp),"[0-9]+\\.+[0-9]*")))
# # impmean <- impmean[order(-impmean$Overall), ]
# # impmean$wave <- round(impmean$wave,0)
# 
# rawdata <- x
# filterlist <- round(impmeannnet[1:200,]$wave,0) # nnet 50 works best
# rawdata <- keepwave(rawdata, filterlist)
# b2 <- rawdata[98:nrow(rawdata),] # Batch 2
# b1 <- rawdata[1:97,] # Batch 1
# 
# x <- model.matrix(TVC~.,rawdata)[,-1]
# y <- rawdata$TVC
# model <- cv.glmnet(x=as.matrix(x),y,alpha=1)
# 
# CF <- as.matrix(coef(model, model$lambda.1se))
# CF <- CF[CF!=0,,drop=F]
# CF <- CF[-1,,drop=F]
# CF <- transform(CF, wave=as.numeric(stringr::str_extract(rownames(CF),"[0-9]+\\.+[0-9]*")))
# CF$wave <- round(CF$wave,0)
# filterlist <- round(CF$wave,0) # nnet 50 works best
#----(NOT USED IN PAPER)----
# # MCUVE-PLS
# pls2 <- mcuve_pls(data$TVC, as.matrix(data[-which(colnames(data) %in% "TVC")]), ncomp = 10, N = 20)
#----(NOT USED IN PAPER)----
# # Bootstrap
# # Bootstrap 95% CI for R-Squared
# library(boot)
# # function to obtain R-Squared from the data 
# rsq <- function(formula, data, indices) {
#   d <- data[indices,] # allows boot to select sample 
#   fit <- pls(formula, data=d, ncomp=2)
#   return(summary(fit)$r.square)
# } 
# # bootstrapping with 1000 replications 
# results <- boot(data=rawdata, statistic=rsq, 
#                 R=10, formula=TVC~.)
# 
# # view results
# results 
# plot(results)
# 
# # get 95% confidence interval 
# boot.ci(results, type="bca")
#----(NOT USED IN PAPER)----
# # Stepwise exhaustive
# leap <- leaps::leaps(as.matrix(data[,-which(colnames(data) %in% "TVC")]),
#               data$TVC,names=colnames(data[,-which(colnames(data) %in% "TVC")]),method = "adjr2")
# View(leap$adjr2)
# 
# leap <- leaps::regsubsets(x=as.matrix(data[,-which(colnames(data) %in% "TVC")]),
#                      y=data$TVC,method = "exhaustive", really.big=T, all.best=T)
#----Save plot with anti-aliasing (prettier)----
ggsave(plot = plot, "plot.png",type = "cairo-png")