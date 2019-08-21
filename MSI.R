#----Imports----
library(caret)
library(plyr)
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
library(scales)
options(scipen=999)
#----Read----
rawdata <-  readMSI(outliersRemoved =  F, batch = 0, centerandscale = F)
b1 <- readMSI(outliersRemoved = T, batch = 1, centerandscale = T)
b2 <- readMSI(outliersRemoved = F, batch = 2, centerandscale = T)
#----Outlier test----
plot <- cooksdistancetest(b1)
print(plot)
#----Single iteration----
# Split if going for train/test split
# 0 = TVC, 1 = time, 2 = temp, 3 = TVC, 4=Ken-stone
return <- split(data = rawdata, option = 0)

train <- return$train
test <- return$test
# Train models
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
# Set CV method
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3
                           ,savePredictions = "final"
                           ,allowParallel = T
                           )
gridlist <- gridList() # Get gridlist
gridlist$plsgrid <- expand.grid(ncomp=1:14) # Change pls grid (MSI data can only go up to ~15)

cl <- makePSOCKcluster(3) # Input cores to use for parallel procesisng
registerDoParallel(cl) # Start parallel processing

# Test all the models listed with their respective grid lists in CV
# Change data parameter of caretList to change training data
models <-  caretEnsemble::caretList(TVC~., data=b1, trControl = fitControl, metric = "RMSE"
                                    ,tuneList = list(
                                      pls=caretModelSpec(method="pls",tuneGrid=gridlist$plsgrid)
                                      ,ridge=caretModelSpec(method="ridge",tuneGrid=gridlist$ridgegrid)
                                      ,svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                                      ,svmRadial=caretModelSpec(method="svmRadial", tuneGrid=gridlist$svmRadialgrid)
                                      ,rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid, ntree=500)
                                      ,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                                      ,pcr=caretModelSpec(method="pcr",tuneGrid=gridlist$plsgrid)
                                      ,lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid, use.Gram=F)
                                      ,nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, linout = TRUE, maxit=1000, MaxNWts=7000)
                                    )
)
stopCluster(cl) # Stop parallel processing
registerDoSEQ()

# Predict test data for all models
multiplemodelpredict(models = models, test = b2, plot = F)

# Predict test data for a selected model from caretList and plot predicted vs actual
custompredict(model = models$rf, test = b2)

#----(NOT USED IN PAPER)----
# # Print CV resample results & boxplot for model comparison
# results <- resamples(models) # RMSE results over 25 resamples from CV phase
# summary(results, metric = "RMSE")
# bwplot(results,metric = "RMSE", main="RMSE over 25 resamples")
#----Multiple iterations----
iters <- 100
results <- multipleiterations(data = rawdata, iters = iters)

# Get performance metrics from the iteration
rmseResults <- results$rmseResults
maeResults <- results$maeResults
r2Results <- results$r2Results
accResults <- results$accResults

# Combine the results
finalresults <- getfinalresults(rmseResults, maeResults, r2Results, accResults)
View(finalresults) 
# Add finalresults to resultlist (list(finalresults1, finalresults2...)) from different workspaces
# for heatmap that will happen later on
#----(NOT USED IN PAPER)----
# Box plot for all rmse iters
# boxplot(t(rmseResults[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(8), 
#         main=paste("RMSE over", iters,"iterations"), cex.axis=0.7)
# 
# # Box plot for all accuracy iters
# boxplot(t(predAcc[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(8), 
#         main=paste("Accuracy over", iters,"iterations"), cex.axis=0.7)
# 
# # Measuring how many times a sample has appeared in best models' training sets
# # Extract time/temp and AorB as new columns
# countbest <- transform(countbest, Temperature=stringr::str_extract(name,"[0-9]*((?=C)|(?=B)|(?=A))"))
# countbest <- transform(countbest, AorB=stringr::str_extract(name,"[A-B]"))
# countbest <- transform(countbest, Time=as.numeric(stringr::str_extract(name,"[0-9]*(?=h)")))
# 
# # Count how many times a storage time has appeared
# timecount <- aggregate(countbest$Count, by=list(Time=countbest$Time), FUN=sum)
# colnames(timecount)[2] <- "Count"
# 
# # Count how many times a storage temperature has appeared
# tempcount <- aggregate(countbest$Count, by=list(Temperature=countbest$Temperature), FUN=sum)
# colnames(tempcount)[2] <- "Count"
# 
# # Count how many times A or B has appeared
# abcount <- aggregate(countbest$Count, by=list(AorB=countbest$AorB), FUN=sum)
# colnames(abcount)[2] <- "Count"

#----(NOT USED IN PAPER)----
# # Calculate correlation between features
# corr <- cor(rawdata)
# highlyCorDescr <- findCorrelation(corr, cutoff = .99) # Find correlation above threshold
# corr <- as.data.frame(corr[,-highlyCorDescr]) # Filter out correlated waves
# 
# rawdata <- rawdata[colnames(rawdata) %in% colnames(corr)] # Filter features that are correlated
#----(NOT USED IN PAPER)----
# # Ensemble
# as.data.frame(modelCor(resamples(models))) # Check correlation between models
# 
# # Use models from caretlist for ensemble (re-run caretlist to filter correlated models etc.)
# ensemble <- caretEnsemble::caretStack(models,trControl = fitControl,metric="RMSE")
# # Predict
# custompredict(ensemble,b2)
#----Plot TVC vs Time grouped by Temperature----
plot <- plotTVC_time_temp(data = rawdata, title = "rawdata")
print(plot)
#----PCA----
plot <- doPca(data = rawdata, ncomp = c(1,2))
print(plot)
#----(NOT USED IN PAPER)----
# # Sample/variable selection by PCA
# res.pca <- PCA(data[,-1], graph=F) # PCA
# eig.val <- get_eigenvalue(res.pca) # get eigens
# eig.val # variance.percent > 1 means its useful
# fviz_pca_var(res.pca, col.var = "black") # Plot arrows and stuff
# 
# var <- get_pca_var(res.pca)
# var
# corrplot(var$cos2,is.corr=F)
# fviz_cos2(pca1, choice = "var", axes = 1:5) # Contribtuion of variables
# fviz_pca_var(res.pca, col.var = "cos2",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE # Avoid text overlapping
# )
# 
# contrib <- fviz_contrib(pca1, choice = "ind", axes = 1:2) # See contributions to pca 1 and 2 by samples
# 
# # Get Clean data
# data <- read.table("VM R1 AND R3.csv",sep = ",", header = TRUE)
# data <- data[,1:20]
# rownames(data) <- data[,1]
# data <- data[,-1]
# 
# # Find contributors above 0.2
# res.pca <- PCA(data[,-1], graph=F) # PCA
# contrib <- fviz_contrib(res.pca, choice = "ind", axes = 1:2) # See contributions to pca 1 and 2 by samples
# contribList <- as.data.frame(contrib[["plot_env"]][["data"]][["contrib"]])
# contribList <- cbind(as.character(contrib[["plot_env"]][["data"]][["name"]]), contribList)
# contribList[,1] <- as.character(contribList[,1])
# contribList <- contribList[order(contribList[,2]), ]
# rownames(contribList) <- contribList[,1]
# arcik <- rownames(contribList)
# contribList <- as.data.frame(contribList[,-1])
# rownames(contribList) <- arcik
# colnames(contribList) <- "Contribution"
# contribTrain <- contribList[contribList[,2] > 0.15, 1] # Threshold
# contribTest <- contribList[contribList[,2] < 0.15, 1]

# train <- subset(data,rownames(data)%in%contribList[54:nrow(contribList), 1]) # Subset data
# test <- subset(data,rownames(data)%in%contribList[1:53, 1] ) # Subset data
#----Plot heatmap----
plot <- finalheatmap(resultlist = resultlist, names = c("Batch 1","Batch 2","Combined"), 
                     limits = c(65,90))
print(plot)
#----Plot comparison of two samples' intensities over wavelengths----
plot <- msi_compare_two_samples(data, "96h_0A", "72h_15B")
print(plot)
#----Save plot with anti-aliasing (prettier)----
ggsave(plot = plot, "plot.png",type = "cairo-png")
#----Dunno---- Might be mean(avg) intensity for batches
plot <- ggplot(data,aes(x=wave,y=data[,1])) +
  geom_line() +
  geom_line(aes(data=wave,y=y[,1],col="red")) +
  theme_bw(base_family = "DMOFCB+AdvGulliv-R") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5), legend.position = "none") +
  xlab("Mean") +
  ylab("Intensity") +
  scale_x_continuous(breaks=c(1:18))