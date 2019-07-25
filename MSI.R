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
library(rms)
library(factoextra)
library(FactoMineR)
library(tidyverse)
library("viridis")
source("Util.R")
library(prospectr)
library(hyperSpec)
library(car)
options(scipen=999)
#----Read----
rawdata <- read.table("VM R1 AND R3.csv",sep = ",", header = TRUE)

rawdata <- rawdata[1:110,] # Batch 1

rawdata <- rawdata[111:nrow(rawdata),] # Batch 2

rownames(rawdata) <- rawdata[,1]
rawdata <- rawdata[,-1]

data <- rawdata
# Outlier test
lm <- lm(TVC~., data=data)
cooksd <- cooks.distance(lm)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 5*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>5*mean(cooksd, na.rm=T),names(cooksd),""), col="red")
#----Split to train and test----
# 0 = naive, 1 = time, 2 = temp, 3 = TVC, 4=Ken-stone
return <- split(data, 0)

train <- return$train
test <- return$test

# Get good train
# train <- data[rownames(data) %in% rownames(RMSEbestpls$trainingData),]
# test <- data[!rownames(data) %in% rownames(RMSEbestpls$trainingData),]
#----Train models----
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5
                           ,savePredictions = "final"
                           #,index = createResample(data$TVC, 25)
                           ,allowParallel = T
                           )
gridlist <- gridListVM(data)
cl <- makePSOCKcluster(3)
registerDoParallel(cl) # Start parellel processing
models <-  caretEnsemble::caretList(TVC~., data=data, trControl = fitControl, metric = "RMSE"
                                    ,preProcess = c("center","scale")
                                    ,tuneList = list(
                                      pls=caretModelSpec(method="pls", tuneGrid=gridlist$plsgrid)
                                      ,svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                                      ,rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid)
                                      ,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                                      ,pcr=caretModelSpec(method="pcr", tuneGrid=gridlist$plsgrid)
                                      ,ridge=caretModelSpec(method="ridge")
                                      ,lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid)
                                      ,nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, linout = TRUE, maxit=1000)
                                    )
)
stopCluster(cl) # Stop parellel processing
registerDoSEQ()

custompredict(models$pls, data) # Predict a single model and plot

plotvariableimportance(models$pls) # Plot variable importance for a single model

results <- resamples(models) # RMSE results over 25 resamples from CV phase
summary(results, metric = "RMSE")
print(bwplot(results,metric = "RMSE", main="RMSE over 25 resamples"))

testresults <- multiplemodelpredict(models) # Model results on test set plotted
testresults
#----Multiple iterations----
iters <- 100
results <- multipleiterations(data, iters, "VM")

rmseResults <- results$rmseResults # RMSE results for each iteration for each model
predAcc <- results$predAcc # Accuracy results for each iteration for each model
RMSEbestlist <- results$RMSEbestlist # Best models by RMSE over all iterations
ACCbestlist <- results$ACCbestlist # Best models by accuracy over all iterations
trainsamplecount <- results$trainsamplecount # Training samples appereance count on best models

# RMSE results
rmseResults$min <- apply(rmseResults,1,min) # Min as extra col
rmseResults$mean <- apply(rmseResults,1,FUN=mean) # Median as extra col
rmseResults$max <- apply(rmseResults,1,max) # Max as extra col
rmseResults <- rmseResults[order(rmseResults$mean),] # Order by median

# Box plot for all rmse iters
boxplot(t(rmseResults[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(8), 
        main=paste("RMSE over", iters,"iterations"), cex.axis=0.7)

# Accuracy results
predAcc$min <- apply(predAcc,1,min) # Same as above
predAcc$mean <- apply(predAcc,1,FUN=mean)
predAcc$max <- apply(predAcc,1,max)
predAcc <- predAcc[order(-predAcc$mean),]

# Box plot for all accuracy iters
boxplot(t(predAcc[,-which(colnames(rmseResults) %in% c("min","median","max"))]), col=viridis(8), 
        main=paste("Accuracy over", iters,"iterations"), cex.axis=0.7)

# Combined final results for min/max
finalresults <- data.frame(matrix(nrow=8,ncol=2)) # Create results DF
rownames(finalresults) <- rownames(predAcc)
colnames(finalresults) <- c("RMSE", "Accuracy")
finalresults[,1] <-  round(rmseResults$min,2) # Minimum RMSE for models
finalresults[,2] <- round(predAcc$max,2) # Max accuracy for models
finalresults$name <- rownames(finalresults)

# Mean combined final results
finalresults <- data.frame(matrix(nrow=8,ncol=2)) # Create results DF
rownames(finalresults) <- rownames(predAcc)
colnames(finalresults) <- c("RMSE", "Accuracy")
finalresults[,1] <-  round(rmseResults$mean,2) # Minimum RMSE for models
finalresults[,2] <- round(predAcc$mean,2) # Max accuracy for models
finalresults$name <- rownames(finalresults)

# Measuring how many times a sample has appeared in best models' training sets
# Extract time/temp and AorB as new columns
countbest <- transform(countbest, Temperature=stringr::str_extract(name,"[0-9]*((?=C)|(?=B)|(?=A))"))
countbest <- transform(countbest, AorB=stringr::str_extract(name,"[A-B]"))
countbest <- transform(countbest, Time=as.numeric(stringr::str_extract(name,"[0-9]*(?=h)")))

# Count how many times a storage time has appeared
timecount <- aggregate(countbest$Count, by=list(Time=countbest$Time), FUN=sum)
colnames(timecount)[2] <- "Count"

# Count how many times a storage temperature has appeared
tempcount <- aggregate(countbest$Count, by=list(Temperature=countbest$Temperature), FUN=sum)
colnames(tempcount)[2] <- "Count"

# Count how many times A or B has appeared
abcount <- aggregate(countbest$Count, by=list(AorB=countbest$AorB), FUN=sum)
colnames(abcount)[2] <- "Count"

#----Ensemble----
# modelCor <- as.data.frame(modelCor(resamples(models)))
# ensemble <- caretEnsemble::caretEnsemble(models,trControl = fitControl, metric="RMSE")
#----Plot TVC vs Time grouped by Temperature----
# Extract temperature, AorB, Time from sample names
data$Sample <- rownames(data)
data <- transform(data, Temperature=stringr::str_extract(Sample,"[0-9]*((?=C)|(?=B)|(?=A))"))
data <- transform(data, Time=as.numeric(stringr::str_extract(Sample,"[0-9]*(?=h)")))

data <- as.data.frame(apply(data, 2, function(x) gsub("^$|^ $", 0, x))) # Change empty temperatures to 0 

timeandtemp <- data %>% dplyr::select("TVC", "Time", "Temperature")
timeandtemp$Time <- as.numeric(as.character(timeandtemp$Time))
timeandtemp$TVC <- as.numeric(as.character(timeandtemp$TVC))
timeandtemp$Temperature <- factor(timeandtemp$Temperature, levels=c(15,10,5,0), ordered=T)
timeandtemp <- timeandtemp[order(timeandtemp$Time),]

ggplot(data=timeandtemp, aes(x=Time, y=TVC, group=Temperature)) +
  geom_line(aes(color=Temperature)) +
  scale_color_manual(values=c("red", "black", "blue","green")) +
  theme_classic()
#----PCA----
# Create new column for TVC intervals
pcadata <- data
quantiles <- quantile(pcadata$TVC, c(.33, .66)) # Find quantiles for TVC values
pcadata$TVCinterval <- ifelse(pcadata$TVC <= quantiles[1], 'low',
                           ifelse(pcadata$TVC > quantiles[1] & pcadata$TVC <= quantiles[2], 'med',
                                  ifelse(pcadata$TVC > quantiles[2], 'high', '')))

pca <- prcomp(pcadata[, c(-1,-20)])
print(ggbiplot::ggbiplot(pca, choices = c(1,2),
                         groups = pcadata$TVCinterval, ellipse = TRUE, var.axes = FALSE))

pca <- PCA(pcadata, graph=F)
plot(pca)

#----Sample/variable selection by PCA----
# res.pca <- PCA(data[,-1], graph=F) # PCA
# eig.val <- get_eigenvalue(res.pca) # get eigens
# eig.val # variance.percent > 1 means its useful
# fviz_pca_var(res.pca, col.var = "black") # Plot arrows and stuff
# 
# var <- get_pca_var(res.pca)
# var
# corrplot(var$cos2,is.corr=F)
# fviz_cos2(res.pca, choice = "var", axes = 1:2) # Contribtuion of variables
# fviz_pca_var(res.pca, col.var = "cos2",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE # Avoid text overlapping
# )
# 
# contrib <- fviz_contrib(res.pca, choice = "ind", axes = 1:4) # See contributions to pca 1 and 2 by samples
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

