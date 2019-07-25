# option 0=naive, 1=time, 2=temperature, 3=TVC
split <- function(data, option) 
{
  set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)) # Set seed
  
  if(option == 0)
  {
    train_index = caret::createDataPartition(y=data$TVC, p=0.7, list=FALSE, times=1)
    
    train <- data[train_index,]
    test <- data[-train_index,]
  }
  else if(option == 1)
  {
    data <- data[order(data$Time), ] # Order by storage times
    
    low <- data[70, 'Time' ]
    med <- data[140, 'Time']
    
    # Create new column for time intervals
    data$Timeinterval <- ifelse(data$Time <= low, 'low',
                                ifelse(data$Time > low & data$Time <= med, 'med',
                                       ifelse(data$Time > med, 'high', '')))
    
    # Find smallest dataset and assign k from that (taking equal from all intervals)
    k <- min(length(which(data$Timeinterval == 'low')),
             length(which(data$Timeinterval == 'med')),
             length(which(data$Timeinterval == 'high'))) * 0.7
    
    train <- data.frame()
    test <- data.frame()
    for(level in c('low','med','high'))
    {
      dataLevel <- data[data$Timeinterval == level,] # Only select current interval
      
      p <-  k / nrow(dataLevel) # Assign p for caret createDataPartition function
      
      train_index <- caret::createDataPartition(y=dataLevel$TVC, p=p, list=FALSE, times=1)
      
      # Append the new sample rows
      train <- rbind(train, dataLevel[train_index,])
      test <- rbind(test, dataLevel[-train_index,])
    }
  }
  else if(option == 2)
  {
    # Find smallest dataset and assign k from that (taking equal from all intervals)
    k <- min(length(which(data$Temperature == 0)),
             length(which(data$Temperature == 5)),
             length(which(data$Temperature == 10)),
             length(which(data$Temperature == 15))) * 0.7 + 6
    
    train <- data.frame()
    test <- data.frame()
    for(temp in c(0,5,10,15))
    {
      dataTemperature = data[data$Temperature == temp,] # Only select current temp
      
      p =  k / nrow(dataTemperature) # Assign p for caret createDataPartition function
      
      train_index =  caret::createDataPartition(y=dataTemperature$TVC, p=p, list=FALSE, times=1)
      
      # Append the new sample rows
      train <- rbind(train, dataTemperature[train_index,])
      test <- rbind(test, dataTemperature[-train_index,])
    }
  }
  else if(option == 3)
  {
    quantiles <- quantile(data$TVC, c(.33, .66)) # Find quantiles for TVC values
    
    # Create new column for TVC intervals
    data$TVCinterval <- ifelse(data$TVC <= quantiles[1], 'low',
                               ifelse(data$TVC > quantiles[1] & data$TVC <= quantiles[2], 'med',
                                      ifelse(data$TVC > quantiles[2], 'high', '')))
    
    # Find smallest dataset and assign k from that (taking equal from all intervals)
    k <- min(length(which(data$TVCinterval == 'low')),
             length(which(data$TVCinterval == 'med')),
             length(which(data$TVCinterval == 'high'))) * 0.7
    
    train <- data.frame()
    test <- data.frame()
    for(level in c('low','med','high'))
    {
      dataLevel <- data[data$TVCinterval == level,] # Only select current interval
      
      p <-  k / nrow(dataLevel) # Assign p for caret createDataPartition function
      
      train_index <- caret::createDataPartition(y=dataLevel$TVC, p=p, list=FALSE, times=1)
      
      # Append the new sample rows
      train <- rbind(train, dataLevel[train_index,])
      test <- rbind(test, dataLevel[-train_index,])
    }
    
    # Remove interval column
    train <- train[,-ncol(train)]
    test <- test[,-ncol(test)]
  }
  else if(option == 4)
  {
    ken <- kenStone(data, k=nrow(data)*0.7, metric="euclid")
    train <- data[ken$model,]
    test <- data[ken$test,]
  }
  
  return(list("train" = train, "test" = test))
}

custompredict <- function(model, set, split=0)
{
  # For cross-validation set
  if(is.null(set))
  {
    prediction <- model$finalModel$predicted
    
    rmse <- getTrainPerf(model)[,1]
    
    difference <- as.data.frame(abs(model$trainingData$.outcome-prediction))
    accuracy <- sum(difference < 1) / nrow(difference)
    
    x <- cbind(prediction,model$trainingData$.outcome)
    x <- as.data.frame(x[order(x[,2]),])
    
    setName <- "CV"
  }
  else
  {
    # Predict
    prediction <- stats::predict(model, set)
    
    # RMSE calc
    rmse <- Metrics::rmse(set$TVC, prediction)
    
    # Within 1 cfu accuracy calc
    difference <- as.data.frame(abs(set$TVC-prediction))
    accuracy <- sum(difference <= 1) / nrow(difference)
    
    x <- cbind(prediction,set$TVC)
    x <- as.data.frame(x[order(x[,2]),])

    setName <- deparse(substitute(set))
  }
  
  if(split == 0) {
    split <- "Naive split"
  }else if(split == 1) {
    split <- "Time split"
  }else if(split == 2) {
    split <- "Temp split"
  }else if(split == 3) {
    split <- "TVC split"
  }
  
  title <- paste(split,",",model$method,",",setName)
  colnames(x) <- c("Predicted", "Actual")
  plot <- ggplot(data=x, aes(x=Actual,y=Predicted)) +
    geom_point() +
    ggtitle(title) +
    ylim(c(0,9)) +
    xlim(c(0,9)) +
    geom_abline(col="blue") +
    geom_abline(intercept = 1, col="red") +
    geom_abline(intercept = -1, col="red") +
    geom_label(label=paste("Accuracy: ", round(accuracy, digits = 2)), x=7, y=3, size=6) +
    geom_label(label=paste("RMSE: ", round(rmse, digits = 2)), x=7, y=3.5, size=6)
  
  print(plot)
}

customSVM <- function(train)
{
  tune <- tune.svm(train[,-1], train[,1],cost = seq(1, 50, 10),
                   gamma = seq(0, 1, 0.1),
                   epsilon = seq(0, 1, 0.1))
  
  model <- e1071::svm(TVC~., data=train, cost=tune$best.parameters$cost,
                      gamma = tune$best.parameters$gamma,
                      epsilon = tune$best.parameters$epsilon)
  
  return(model)
}

customRF <- function(train)
{
  tune <- tune.randomForest(train[, -1], train[, 1], mtry = seq(2,18,1)
                            , nodesize = seq(2,10,1))
  
  model <- randomForest(TVC~., data=train, mtry = tune$best.parameters$mtry, 
                        ntree = 1000, 
                        nodesize = tune$best.parameters$nodesize, importance=T) 
  
  return(model)
}

keepwave <- function(data, list)
{
  # Transpose to filter wave-lengths
  data_T <- as.data.frame(t(data))
  # Add waves as a column
  data_T$wave <- as.numeric(stringr::str_extract(rownames(data_T[,-1]),"[0-9]*(?=\\.)"))
  # Add temporary 0 to TVC wave
  data_T["TVC","wave"] <- 0
  
  list <- c(list, 0)
  # Filter wavelengths (include TVC)
  x <- data_T[data_T$wave %in% list, ]
  
  # Transpose back
  x <- as.data.frame(t(x))
  
  return(x[-nrow(x),])
}

gridListVM <- function(data)
{
  plsgrid <- expand.grid(.ncomp = c(1:12)) # PLS
  
  svmRadialgrid <- expand.grid(sigma = 2^c(-25, -20, -15, -10, -5, -0), C = 2^c(0:5)) # svmRadial
  
  svmLineargrid <- expand.grid(C = 2^c(0:5)) # svmLinear
  
  knngrid <- expand.grid(k = seq(3,7,2)) # KNN
  
  nnetgrid <- expand.grid(size = c(1,3,5),
                          decay = c(0.1,0.5,1)) # Nnet
  
  brnngrid <- expand.grid(neurons=c(1:1)) # Brnn
  
  larsgrid <- expand.grid(fraction = seq(from = 0, to = 1, length = 100)) # Lars
  
  gbmgrid <- expand.grid(n.trees = c(1000,1500),
                         interaction.depth=c(1:3), shrinkage=c(0.01,0.05,0.1),
                         n.minobsinnode=c(20)) # gbm
  
  tune <- tuneRF(data[,-1], data[,1], trace = FALSE, plot = FALSE)
  mtry <- as.numeric(row.names(tune)[(which(tune[,2]==min(tune[,2])))])
  rfgrid <- expand.grid(.mtry = mtry)
  
  glmnetgrid <- expand.grid(lambda = seq(0.001,0.1,by = 0.001), alpha=1) #glmnet
  
  xgbTreegrid <- expand.grid(nrounds = 1, eta = 0.3, max_depth = 5,gamma = 0,
                             colsample_bytree=1, min_child_weight=1) # xgbTree
  
  xgbLineargrid <- expand.grid(nrounds= 2400,lambda = 1,alpha=0) # xgbLinear
  
  # Create a list of grids
  gridlist <- list(plsgrid=plsgrid, svmRadialgrid=svmRadialgrid, svmLineargrid=svmLineargrid,
                   knngrid=knngrid,nnetgrid=nnetgrid, larsgrid=larsgrid, gbmgrid=gbmgrid,
                   rfgrid=rfgrid, glmnetgrid=glmnetgrid, nnetgrid=nnetgrid, 
                   xgbTreegrid=xgbTreegrid, xgbLineargrid=xgbLineargrid, brnngrid=brnngrid)
  
  return(gridlist)
}

gridListFTIR <- function(data)
{
  plsgrid <- expand.grid(.ncomp = seq(from=1,to=15,by=1))
  
  svmRadialgrid <- expand.grid(sigma = 2^c(-25, -20, -15, -10, -5, -0), C = 2^c(0:5))
  
  svmLineargrid <- expand.grid(C = 2^c(0:5))
  
  knngrid <- expand.grid(k = 1:40)
  
  nnetgrid <- expand.grid(size = 1, decay = c(0.5,1))
  
  larsgrid <- expand.grid(fraction = seq(from = 0, to = 1, by = 0.1))
  
  tune <- tuneRF(data[,-1], data[,1], trace = FALSE, plot = FALSE)
  mtry <- as.numeric(row.names(tune)[(which(tune[,2]==min(tune[,2])))])
  rfgrid <- expand.grid(.mtry = mtry) 
  
  # Create a list of grids
  gridlist <- list(plsgrid=plsgrid, svmRadialgrid=svmRadialgrid, svmLineargrid=svmLineargrid,
                   knngrid=knngrid, nnetgrid=nnetgrid, larsgrid=larsgrid,
                   rfgrid=rfgrid, nnetgrid=nnetgrid)
  
  return(gridlist)
}

findandkeepwavesRF <- function(data, topx=0, threshold=0) # no topx argument returns the list
{
  cl <- makePSOCKcluster(3)
  registerDoParallel(cl)
  tune <- tuneRF(data[,-1], data[,1], trace = FALSE, plot = FALSE)
  mtry <- as.numeric(row.names(tune)[(which(tune[,2]==min(tune[,2])))])
  
  model <- randomForest(TVC~., data=data, mtry = mtry, ntree=1000, importance=T)
  importance <- as.data.frame(model$importanceSD)
  importance <- transform(importance, wave=stringr::str_extract(rownames(importance),"([0-9]*)(?=\\.)"))
  importance$wave <- as.numeric(as.character(importance$wave))
  
  stopCluster(cl)
  registerDoSEQ()
  
  if(topx != 0) {
    keeplist <- top_n(importance, topx, importance$model.importanceSD)$wave
  } else if(threshold !=0){
    keeplist <- importance[importance$model.importanceSD > threshold,]$wave
  } else {
    return(importance)
  }
  
  data <- keepwave(data, keeplist)
  return(data)
}

compareIntensity <- function(data, list)
{
  ar <- data[rownames(data) %in% list,]
  ar <- as.data.frame(t(ar))
  ar <- transform(ar, wave=stringr::str_extract(rownames(ar),"[0-9]+\\.+[0-9]*"))
  ar$wave <- as.numeric(as.character(ar$wave)) # Set as numeric
  difference <- (ar[,1] - ar[,2])[1]
  ar <- ar[-nrow(ar), ]
  
  if(difference < 0) {
    fresh <- ar[, 1]
    spoiled <- ar[ ,2]
  } else {
    fresh <- ar[, 2]
    spoiled <- ar[ ,1]
  }
  
  ggplot(data=ar, aes(x=ar[,3])) +
    geom_line(aes(y=fresh, color="Fresh")) +
    geom_line(aes(y=spoiled, color="Spoiled")) +
    scale_color_manual(values=c("green","red")) +
    labs(color="Level") +
    theme_classic()
  
  # plot(ar[-1, 3], abs(ar[-1,4] / ar[1,4]), type="l", main="Intensity difference / TVC difference")
}

meanintensitygraph <- function(data, low, high)
{
  dataF <- data[data$TVC <= low,]
  dataS <- data[data$TVC >= high, ]
  
  dataF <- as.data.frame(t(dataF))
  dataF <- dataF[-1,]
  dataFsums <- data.frame(matrix(ncol=2,nrow=ncol(data)-1))
  dataFsums[,1] <- rowSums(dataF) / ncol(dataF)
  rownames(dataFsums) <- rownames(dataF)
  dataFsums <- transform(dataFsums, wave=stringr::str_extract(rownames(dataF),"([0-9]*)(?=\\.)"))
  dataFsums <- dataFsums[,-2]
  dataFsums$wave <- as.numeric(as.character(dataFsums$wave))
  
  dataS <- as.data.frame(t(dataS))
  dataS <- dataS[-1,]
  dataSsums <- data.frame(matrix(ncol=2,nrow=ncol(data)-1))
  dataSsums[,1] <- rowSums(dataS) / ncol(dataS)
  rownames(dataSsums) <- rownames(dataS)
  dataSsums <- transform(dataSsums, wave=stringr::str_extract(rownames(dataS),"([0-9]*)(?=\\.)"))
  dataSsums <- dataSsums[,-2]
  dataSsums$wave <- as.numeric(as.character(dataSsums$wave))
  
  dataSub <- data.frame(matrix(ncol=4,nrow=ncol(data)-1))
  dataSub[,1] <- dataFsums$X1
  dataSub[,2] <- dataSsums$X1
  dataSub[,3] <- dataFsums$wave
  dataSub[,4] <- abs(dataSub[,2] - dataSub[,1])
  dataSub <- dataSub[-nrow(dataSub), ]
  colnames(dataSub) <- c("FreshSums", "SpoiledSums", "Wave","Difference")
  
  plot(dataSub[,3], dataSub[,1], type="l", xlab="Wave", ylab="Intensity", col="green")
  legend(x=1000, y=1.5, c("Fresh","Spoiled"),cex=1,col=c("green","red"),pch=c(15,15), bty = "n")
  lines(dataSub[,3], dataSub[,2], col="red")
  
  plot(dataSub[,3], abs(dataSub[,2] - dataSub[,1]), type="l", xlab="wave",ylab = "Abs intensity difference")
  
  return(dataSub)
}

plotSTD <- function(data)
{
  data_T <- as.data.frame(t(data)) # Transpose
  data_T <- data_T[-which(rownames(data_T) %in% "TVC"),] # Remove TVC row
  data_T <- transform(data_T,SD=apply(data_T,1,sd)) # SD of every row
  data_T$waves <- as.numeric(stringr::str_extract(rownames(data_T),"[0-9]*(?=\\.)")) # Create waves as column
  
  plot <- ggplot(data=data_T, aes(x=waves,y=SD)) +
    geom_line()+
    xlab("Wave length")+
    ylab("Std")+
    ggtitle("Standard deviation of intensity over wave lengths") +
    scale_x_continuous("Wave", breaks=seq(500,4000,500)) +
    theme_classic()
  
  print(plot)
  
  return(data_T[,c("SD","waves")])
}

filterwavesfromSTD <- function(data, std, threshold)
{
  std <- std[which(std$SD > threshold),]
  names <- c(rownames(std),'TVC')
  data <- data[, colnames(data) %in% names]
  
  return(data)
}

prep <- function(data, methodlist, pcaComp=0)
{
  data <- data[, -1]
  
  if(pcaComp > 0) {
    pp <- preProcess(data, method = methodlist, pcaComp = pcaComp)
  } else {
    pp <- preProcess(data, method = methodlist)
  }
  
  data <- predict(pp, data)
  
  return(data)
}

createspc <- function(data)
{
  data_T <- as.data.frame(t(data))
  data_T <- data_T[-1,] # remove TVC row
  data_T$waves <- as.numeric(stringr::str_extract(rownames(data_T),"[0-9]+\\.+[0-9]*"))
  
  data <- as.data.frame(t(data_T))
  spc <- new("hyperSpec", spc=data[-nrow(data), ], wavelength=as.numeric(data["waves",]))
  
  return(spc)
}

applypreprocess <- function(data, preproc, sgwindow=0, deriv=1)
{
  spc <- createspc(data)
  if(preproc == "baseline") {
    bl <- spc.fit.poly.below(spc)
    spc <- spc - bl
    data <- as.data.frame(spc$spc)
  } else if(preproc == "snv") {
    data <- as.data.frame(standardNormalVariate(data)) 
  } else if(preproc == "snv-d") {
    data <- as.data.frame(detrend(X=spc$spc, wav=attributes(spc)$wavelength))
  } else if(preproc == "s-g"){
    data <- as.data.frame(savitzkyGolay(data, m=2,p=3,w=sgwindow))
  } else if(preproc == "derivation"){
    data <- as.data.frame(t(diff(t(spc$spc), differences =deriv)))
  } else if(preproc == "normal"){
    spc <- sweep(spc,1,mean,"/")
    data <- as.data.frame(spc$spc)
  }

  return(data)
}

plotvariableimportance <- function(model)
{
  imp <- as.data.frame(varImp(model)$importance)
  imp <- transform(imp, wave=as.numeric(stringr::str_extract(rownames(imp),"[0-9]+\\.+[0-9]*")))
  imp <- imp[order(imp$wave), ]
  ggplot(data=imp, aes(x=wave, y=Overall)) +
    geom_line()
}

multiplemodelpredict <- function(models, plot=F)
{
  rmseResults <- data.frame(matrix(nrow=length(models),ncol=3))
  colnames(rmseResults) <- c("x","RMSE", "Accuracy")
  for(i in 1:length(models))
  {
    model <- models[[i]]
    
    # RMSE
    prediction <- predict(model, test)
    rmse <- Metrics::rmse(test$TVC, prediction)
    
    # Within 1 cfu accuracy
    difference <- as.data.frame(abs(test$TVC-prediction))
    accuracy <- sum(difference <= 1) / nrow(difference)
    
    if(plot)
    {
      custompredict(model, test, 0)
    }
    rmseResults[i,] <- list(model$method, rmse, accuracy)
  }
  rownames(rmseResults) <- rmseResults[,1]
  rmseResults <- rmseResults[, -1, drop=F]
  
  return(round(rmseResults,2))
}

lmcoefficient <- function(data)
{
  model <- lm(TVC~., data=data)
  
  coef <- as.data.frame(model$coefficients)
  coef <- transform(coef, wave=as.numeric(stringr::str_extract(rownames(coef),"([0-9]*)(?=\\.)")))
  coef <- coef[-1,]
  coef <- na.omit(coef)
  
  plot <- ggplot(data=coef, aes(x=wave, y=model.coefficients)) +
    geom_line()
  
  print(plot)
  
  return(coef)
}

filterwavesbyPCAcontrib <- function(data, contribution, threshold)
{
  filterlist <- rownames(contribution[contribution$sum > threshold,]) # Filter contribs below threshold
  filterlist <- str_extract(filterlist, "[0-9]*(?=\\.)") # Convert row names to wave numbers
  data <- keepwave(data, filterlist) # Filter the data 
  
  return(data)
}

findpatternbestmodels <- function(RMSEbestlist)
{
  RMSEbesttrain <- cbind(rownames(RMSEbestlist[[1]]$trainingData), rownames(RMSEbestlist[[2]]$trainingData))
  for(i in 3:length(RMSEbestlist))
  {
    RMSEbesttrain <- cbind(RMSEbesttrain, rownames(RMSEbestlist[[i]]$trainingData))
  }
  
  countbest <- data.frame(matrix(ncol=2)) # create DF to store count for each sample appeared in training sets
  for(rowname in rownames(data)) # Check each sample 
  {
    countbest <- rbind(countbest, length(which(RMSEbesttrain == rowname)))
  }
  
  countbest <- countbest[-1,] # remove NA row
  rownames(countbest) <- rownames(data) # set row names as sample names
  countbest <- countbest[,-2, drop=F] # drop unnecessary
  colnames(countbest)[1] <- "Count"
  countbest$pct <- round(countbest$Count / length(RMSEbestlist), 2) # What percent did the sample got used in all iters
  
  return(countbest)
}

multipleiterations <- function(data, iters, datatype = "VM")
{
  x <- iters
  methodlist <- c("pls","svmLinear","rf", "knn", "pcr","ridge","lars","nnet")
  
  rmseResults <- data.frame(matrix(nrow=length(methodlist), ncol=1)) # DF to store RMSE results
  rownames(rmseResults) <- methodlist
  
  predAcc <- data.frame(matrix(NA, nrow=length(methodlist), ncol=1)) # DF to save acc results
  rownames(predAcc) <- methodlist
  
  max <- 0
  min <- 50
  cl <- makePSOCKcluster(3) # Parelell processing with 3 cores
  registerDoParallel(cl)
  
  for(i in 1:x)
  {
    print(paste("Starting iteration", i ,"..."))
    starttime <- Sys.time()
    set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
    
    # Split to train and test
    return <- split(data, 0) 
    train <- return$train
    test <- return$test
    
    if(datatype == "VM") {
      gridlist <- gridListVM(data)
    } else if(datatype == "FTIR") {
      gridlist <- gridListFTIR(data)
    }
    
    fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                               savePredictions = "final"
                               ,index = createResample(train$TVC, 25)
                               ,allowParallel = T)
    
    models <-  caretEnsemble::caretList(TVC~., data=train, trControl = fitControl, metric = "RMSE", continue_on_fail = T
                                        ,preProcess = c("center","scale")
                                        ,tuneList = list(
                                          pls=caretModelSpec(method="pls", tuneGrid=gridlist$plsgrid)
                                          ,svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                                          ,rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid, ntree=1000)
                                          ,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                                          ,pcr=caretModelSpec(method="pcr", tuneGrid=gridlist$plsgrid)
                                          ,ridge=caretModelSpec(method="ridge")
                                          ,lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid, use.Gram=F)
                                          ,nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, linout = TRUE, maxit=1000)
                                        )
    )
    
    # Predictions for test set
    predictionsDF <- as.data.frame(predict(models, newdata = test))
    
    # Calculate prediction difference
    predDiff <- abs(predictionsDF - test$TVC)
    
    # Calculate rmse for each model
    rmseVector <- c()
    for(k in 1:length(models))
    {
      # RMSE
      rmse <- Metrics::rmse(test$TVC,predictionsDF[,k]) # Calculate rmse for the current model and iter
      rmseVector <- c(rmseVector, rmse) # Add to the rmselist for all rmse on this iter
      
      model <- models[[k]] # specify current model
      
      # If it's the first iteration, should save the model
      if(i == 1) {
        min <- 50
      }
      else {
        temp <- as.data.frame(rmseResults[,-1]) # find minimum rmse of this model so far
        min <- min(temp[k,])
      }
      
      # Capture best model by RMSE
      if(rmse < min) # If the current model has lower RMSE than all the previous ones, save it
      {
        if(k==1) {
          RMSEbestpls <- model
        } else if(k==2) {
          RMSEbestsvmLinear <- model
        } else if(k==3) {
          RMSEbestrf <- model
        } else if(k==4) {
          RMSEbestknn <- model
        } else if(k==5) {
          RMSEbestpcr <- model
        } else if(k==6) {
          RMSEbestridge <- model
        } else if(k==7) {
          RMSEbestlars <- model
        } else if(k==8) {
          RMSEbestnnet <- model
        } 
      }
      
      # ACCURACY
      # Check how many predictions are within 1 unit
      table <- table(predDiff[,k] < 1)
      pred <- as.numeric(table[2]) / as.numeric(table[1] + table[2])
      predAcc[k,i] <- pred # Store current model+iter accuracy on the DF
      
      max <- max(predAcc[k, ]) # Find max accuracy so far on this model type
      
      # Capture best model by Accuracy
      if(i == 1 | pred >= max) # If first iter or the accuracy is better than any before it, save it
      {
        if(k==1) {
          ACCbestpls <- model
        } else if(k==2) {
          ACCbestsvmLinear <- model
        } else if(k==3) {
          ACCbestrf <- model
        } else if(k==4) {
          ACCbestknn <- model
        } else if(k==5) {
          ACCbestpcr <- model
        } else if(k==6) {
          ACCbestridge <- model
        } else if(k==7) {
          ACCbestlars <- model
        } else if(k==8) {
          ACCbestnnet <- model
        } 
      }
    }
    
    # Add RMSEs together for this iteration
    rmseResults <- cbind(rmseResults,rmseVector)
    
    print(paste("This iteration took:", round(Sys.time() - starttime,3), "minutes"))
  }
  stopCluster(cl) # Stop parellel processing
  registerDoSEQ()
  
  rmseResults <- rmseResults[,-1, drop=F]
  colnames(predAcc) <- 1:x 
  colnames(rmseResults) <- 1:x
  
  rmseResults <- round(rmseResults,2)
  predAcc <- round(predAcc,2)
  
  RMSEbestlist <- list(RMSEbestpls=RMSEbestpls,RMSEbestsvmLinear=RMSEbestsvmLinear,
                       RMSEbestrf=RMSEbestrf, RMSEbestknn=RMSEbestknn,
                       RMSEbestpcr=RMSEbestpcr, RMSEbestridge=RMSEbestridge,
                       RMSEbestlars=RMSEbestlars, RMSEbestnnet=RMSEbestnnet)
  
  ACCbestlist <- list(ACCbestpls=ACCbestpls,ACCbestsvmLinear=ACCbestsvmLinear,
                      ACCbestrf=ACCbestrf, ACCbestknn=ACCbestknn,
                      ACCbestpcr=ACCbestpcr, ACCbestridge=ACCbestridge,
                      ACCbestlars=ACCbestlars, ACCbestnnet=ACCbestnnet)
  
  trainsamplecount <- findpatternbestmodels(RMSEbestlist)
  
  returnlist <- list(RMSEbestlist=RMSEbestlist, ACCbestlist=ACCbestlist, rmseResults=rmseResults,
                     predAcc=predAcc, trainsamplecount=trainsamplecount)
  
  return(returnlist)
}