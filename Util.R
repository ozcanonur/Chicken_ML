# Reads the MSI file
# outliers = read file with outliers removed or not, batch = batch 1, 2 or 0 (all data)
# centerandscale = preprocess centerscale 
readMSI <- function(outliersRemoved = F, batch = 0, centerandscale = T)
{
  if(outliersRemoved) {
    index <- 110
    data <- read.table("VM R1 AND R3_NEWNEW.csv",sep = ",", header = TRUE, row.names = 1)
  } else {
    index <- 115
    data <- read.table("VM R1 AND R3_Original.csv",sep = ",", header = TRUE, row.names = 1)
  }
  
  if(centerandscale) {
    data <- applypreprocess(data,c("center","scale"))
  }
  
  if(batch == 1){
    return(data[1:index,])
  } else if(batch == 2){
    return(data[(index+1):nrow(data),])
  } else {
    return(data)
  }
}

# Same as MSI, except you can keep specific waves with waves parameter
readFTIR <- function(outliersRemoved = F, batch = 0, waves = c(400:4000))
{
  if(outliersRemoved) {
    index <- 110
    data <- read.table("FTIR_newnew.csv",sep = ",", header = TRUE, row.names = 1)
  } else {
    index <- 116
    data <- read.table("FTIR_R1-R3.csv",sep = ",", header = TRUE, row.names = 1)
  }
  
  # Keep only specific waves
  data <- keepwave(data, waves)
  
  if(batch == 1){
    return(data[1:index,])
  } else if(batch == 2){
    return(data[(index+1):nrow(data),])
  } else {
    return(data)
  }
}

# Applies baseline correction, binning, normalisation and center and scale to a given data
apply_paper_preprocess <- function(datalist, bin = 4, normalarea = F, centerandscale = T)
{
  for(i in 1:length(datalist)){
    datalist[[i]] <- applyspectrapreprocess(datalist[[i]], "baseline") # Baseline correction
    datalist[[i]] <- applyspectrapreprocess(datalist[[i]], "bin", bin=bin) # Binning
    datalist[[i]] <- applyspectrapreprocess(datalist[[i]], "normal", normalarea = normalarea) # Normalisation on peak on 1500:1700
    
    if(centerandscale){
      datalist[[i]] <- applypreprocess(datalist[[i]],c("center","scale")) # Center and scale
    }
  }
  
  return(datalist)
}

# Performs the cooks distance test and plots it on a data frame
cooksdistancetest <- function(data)
{
  glm <- glm(TVC~., data=data)
  cooksd <- cooks.distance(glm)
  plot <- plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
  abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
  text(x=1:length(cooksd)+1, y=cooksd, col="red", labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""))
  
  
  return(plot)
}

# Plots time vs temperature grouped by temperature for a given data
plotTVC_time_temp <- function(data, title = "")
{
  data$Temperature <- as.numeric(as.character(stringr::str_extract(rownames(data),"[0-9]*((?=C)|(?=B)|(?=A))")))
  data$Time <- as.numeric(as.character(stringr::str_extract(rownames(data),"[0-9]*(?=h)")))
  
  # Extract tvc, time and temperature to another data frame
  data <- data[,c("TVC", "Time", "Temperature")]
  
  data <- as.data.frame(apply(data, 2, function(x) gsub("^$|^ $", 0, x))) # Change empty temperatures to 0 
  data[is.na(data)] <- 0 # Change NA temperatures to 0
  
  # Change TVC and time to numeric from factor
  data$TVC <- as.numeric(as.character(data$TVC))
  data$Time <- as.numeric(as.character(data$Time))
  
  data$Temperature <- factor(data$Temperature, levels=c(15,10,5,0), ordered=T) # Change order of temperatures
  data <- data[order(data$Time),] # Sort by time
  
  plot <- ggplot(data=data, aes(x=Time, y=TVC, group=Temperature)) +
    geom_line(aes(color=Temperature)) +
    ggtitle(title) +
    scale_color_manual(values=c("red", "black", "blue","green")) +
    theme_bw(base_family = "DMOFCB+AdvGulliv-R") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.text = element_text(size=15), plot.title = element_text(hjust=0.5)) +
    xlab("Storage time") +
    ylab(expression(paste("TVC log"[10]," cfu g"^-1)))
  
  return(plot)
}

# Does PCA on given data, prints summary and plots selected ncomps
doPca <- function(data, ncomp = c(1,2))
{
  # Do PCA
  pca <- PCA(data[,-which(colnames(data) %in% "TVC")], graph = F)
  summary(pca)
  
  # Create new column for TVC intervals
  quantiles <- quantile(data$TVC, c(.33, .66)) # Find quantiles for TVC values
  data$TVCinterval <- ifelse(data$TVC <= quantiles[1], 'low',
                             ifelse(data$TVC > quantiles[1] & data$TVC <= quantiles[2], 'med',
                                    ifelse(data$TVC > quantiles[2], 'high', '')))
  
  plot <- ggbiplot::ggbiplot(pca, choices = ncomp,
                             groups = data$TVCinterval, ellipse = TRUE, var.axes = FALSE) +
    theme_bw(base_family = "DMOFCB+AdvGulliv-R") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
          legend.text = element_text(size=15))
  
  return(plot)
}

# Plots intensities of two samples, first sample is colored green, second red
msi_compare_two_samples <- function(data, sample_one, sample_two)
{
  data <- as.data.frame(t(data)) # Transpose
  data <- data[-which(rownames(data) %in% "TVC"),] # Remove TVC row
  data$wave <- c(1:18) # Name the waves
  
  plot <- ggplot(data,aes(x=wave,y=data[,sample_one])) +
    geom_line(col="green") +
    geom_line(aes(x=wave,y=data[,sample_two]),col="red") +
    theme_bw(base_family = "DMOFCB+AdvGulliv-R") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust=0.5)) +
    xlab("Mean number at certain wavelengths") +
    ylab("Intensity") +
    scale_x_continuous(breaks = data$wave) +
    scale_fill_discrete(name="Spoilage level")
  
  return(plot)
}

# Compares intensities of two samples over wavelengths
ftir_compare_two_samples <- function(data, sample_one, sample_two)
{
  data <- data[rownames(data) %in% c(sample_one, sample_two),]
  data <- as.data.frame(t(data))
  data <- transform(data, wave=stringr::str_extract(rownames(data),"[0-9]+\\.+[0-9]*"))
  data$wave <- as.numeric(as.character(data$wave)) # Set as numeric
  difference <- (data[,1] - data[,2])[1]
  data <- data[-nrow(data), ]
  
  if(difference < 0) {
    fresh <- data[, 1]
    spoiled <- data[ ,2]
  } else {
    fresh <- data[, 2]
    spoiled <- data[ ,1]
  }
  
  plot <- ggplot(data=data, aes(x=data[,3])) +
    geom_line(aes(y=fresh, color="Fresh")) +
    geom_line(aes(y=spoiled, color="Spoiled")) +
    scale_color_manual(values=c("green","red")) +
    theme_bw(base_family = "DMOFCB+AdvGulliv-R") +
    labs(color="Level") +
    ylab("Wave") +
    xlab("Intensity")
  
  return(plot)
}

# Custom partioning function for train/test split, option 0=naive, 1=time, 2=temperature, 4=Kennard-stone
split <- function(data, option) 
{
  set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)) # Set seed
  
  if(option == 0) # Naive with caret TVC
  {
    train_index = caret::createDataPartition(y=data$TVC, p=0.7, list=FALSE, times=1)
    
    train <- data[train_index,]
    test <- data[-train_index,]
  }
  else if(option == 1) # Splitted equally with time, <70 is low, 70<time<140 is med, >140 is high
  {
    data$Time <- as.numeric(as.character((stringr::str_extract(rownames(data),"[0-9]*(?=h)"))))
    
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
             length(which(data$Timeinterval == 'high'))) * 0.7 + 6
    
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
    
    train <- train[,-which(colnames(train)%in%"Time")]
    test <- test[,-which(colnames(test)%in%"Time")]
  }
  else if(option == 2) # Split equally by temperature
  {
    data$Temperature <- as.numeric(as.character(stringr::str_extract(rownames(data),"[0-9]*((?=C)|(?=B)|(?=A))")))
    
    data[is.na(data)] <- 0
    
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
    
    train <- train[,-which(colnames(train)%in%"Temperature")]
    test <- test[,-which(colnames(test)%in%"Temperature")]
  }

  else if(option == 3) # Kennard-stone algorithm
  {
    ken <- kenStone(data, k=nrow(data)*0.7, metric="euclid")
    train <- data[ken$model,]
    test <- data[ken$test,]
  }
  
  return(list("train" = train, "test" = test))
}

# Custom predict functi10on which prints & returns prediction vs actual plot and displays RMSE/ACC for a single model
custompredict <- function(model, test)
{
  # Predict
  prediction <- stats::predict(model, test)
  
  # RMSE calc
  rmse <- Metrics::rmse(test$TVC, prediction)
  
  # Within 1 cfu accuracy calc
  difference <- as.data.frame(abs(test$TVC-prediction))
  accuracy <- sum(difference <= 1) / nrow(difference)
  
  x <- cbind(prediction,test$TVC)
  x <- as.data.frame(x[order(x[,2]),])

  # Plot
  title <- paste(model$method)
  colnames(x) <- c("Predicted", "Actual")
  plot <- ggplot(data=x, aes(x=Actual,y=Predicted)) +
    geom_point() +
    ggtitle(title) +
    ylim(c(0,9)) +
    xlim(c(0,9)) +
    geom_abline(col="blue") +
    geom_abline(intercept = 1, col="red") +
    geom_abline(intercept = -1, col="red") +
    geom_label(label=paste("Accuracy: ", round(accuracy*100, digits = 2),"%"), 
               x=7, y=3, size=5, family="DMOFCB+AdvGulliv-R",label.size = 0) +
    geom_label(label=paste("RMSE: ", round(rmse, digits = 4)), 
               x=7, y=3.5, size=5, family="DMOFCB+AdvGulliv-R", label.size = 0) +
    theme_bw(base_family = "DMOFCB+AdvGulliv-R") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust=0.5), legend.position = "none") +
    xlab(expression(paste("Sample Bacterial Count log"[10]," cfu g"^-1))) +
    ylab(expression(paste("Predicted Bacterial Count log"[10]," cfu g"^-1)))
  
  print(plot)
}

# Filter wavelengths from FTIR, list = waves to keep as integer ranges (i.e. 700:1000, 2000:2500)
keepwave <- function(data, list)
{
  # Transpose to filter wave-lengths
  data_T <- as.data.frame(t(data))
  # Add waves as a column
  data_T$wave <- round(as.numeric(stringr::str_extract(rownames(data_T[,-1]),"[0-9]+\\.+[0-9]*")),0)
  # Add temporary 0 to TVC wave
  data_T["TVC","wave"] <- 0
  
  list <- c(list, 0)
  # Filter wavelengths (include TVC)
  x <- data_T[data_T$wave %in% list, ]
  
  # Transpose back
  x <- as.data.frame(t(x))
  
  return(x[-nrow(x),])
}

# Set up a grid list for CV
gridList <- function()
{
  plsgrid <- expand.grid(ncomp = 1:40)
  
  svmRadialgrid <- expand.grid(sigma = 2^c(-25, -20, -15, -10, -5, 0), C = 2^c(0:5))
  
  svmLineargrid <- expand.grid(C = 2^c(0:5))
  
  knngrid <- expand.grid(k = 1:40)
  
  nnetgrid <- expand.grid(size = c(1:2), decay = seq(from = 0.1, to = 1, by = 0.1))
  
  larsgrid <- expand.grid(fraction = seq(from = 0, to = 1, by = 0.01))
  
  ridgegrid <- expand.grid(lambda=seq(0.00001,0.0001,0.00001))
  
  rfgrid <- expand.grid(mtry=c(2,4,6,12,24,48))
  
  # Create a list of grids
  gridlist <- list(plsgrid=plsgrid, svmRadialgrid=svmRadialgrid, svmLineargrid=svmLineargrid,
                   knngrid=knngrid, nnetgrid=nnetgrid, larsgrid=larsgrid,
                   rfgrid=rfgrid, nnetgrid=nnetgrid, ridgegrid=ridgegrid)
  
  return(gridlist)
}

# Splits the data depending on low and high parameters, plots the mean intensities of the split data over wavelengths
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

# Applies caret pre-process methods in methodlist to the data
# If pca pre-process is selected, pcaComp = components to keep
applypreprocess <- function(data, methodlist, pcaComp=0)
{
  TVC <- data$TVC
  data <- data[, -which(colnames(data)%in%"TVC")]
  
  if(pcaComp > 0) {
    pp <- preProcess(data, method = methodlist, pcaComp = pcaComp)
  } else {
    pp <- preProcess(data, method = methodlist)
  }
  
  data <- predict(pp, data)
  data$TVC <- TVC
  return(data)
}

# Creates a hyperspec object, type is either "FTIR" or else ("MSI")
createspc <- function(data, type="FTIR") # Removes TVC
{
  data <- data[,-which(colnames(data) %in% "TVC")]
  
  if(type == "FTIR") {
    waves <- as.numeric(stringr::str_extract(colnames(data),"[0-9]+\\.+[0-9]*"))
    spc <- new("hyperSpec", spc=data, wavelength=waves)
  }
  else {
    spc <- new("hyperSpec", spc=data, wavelength=as.numeric(colnames(data)))
  }
  
  return(spc)
}

# Applies spectral pre-process from the hyperspec package
# sgwindow = Savitzky-Golay filtering window, deriv = derivation factor, bin = binning factor
# type = "FTIR" or "VM", normalarea = normalise under area or a point
applyspectrapreprocess <- function(data, preproc, sgwindow=0, deriv=1, bin=4, type="FTIR", normalarea=F)
{
  TVC <- data$TVC
  
  if(preproc == "baseline") {
    spc <- createspc(data,type)
    bl <- spc.fit.poly.below(spc)
    spc <- spc - bl
    data <- as.data.frame(spc$spc)
    data$TVC <- TVC
  } else if(preproc == "snv") {
    data <- data[,-which(colnames(data) %in% "TVC")]
    data <- as.data.frame(standardNormalVariate(data)) 
    data$TVC <- TVC
  } else if(preproc == "snv-d") {
    data <- as.data.frame(detrend(X=spc$spc, wav=attributes(spc)$wavelength))
  } else if(preproc == "s-g"){
    data <- data[,-which(colnames(data) %in% "TVC")]
    data <- as.data.frame(savitzkyGolay(data, m=2,p=3,w=sgwindow))
    data$TVC <- TVC
  } else if(preproc == "deriv"){
    spc <- createspc(data,"VM")
    data <- as.data.frame(t(diff(t(spc$spc), differences=deriv)))
    data$TVC <- TVC
  } else if(preproc == "normal"){
    spc <- createspc(data,type)
    if(normalarea){
      spc <- spc/ rowMeans(spc[, ,1500~1700])
    } else if(type=="FTIR" & !normalarea){
      factors <- 1/apply(spc[, ,1500~1700],1,mean)
      spc <- sweep(spc,1,factors,"*")
    } else{
      spc <- sweep(spc,1,mean,"/")
    }
    data <- as.data.frame(spc$spc)
    data$TVC <- TVC
  } else if(preproc == "bin"){
    spc <- createspc(data,type)
    spc <- spc.bin(spc,bin)
    data <- as.data.frame(spc$spc)
    data$TVC <- TVC
  }
  
  return(data)
}

# Predicts and calculates performances of multiple models on a test data
# models = caret model list, test = test data
# plot = T would plot each model's predicted vs actual with custompredict()
multiplemodelpredict <- function(models, test, plot=F)
{
  rmseResults <- data.frame(matrix(nrow=length(models),ncol=5))
  colnames(rmseResults) <- c("RMSE", "MAE","R2", "Accuracy","name")
  for(i in 1:length(models))
  {
    model <- models[[i]]
    
    # RMSE
    prediction <- predict(model, test)
    rmse <- Metrics::rmse(test$TVC, prediction)
    
    # R2
    r2 <- MLmetrics::R2_Score(test$TVC,prediction)
    
    # MAE
    mae <- Metrics::mae(test$TVC,prediction) # Calculate mae for the current model and iter

    # Within 1 cfu accuracy
    difference <- as.data.frame(abs(test$TVC-prediction))
    accuracy <- sum(difference <= 1) / nrow(difference)
    
    if(plot == T)
    {
      custompredict(model, test)
    }
    rmseResults[i,] <- list(rmse, mae, r2, accuracy, model$method)
  }
  
  rownames(rmseResults) <- rmseResults[,5]
  modelmethods <- rmseResults$name
  rmseResults <- round(rmseResults[,1:4],4)
  rmseResults$name <- modelmethods
  return(rmseResults)
}

# A.k.a monte-carlo CV, does multiple iterations with multiple models, multiple splits and calculate average performance metrics
multipleiterations <- function(data, iters, parallelcores = 3)
{
  x <- iters
  methodlist <- c("pls", "svmLinear", "svmRadial","rf","knn", "pcr","lars", "ridge", "nnet")
  
  rmseResults <- data.frame(matrix(nrow=length(methodlist), ncol=1)) # DF to store RMSE results
  rownames(rmseResults) <- methodlist
  
  maeResults <- data.frame(matrix(nrow=length(methodlist), ncol=1)) # DF to store RMSE results
  rownames(maeResults) <- methodlist
  
  r2Results <- data.frame(matrix(nrow=length(methodlist), ncol=1)) # DF to store RMSE results
  rownames(r2Results) <- methodlist
  
  accResults <- data.frame(matrix(NA, nrow=length(methodlist), ncol=1)) # DF to save acc results
  rownames(accResults) <- methodlist
  
  max <- 0
  min <- 50
  cl <- makePSOCKcluster(parallelcores) # Parelell processing with 3 cores
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
    
    gridlist <- gridList()
    
    fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                               savePredictions = "final"
                               ,allowParallel = T)
    
    models <-  caretEnsemble::caretList(TVC~., data=train, trControl = fitControl, metric = "RMSE", continue_on_fail = T
                                        ,tuneList = list(
                                          pls=caretModelSpec(method="pls", tuneGrid=gridlist$plsgrid)
                                          ,svmLinear=caretModelSpec(method="svmLinear", tuneGrid=gridlist$svmLineargrid)
                                          ,svmRadial=caretModelSpec(method="svmRadial", tuneGrid=gridlist$svmRadialgrid)
                                          ,rf=caretModelSpec(method="rf", tuneGrid=gridlist$rfgrid, ntree=500)
                                          ,knn=caretModelSpec(method="knn", tuneGrid=gridlist$knngrid)
                                          ,pcr=caretModelSpec(method="pcr", tuneGrid=gridlist$plsgrid)
                                          ,lars=caretModelSpec(method="lars", tuneGrid=gridlist$larsgrid)
                                          ,ridge=caretModelSpec(method="ridge", tuneGrid=gridlist$ridgegrid)
                                          ,nnet=caretModelSpec(method="nnet", tuneGrid=gridlist$nnetgrid, linout = TRUE, maxit=1000)
                                        )
    )
    
    # Predictions for test set
    predictionsDF <- as.data.frame(predict(models, newdata = test))
    
    # Calculate prediction difference
    predDiff <- abs(predictionsDF - test$TVC)
    
    # Calculate rmse for each model
    rmseVector <- c()
    maeVector <- c()
    r2Vector <- c()
    for(k in 1:length(models))
    {
      # RMSE
      rmse <- Metrics::rmse(test$TVC,predictionsDF[,k]) # Calculate rmse for the current model and iter
      rmseVector <- c(rmseVector, rmse) # Add to the rmselist for all rmse on this iter
      
      # MAE
      mae <- Metrics::mae(test$TVC,predictionsDF[,k]) # Calculate mae for the current model and iter
      maeVector <- c(maeVector, mae) # Add to the rmselist for all mae on this iter
      
      # R-squared
      r2 <- MLmetrics::R2_Score(test$TVC,predictionsDF[,k]) # Calculate rsquared for the current model and iter
      r2Vector <- c(r2Vector, r2) # Add to the rmselist for all rsquared on this iter
      
      # ACCURACY
      # Check how many predictions are within 1 cfu/log
      table <- table(predDiff[,k] < 1)
      if(is.na(table[2])){
        pred <- 1
      } else {
        pred <- as.numeric(table[2]) / as.numeric(table[1] + table[2])
      }
      accResults[k,i] <- pred # Store current model+iter accuracy on the DF
    }
    
    # Add RMSEs together for this iteration
    rmseResults <- cbind(rmseResults,rmseVector)
    # Add MAEs together for this iteration
    maeResults <- cbind(maeResults,maeVector)
    # Add R-squareds together for this iteration
    r2Results <- cbind(r2Results,r2Vector)
    
    print(paste("This iteration took:", round(Sys.time() - starttime,3)))
  }
  stopCluster(cl) # Stop parellel processing
  registerDoSEQ()
  
  rmseResults <- rmseResults[,-1, drop=F]
  maeResults <- maeResults[,-1, drop=F]
  r2Results <- r2Results[,-1, drop=F]
  
  colnames(rmseResults) <- 1:x
  colnames(maeResults) <- 1:x
  colnames(r2Results) <- 1:x
  colnames(accResults) <- 1:x 
  
  rmseResults <- round(rmseResults,4)
  maeResults <- round(maeResults,4)
  r2Results <- round(r2Results,4)
  accResults <- round(accResults,4)
  
  returnlist <- list(rmseResults=rmseResults, maeResults=maeResults, 
                     r2Results=r2Results, accResults=accResults)
  
  return(returnlist)
}

# Combines performance metrics(recieved from multipleiterations())into one data frame
getfinalresults <- function(rmseResults, maeResults, r2Results, accResults)
{
  # RMSE results
  rmseResults$mean <- apply(rmseResults,1,FUN=mean)
  
  # MAE results
  maeResults$mean <- apply(maeResults,1,FUN=mean)
  
  # R2 results
  r2Results$mean <- apply(r2Results,1,FUN=mean)
  
  # Accuracy results
  accResults$mean <- apply(accResults,1,FUN=mean)

  # Mean combined final results
  finalresults <- data.frame(matrix(nrow=nrow(rmseResults),ncol=4)) # Create results DF
  rownames(finalresults) <- rownames(accResults)
  colnames(finalresults) <- c("RMSE", "MAE","R2","Accuracy")
  
  finalresults[,1] <-  round(rmseResults$mean,4)
  finalresults[,2] <- round(maeResults$mean,4)
  finalresults[,3] <- round(r2Results$mean,4)
  finalresults[,4] <- round(accResults$mean,4)
  finalresults$name <- rownames(finalresults)
  
  return(finalresults)
}

# Does multiple iterations with a model on a data frame and averages variable importance over iterations
# Only pls and nnet supported
variableimportancetest <- function(data, iters, model, parallelcores = 3)
{
  cl <- makePSOCKcluster(parallelcores) # Set parallel processing
  registerDoParallel(cl)
  result <- data.frame(matrix(ncol=1,nrow=ncol(data)-1))
  
  for(i in 1:iters)
  {
    print(paste("Iteration",i))
    starttime <- Sys.time()
    
    # Split to train and test
    return <- split(data, 0) 
    train <- return$train
    test <- return$test
    
    set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
    
    # Set grids
    nnetgrid <- expand.grid(size = c(1,2), decay = seq(0.1, 1, 0.1))
    plsgrid <- expand.grid(ncomp=seq(1,40,1))
    
    # Set CV method
    fitControl <- trainControl(method = "repeatedcv", number = 5, repeats=3
                               ,savePredictions = "final"
                               ,allowParallel = T)
    
    # Run the model
    if(model == "pls"){
      models <-  caretEnsemble::caretList(TVC~., data=data, trControl = fitControl, metric = "RMSE"
                                          ,tuneList = list(
                                            pls=caretModelSpec(method="pls", tuneGrid=plsgrid)
                                          )
      )
    } else if(model == "nnet"){
      models <-  caretEnsemble::caretList(TVC~., data=data, trControl = fitControl, metric = "RMSE"
                                          ,tuneList = list(
                                            nnet=caretModelSpec(method="nnet", tuneGrid=nnetgrid, 
                                                                linout = TRUE, maxit=1000, MaxNWts=7000)
                                          )
      )
    }
    
    # Get variable importance
    imp <- as.data.frame(varImp(models[[1]])$importance)
    result <- cbind(result,imp)
    
    print(paste("Iteration took", Sys.time()-starttime))
  }
  
  stopCluster(cl)
  registerDoSEQ()
  result <- result[,-1]
  
  # Put mean as extra col
  result$mean <- apply(result,1,FUN=mean)
  
  # Extract mean importance corresponding to integer waves
  impmean <- result[,"mean",drop=F]
  impmean <- transform(impmean, wave=as.numeric(stringr::str_extract(rownames(impmean),"[0-9]+\\.+[0-9]*")))
  impmean <- impmean[order(-impmean$mean), ]
  impmean$wave <- round(impmean$wave,0)
  
  return(list(iterresult = result, meanimp = impmean))
}

# Plots a heatmap with multiple results from getfinalresults()
# resultlist = list(finalresults1, finalresult2,..), names = labels for each result (B1, B2..)
# excludeModel = exclude these models by name ("knn", "rf"..)
# includeModel = include these models by name ("knn", "rf"..)
finalheatmap <- function(resultlist, names, excludeModel = c(), includeModel = c(), limits=c(50,90))
{
  getperformance <- function(resultlist, names, perfname)
  {
    for(i in 1:length(resultlist))
    {
      resultlist[[i]] <- resultlist[[i]][,perfname, drop=F]
      resultlist[[i]]$model <- rownames(resultlist[[i]])
      resultlist[[i]]$batch <- rep(names[i],nrow(resultlist[[i]]))
    }
    
    return(resultlist)
  }
  
  rmselist <- getperformance(resultlist, names, "RMSE")
  maelist <- getperformance(resultlist, names, "MAE")
  rsquaredlist <- getperformance(resultlist, names, "R2")
  accuracylist <- getperformance(resultlist, names, "Accuracy")
  
  forheatmapRMSE <- ldply(rmselist, data.frame)
  forheatmapMAE <- ldply(maelist, data.frame)
  forheatmapR2 <- ldply(rsquaredlist, data.frame)
  forheatmapAcc <- ldply(accuracylist, data.frame)
  forheatmapAcc$Accuracy <- forheatmapAcc$Accuracy * 100
  
  for(model in excludeModel)
  {
    forheatmapRMSE <- forheatmapRMSE[which(!forheatmapRMSE$model %in% excludeModel),]
    forheatmapMAE <- forheatmapMAE[which(!forheatmapMAE$model %in% excludeModel),]
    forheatmapR2 <- forheatmapR2[which(!forheatmapR2$model %in% excludeModel),]
    forheatmapAcc <- forheatmapAcc[which(!forheatmapAcc$model %in% excludeModel),]
  }
  
  for(model in includeModel)
  {
    forheatmapRMSE <- forheatmapRMSE[which(forheatmapRMSE$model %in% includeModel),]
    forheatmapR2 <- forheatmapR2[which(forheatmapR2$model %in% excludeModel),]
    forheatmapAcc <- forheatmapAcc[which(forheatmapAcc$model %in% excludeModel),]
    forheatmapAcc <- forheatmapAcc[which(forheatmapAcc$model %in% includeModel),]
  }
  
  forheatmapAcc <- forheatmapAcc[order(-forheatmapAcc$Accuracy),]
  plot <- ggplot(forheatmapAcc, aes(batch, reorder(model,Accuracy))) +
    geom_tile(aes(fill=Accuracy)) +
    scale_fill_gradientn(colors=colorRampPalette(c("tomato2","green"))(3), limits=limits, oob=squish) +
    theme_bw(base_size = 20) +
    labs(x="", y="") +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    geom_text(data=forheatmapAcc,aes(batch,model, label=paste("Acc:",Accuracy,"%")), 
              vjust=-1, hjust=1, size=4.5, family="DMOFCB+AdvGulliv-R") +
    geom_text(data=forheatmapRMSE,aes(batch,model,label=paste("RMSE:",RMSE)),
              vjust=2, hjust=1, size=4.5, family="DMOFCB+AdvGulliv-R") +
    geom_text(data=forheatmapMAE,aes(batch,model,label=paste("MAE:",MAE)),
              vjust=-1, hjust=-0.2, size=4.5, family="DMOFCB+AdvGulliv-R") +
    geom_text(data=forheatmapR2,aes(batch,model,label=paste("R2:",R2)),
              vjust=2, hjust=-0.2, size=4.5, family="DMOFCB+AdvGulliv-R") +
    theme(legend.position = "none", text=element_text(family="DMOFCB+AdvGulliv-R")) 
    

  print(plot)
  
  return(plot)
}

#----(NOT USED IN PAPER)----
# Plots the importance of variables over wavelengths
plotvariableimportance <- function(model)
{
  imp <- as.data.frame(varImp(model)$importance)
  imp <- transform(imp, wave=as.numeric(stringr::str_extract(rownames(imp),"[0-9]+\\.+[0-9]*")))
  imp <- imp[order(imp$wave), ]
  ggplot(data=imp, aes(x=wave, y=Overall)) +
    geom_point()
}

# Filter waves according to RF variable importance, returns filtered data
# topx = Return FTIR data with topx waves. threshold = Return FTIR data with filtering values under threshold
# if topx and threshold is 0, returns the importance list
findandkeepwavesRF <- function(data, topx=0, threshold=0) # no topx argument returns the list
{
  cl <- makePSOCKcluster(3)
  registerDoParallel(cl)
  # Find optimal mtry value (RF parameter)
  tune <- tuneRF(data[,-1], data[,1], trace = FALSE, plot = FALSE)
  mtry <- as.numeric(row.names(tune)[(which(tune[,2]==min(tune[,2])))])
  
  # Train model
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

# Plots std of wavelengths
plotSTD <- function(data)
{
  data_T <- as.data.frame(t(data)) # Transpose
  data_T <- data_T[-which(rownames(data_T) %in% "TVC"),] # Remove TVC row
  data_T <- transform(data_T,SD=apply(data_T,1,sd)) # SD of every row
  data_T$waves <- as.numeric(stringr::str_extract(rownames(data_T),"[0-9]+\\.+[0-9]*")) # Create waves as column
  
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

# Filters the data according to the std values
# std = std data frame, threshold = lowest std allowed 
filterwavesfromSTD <- function(data, std, threshold)
{
  std <- std[which(std$SD > threshold),]
  names <- c(rownames(std),'TVC')
  data <- data[, colnames(data) %in% names]
  
  return(data)
}

# Calculates and plots the linear model coefficients of a data
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

# Filters waves in FTIR by contribution to PCA
# contribution = contribution data frame calculated elsewhere, threshold = values under threshold will be filtered
filterwavesbyPCAcontrib <- function(data, contribution, threshold)
{
  filterlist <- rownames(contribution[contribution$sum > threshold,]) # Filter contribs below threshold
  filterlist <- str_extract(filterlist, "[0-9]*(?=\\.)") # Convert row names to wave numbers
  data <- keepwave(data, filterlist) # Filter the data 
  
  return(data)
}

# Tries to find patterns of train/test split from a list of RMSE results for multiple models
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