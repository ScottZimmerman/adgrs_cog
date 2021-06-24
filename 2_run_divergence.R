rm(list=ls())
library("ggplot2")
library("splines")
library("stringr")
library("dplyr")
library("gridExtra")
library("MASS")
library("butcher")

#----------------------------------------------
# Description
#----------------------------------------------
# Calculate age-of-divergence between low- and
#   high-ADGRS cognitive phenome-wide ADGRS
#   association study
#
# By Scott Zimmerman (scott.zimmerman@ucsf.edu)

#----------------------------------------------
# Analysis parameters
#----------------------------------------------
#This script is run twice, setting sn to different values:
# "withAPOE": Main analyses
# "withoutAPOE": Sensitivity analysis omitting APOE alleles from the AD-GRS
sn <- "withAPOE"

#We use string replacement below to make the appropriate functional form for each phenotype,
# using covsStrings and modelStrings as inputs (see "makeModelText" function)
covsStrings <- list(
  covs = paste0("sex_0_0+genoChip_0_0+",paste0("principalComponents_0_",1:10,collapse="+"),"{+ASSESSMENTCENTER}","{+PRACEFFECT}")
)
modelStrings <- list(
  noCovs = "{PHENOTYPE}~I({AGEVARIABLE}^2)+{AGEVARIABLE}+1+{GRSVARIABLE}:gt_t_c:I(({AGEVARIABLE}-{CUTOFF})^3)",
  covs = "{PHENOTYPE}~I({AGEVARIABLE}^2)+{AGEVARIABLE}+{COVSSTRING}+1+{GRSVARIABLE}:gt_t_c:I(({AGEVARIABLE}-{CUTOFF})^3)"
)

#The method to fit the data uses V-fold cross validation to select the threshold cutoff
p <- list(
  modelInfoFile = "modelInfo.csv",
  outputDir = sn,
  cutoff = list( #Range of cutoff ages to try
    min=40,
    max=70,
    step=1
  ),
  folds=list(
    n=10
  ),
  errorType = "L2", #Type of error (L2 or L1)
  grsVariable = "ADGRS_z",
  ordinalAsLinear=TRUE, #Error with prediction if run with FALSE: multicollinearity?
  covsString = covsStrings$covs, #This is left over from trying multiple covariate options
  model = modelStrings$covs      #This is left over from a sensitivity analysis omitting covariates ("noCovs")
)

#Plot predictions while running (not final plots, but useful for debugging, runs faster if FALSE)
plotOn <- TRUE

#----------------------------------------------
# Directories
#----------------------------------------------
dirs <- list(
  output = file.path('Results', 'Divergence', p$outputDir),
  data = file.path("Data")
)
dirs$plots <- file.path(dirs$output,"Plots")

# Create new directories
if(!dir.exists(dirs$output)){
  dir.create(dirs$output)
}
dir.create(dirs$plots)

#----------------------------------------------
# Import data
#----------------------------------------------
data <- read.csv(file.path(dirs$data,sn,"PHEWAS_cog_analysis_data.csv"))

modelInfo <- read.csv(file.path(dirs$data,p$modelInfoFile), #File that describes all of the outcome variables to loop over
                      stringsAsFactors = FALSE)
rownames(modelInfo) <- modelInfo[,2]
modelInfo$variable <- NULL
modelInfo$ageMin <- NA
modelInfo$ageMax <- NA
for(pheno in rownames(modelInfo)){
  d <- data[!is.na(data[[pheno]])& !is.na(data[[p$grsVariable]]),]
  if(modelInfo[pheno,"regressionType"] != "SKIP"){
    ageVar <- modelInfo[pheno,"ageVar"]
    modelInfo[pheno,"ageMin"] <- min(d[[ageVar]], na.rm=TRUE)
    modelInfo[pheno,"ageMax"] <- max(d[[ageVar]], na.rm=TRUE)
  }
}

#Outcome variables
phenotypeNames <- rownames(modelInfo)

#Keep track of some information
trackingCols <- c("ageVar","regressionType","sampleSize","skipped","success","formula")
trackingData <- data.frame(matrix(nrow=nrow(modelInfo),ncol=length(trackingCols)))
rownames(trackingData) <- rownames(modelInfo)
colnames(trackingData) <- trackingCols

#-------------------------------
# Helper functions
#-------------------------------
#Function to set up data frame (cv_df) for storing cross-validation results
make_cv_df <- function(pheno,p){
  minCutoff <- max(modelInfo[pheno,"ageMin"],p$cutoff$min, na.rm=TRUE)
  maxCutoff <- min(modelInfo[pheno,"ageMax"],p$cutoff$max, na.rm=TRUE)
  cv_df <- data.frame(cutoff=seq(minCutoff,
                                 maxCutoff,
                                 p$cutoff$step),
                      err_avg=NA)
  cv_folds <- data.frame(matrix(nrow=nrow(cv_df),ncol=p$folds$n))
  colnames(cv_folds) <- paste0("err_",1:p$folds$n)
  cv_df <- cbind(cv_df,cv_folds)
  return(cv_df)
}

#Function to create the regression formula by
# replacing {TEXT} with variable names or functions
makeModelText <- function(pheno,cutoff,ageVar,modelInfo,p){
  covsString <- p$covsString

  #Add assessment center only if it is an in-person test
  splitPheno <- strsplit(pheno,"_")[[1]]
  ac_or_ol <- splitPheno[[2]]
  instance <- splitPheno[[3]]
  if(ac_or_ol == "ac" & instance != 1){
    #There is no variation in assessment center in instance 1
    acString <- paste0("+assessmentCenter_",instance,"_0")
  }else{
    acString <- ""
  }

  #For online tests, add indicator for whether they took an in-person version of the test (practice effect)
  if(ac_or_ol == "ol"){
    olString <- "+prevTest"
  }else{
    olString <- ""
  }

  #Do string replacements
  covsString <- paste0("(",covsString,")")
  modelText <- str_replace_all(p$model,"\\{CUTOFF\\}",as.character(cutoff))
  modelText <- str_replace_all(modelText,"\\{PHENOTYPE\\}",as.character(pheno))
  modelText <- str_replace_all(modelText,"\\{COVSSTRING\\}",as.character(covsString))
  modelText <- str_replace_all(modelText,"\\{AGEVARIABLE\\}",as.character(ageVar))
  modelText <- str_replace_all(modelText,"\\{GRSVARIABLE\\}",as.character(p$grsVariable))
  modelText <- str_replace_all(modelText,"\\{\\+ASSESSMENTCENTER\\}",acString)
  modelText <- str_replace_all(modelText,"\\{\\+PRACEFFECT\\}",olString)

  return(modelText)
}

#Function to run the appropriate type of regression
run_model <- function(modelText,d,modelInfo,p){
  pheno <- strsplit(modelText,"~")[[1]][[1]]
  regressionType <- modelInfo[pheno,"regressionType"]

  if(regressionType == "linear"){
    model <- lm(modelText,
                data=d,
                model=FALSE,
                x=FALSE,
                y=FALSE)

  }else if(regressionType == "ordinal"){
    if(p$ordinalAsLinear){
      model <- lm(modelText,
                  data=d,
                  model=FALSE,
                  x=FALSE,
                  y=FALSE)
    }else{
      model <- polr(modelText,
                    d,
                    Hess=TRUE,
                    model=FALSE)
    }

  }else if(regressionType == "binary"){
    model <-  glm(modelText,
                 d,
                 family = "binomial",
                 model=FALSE,
                 x=FALSE,
                 y=FALSE)
  }
  return(model)
}

#Function to do V-fold cross-validation
do_xVal <- function(d,pheno,modelInfo,p){
  regressionType <- modelInfo[pheno,"regressionType"]
  print(paste0("Running model (",regressionType,"): ",pheno))
  if(regressionType == "ordinal"){
    if(p$ordinalAsLinear != TRUE){
      d[[pheno]] <- factor(d[[pheno]])
    }
  }
  n <- nrow(d)
  cv_df <- make_cv_df(pheno,p)
  nPerFold <- floor(n/p$folds$n)

  #Add practice effect variable to data
  splitPheno <- strsplit(pheno,"_")[[1]]
  ac_or_ol <- splitPheno[[2]]
  if(ac_or_ol == "ol"){
    pracEffVar <- modelInfo[pheno,"pracEffVar"]
    data$prevTest <- is.na(data[[pracEffVar]])
  }

  #Loop through the folds and calculate the error
  # Use the same folds for each cutoff option tested
  for(f_id in 1:p$folds$n){
    print(paste0(round(100*f_id/p$folds$n,1),"%")) #Progress tracker

    #Set up data
    startIndex <- nPerFold*(f_id-1)+1
    test_indices <- startIndex:(startIndex+nPerFold-1)
    train_indices <- setdiff(1:n,test_indices)

    ageVar <- modelInfo[pheno,"ageVar"]
    for(i in 1:nrow(cv_df)){

      #Set up model
      cutoff <- cv_df[i,"cutoff"]
      d$t <- d[[as.character(ageVar)]]
      ageVar <- "t"
      d$gt_t_c <- as.numeric(d[[as.character(ageVar)]] > cutoff)

      #Split into train and test data sets
      d_train <- d[train_indices,]
      d_test <- d[test_indices,]

      modelText <- makeModelText(pheno,cutoff,ageVar,modelInfo,p)

      #Train model
      model <- run_model(modelText,d_train,modelInfo,p)

      #Predict from the model on the test data
      y_pred <- predict(model, newdata=d_test)

      #Calculate error
      y_true <- d_test[[pheno]]
      if(p$errorType == "L2"){
        err <- sum((y_true-y_pred)^2,na.rm=TRUE)
      }else if (p$errorType == "L1"){
        err <- sum(abs(y_true-y_pred),na.rm=TRUE)
      }else{
        stop("p$errorType must be 'L1' or 'L2'")
      }

      #assign error for this fold to cv_df
      cv_df[i,paste0("err_",f_id)] <- err
    }
  }

  #Calculate mean error across folds for each cutoff value
  cv_df$err_avg <- apply(cv_df[,paste0("err_",1:p$folds$n)],1,mean)
  print(cv_df[,c("cutoff","err_avg")])

  #Which cutoff has minimum cv error?
  best_cutoff <- cv_df[which.min(cv_df$err_avg),"cutoff"]

  #Run the model on the full data with the best cutoff
  modelText <- makeModelText(pheno,best_cutoff,ageVar,modelInfo,p)

  model <- run_model(modelText,d,modelInfo,p)

  #Remove some info from model to save space
  model <- butcher::axe_env(model)

  return(list(cv_df=cv_df,model=model,modelText=modelText))
}

#-------------------------------
# Analysis (Saving thoroughout)
#-------------------------------
saveRDS(p,file.path(dirs$output,"parameters.RDS"))
for(pheno in phenotypeNames){
  regressionType <- modelInfo[pheno,"regressionType"]
  trackingData[pheno,"success"] <- 1

  if(regressionType != "SKIP"){
    tryCatch({
      #Subset the data to only obs with the phenotype observed
      d <- data[!(is.na(data[[pheno]])),]

      #Calculate age of divergence using cross-validation
      # to choose the best-fitting model
      results <- do_xVal(d,pheno,modelInfo,p)

      #plot CV error
      g <- ggplot(results$cv_df,aes(x=cutoff,y=err_avg))+
        geom_line() +
        xlab("Threshold") +
        ylab("Mean CV Error") +
        ggtitle(paste0("CV Error - ",pheno))
      ggsave(file.path(dirs$plots,paste0("CVerror_",pheno,".png")),g)
      results$cv_ggplot <- g

      #Save outputs
      saveRDS(results,file.path(dirs$output,paste0(pheno,".RDS")))
      write.csv(results$cv_df,file.path(dirs$output,paste0(pheno,"_cv_df.csv")))

      #Record model information
      ageVar <- as.character(modelInfo[pheno,"ageVar"])
      trackingData[pheno,"formula"] <- results$modelText
      trackingData[pheno,"ageVar"] <- ageVar
      trackingData[pheno,"sampleSize"] <- nrow(d)
      trackingData[pheno,"regressionType"] <- regressionType
      trackingData[pheno,"skipped"] <- 0

      if(plotOn){
        #Plot predictions
        # This is not the final plot,
        # but just used to check things look ok while it's running
        d <- d[1:1000,]
        cutoff <- results$cv_df[which.min(results$cv_df$err_avg),"cutoff"]
        d$t <- d[[ageVar]]
        d$gt_t_c <- as.numeric(d$t>=cutoff)
        if(p$grsVariable=="ADGRS_z"){
          d[[p$grsVariable]] <- (rbinom(nrow(d),1,0.5)-.5)*4
        }
        d$pred_y <- predict(results$model,newdata=d)
        d$ADGRS_f <- factor(d[[p$grsVariable]])
        g_pred <- ggplot(d,aes(x=t,y=pred_y,color=ADGRS_f))+
          geom_line(size=0.5)+
          ggtitle("Prediction")
        ggsave(file.path(dirs$plots,paste0("pred_",pheno,".png")),g_pred)
      }

    },error= function(e){
      trackingData[pheno,"success"] <- 0
      print(e)
    })
  }else{
    trackingData[pheno,"regressionType"] <- regressionType
    trackingData[pheno,"skipped"] <- 1
  }
}

#Save tracking info
write.csv(trackingData,file.path(dirs$output,"2_run_divergence_tracking.csv"))