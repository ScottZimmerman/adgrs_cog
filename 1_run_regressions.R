rm(list=ls())
library("MASS")
library("butcher")

#----------------------------------------------
# Description
#----------------------------------------------
# Main effect and interaction analyses for
#  cognitive phenome-wide ADGRS association study
#
# By Scott Zimmerman (scott.zimmerman@ucsf.edu)

#----------------------------------------------
# Analysis parameters
#----------------------------------------------
#This script is run twice, setting sn to different values:
# "withAPOE": Main analyses
# "withoutAPOE": Sensitivity analysis omitting APOE alleles from the AD-GRS

sn <- "withAPOE"

p <- list( #parameters
  inScenarioName = sn,
  outScenarioName = sn,
  dataFile = "PHEWAS_cog_analysis_data.csv",
  modelInfoFile = "modelInfo.csv",
  centerAge = 40,
  covariateNames = c("sex_0_0",
                      "genoChip_0_0",
                      paste0("principalComponents_0_",1:10)),
  ordinalAsLinear=TRUE #Treat outcome variables of type "ordinal" in modelInfo.csv file as continuous covariates?
)

#----------------------------------------------
# Import data and set paths
#----------------------------------------------
#Directories
dirs <- list(
  data = file.path("Data"),
  out = file.path("Results","Main")
)

dirs$inScenario <- file.path(dirs$data,p$inScenarioName)
dirs$outScenario <- file.path(dirs$out,p$outScenarioName)
dir.create(dirs$outScenario)

#Model details for each phenotype variable
modelInfo <- read.csv(file.path(dirs$data,p$modelInfoFile), #File that describes all of the outcome variables to loop over
                      stringsAsFactors = FALSE)
rownames(modelInfo) <- modelInfo[,2]
modelInfo$variable <- NULL

#Analysis data
data <- read.csv(file.path(dirs$inScenario,p$dataFile))

#Check for missing phenotype variables
for(i in rownames(modelInfo)){
  if(!(i %in% colnames(data))){
    print(paste0("Missing variable: ",i))
  }
}

#Keep track of some information for each of the regressions
trackingCols <- c("ageVar","regressionType","sampleSize","skipped","success","formula")
trackingData <- data.frame(matrix(nrow=nrow(modelInfo),ncol=length(trackingCols)))
rownames(trackingData) <- rownames(modelInfo)
colnames(trackingData) <- trackingCols

#----------------------------------------------
# Analysis-specific data cleaning
#----------------------------------------------
#Center all age vars (different cognitive tests have different age variables)
ageVars <- c()
for(pheno in rownames(modelInfo)){
  ageVar <- modelInfo[pheno,"ageVar"]
  ageVars <- unique(c(ageVars,ageVar))
}
for(ageVar in ageVars){
  data[,ageVar] <- data[,ageVar] - p$centerAge
}

#----------------------------------------------
# Analysis
#----------------------------------------------
# Main loop: Run the regressions specified in modelInfo
#  for each phenotype as the outcome
models <- list()
for(pheno in rownames(modelInfo)){
  tryCatch({

    print(paste0("Running ", pheno))

    #Get the phenotype-specific information
    ageVar <- modelInfo[pheno,"ageVar"]
    pracEffVar <- modelInfo[pheno,"pracEffVar"]

    regressionType <- modelInfo[pheno,"regressionType"]
    splitPheno <- strsplit(pheno,"_")[[1]]
    ac_or_ol <- splitPheno[[2]] #Assessment center or online
    instance <- splitPheno[[3]] #Wave of the study

    #Construct the formula
    if(ac_or_ol == "ac" & instance != 1){
      acString <- paste0("assessmentCenter_",instance,"_0+")
    }else{
      ##There is no variation in in-person assessment center in instance 1
      ## so do not include it
      acString <- ""
    }

    #For online tests, add indicator for whether they took an in-person version of the test (practice effect)
    if(ac_or_ol == "ol"){
      data$prevTest <- is.na(data[[pracEffVar]])
      olString <- "prevTest"
    }else{
      olString <- ""
    }

    formula <- paste0(pheno,"~",ageVar,"*ADGRS_z+",acString,olString,paste0(p$covariateNames,collapse="+"))

    #Subset the data to only obs with the phenotype observed
    d <- data[!(is.na(data[[pheno]])),]

    #Record model information
    trackingData[pheno,"formula"] <- formula
    trackingData[pheno,"ageVar"] <- ageVar
    trackingData[pheno,"sampleSize"] <- nrow(d)
    trackingData[pheno,"regressionType"] <- regressionType
    trackingData[pheno,"skipped"] <- 0

    #Run the appropriate type of regression for the outcome
    if(regressionType == "linear" | (regressionType=="ordinal" & p$ordinalAsLinear==TRUE)){
      models[[pheno]] <- lm(formula,
                            data=d,
                            model=FALSE,
                            x=FALSE,
                            y=FALSE)

    }else if(regressionType == "ordinal"){
      d[[pheno]] <- factor(d[[pheno]])
      models[[pheno]] <- polr(formula,
                              d,
                              Hess=TRUE,
                              model=FALSE)

    }else if(regressionType == "binary"){
      vals <- sort(unique(d[[pheno]]))

      #Run checks that the variable is actually binary
      if(length(vals)!=2){
        print("Variable is not binary")

      }else if(!all.equal(vals, c(0,1))){
        print("Need to clean the variables so they're 0,1")

      }else{
        models[[pheno]] <- glm(formula,
                               d,
                               family = "binomial",
                               model=FALSE,
                               x=FALSE,
                               y=FALSE)
      }
    }else{
      trackingData[pheno,"skipped"] <- 1
    }

    #Record if success
    trackingData[pheno,"success"] <- 1

  },error=function(e){
    trackingData[pheno,"success"] <- 0
  }) ##End trycatch
} ##End main loop

#----------------------------------------------
# Save results
#----------------------------------------------
#Remove unnecessary variables from model objects to save disk space
for(modelName in names(models)){
  models[[modelName]] <- butcher::axe_env(models[[modelName]])
}

#Outputs
saveRDS(models, file=file.path(dirs$outScenario,"models.rds"))
saveRDS(p, file=file.path(dirs$outScenario,"params.rds"))
write.csv(trackingData,file.path(dirs$outScenario,"1_run_regressions_tracking.csv"))
write.csv(modelInfo,file.path(dirs$outScenario,p$modelInfoFile))