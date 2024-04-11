# Generates summary file for evaluate parameters
#USAGE: Rscript evaluate_parameters/summaryEvaluate.r [tools_params_file] [datasets_params_file]
source("utils/cnvStats.r")  # Load utils functions
suppressPackageStartupMessages(library(yaml))
library(methods)
library("optparse")
library("dplyr")

#Get parameters----
## load tools and datasets from yaml files----
args <- commandArgs(TRUE)
cat("\n"); print(args)
if(length(args)>0) {
  toolsParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  toolsParamsFile <- "tools.yaml"
  datasetsParamsFile <- "datasets.yaml"
}

setwd("evaluate_parameters")
#Load the parameters file
tools <- yaml.load_file(toolsParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# Run summary evaluate for all datasets and tools----
for (dName in names(datasets)) {
  dataset <- datasets[[dName]]
  if (dataset$include){
    ## Calculate metrics----
    for (alg in names(tools)) {
      useAlg <- tools[[alg]]
      if (useAlg){
        params <- list.dirs(paste0(getwd(), "/evaluate_parameters/", alg, "/", dName ), recursive=FALSE) # listdirectories of the parameters to evaluate
        for (param in params){
          values <- list.dirs(param, recursive = FALSE)
          ss <- SummaryStats()
          ss$loadValidatedResults(path = dataset$validated_results_file,
                                  datasetName = dName,
                                  bedFile = dataset$bed_file)
          for (value in values){
            outputFolder <- file.path(value, "output")
            try(ss$loadAlgorithmResults(outputFolder = outputFolder,
                                    algorithmName = paste0(alg,".", basename(param),".", basename(value)),
                                    datasetName = dName,
                                    bedFile = dataset$bed_file))
          }
          #write results in a CSV file----
          ss$writeCSVresults(file.path(param, paste0("results-", dName, "_", alg, "_", basename(param), ".csv")), dName)
        }
      }
    }
  }
}
