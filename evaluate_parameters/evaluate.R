#USAGE: Rscript optimizer.r [-t tools.yaml]  [-d dataset.yaml] (call from evaluate folder)
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
library(methods)
library("optparse")
options(scipen = 999)  # to disable scientific number notation
source("../utils/cnvStats.r")  # Load utils functions
source("../utils/optimizerUtils.r")  # Load opt. utils functions

# Return available options for the different params
getFixedValues = function(param, name){
  vl <- ValuesList()
  vals <- list()
  for (v in param$options){
    newVal <- Value()
    newVal$new(v)
    vals[[length(vals) + 1]] <- newVal
  }
  vl$new(vals, name, param$default, FALSE)
  
  return(vl)
}


# Build options list
option_list <- list(
  make_option(c("-t", "--tools"), type="character", default="tools.yaml",
              help="Path to tools file (yaml)", metavar="character"),
  make_option(c("-d", "--datasets"), type="character", default="datasets.yaml",
              help="Path to datasets file (yaml)", metavar="character")
);
opt_parser <- OptionParser(option_list=option_list);

# Load params
args <- parse_args(opt_parser);
tools <- yaml.load_file(args$tools)
datasets <- yaml.load_file(args$datasets)


#CREATE NEW PARAMS FILE

#For each dataset


for (name in names(datasets)) {
  dataset <- datasets[[name]]

for (i in 1:seq_len(length(tools))){
  
  if(isTRUE(tools[[i]])){
    algName <- names(tools[i])
    
    #create a folder for the dataset
    folder_name <- file.path(algName, name)
    dir.create(folder_name)
    
    
    #default params 
    params_default <- yaml.load_file(paste0("../tools/", algName, "/", algName, "Params.yaml"))
    params_optimizer <- yaml.load_file(paste0(algName, "/", algName, "Params.yaml"))
  
    paramsValues <- list()  # list of values to be tested for each param
    for (i in 1:length(params_optimizer)){
      p <-params_optimizer[i]
      paramName <- names(p)
      #create a folder for the param
      folder_param <- file.path(folder_name, paramName)
      dir.create(folder_param)
      #get params options
      options_param <- p[[paramName]]$options
      #for every_option
      for(j in seq_len(length(options_param))){
        value_param <- options_param[j]%>% unlist()
        print(value_param)
        folder_param_value <- file.path(folder_param, value_param )
        dir.create(folder_param_value , showWarnings = FALSE)
        
        # create logs and input folder to be used by jobs
        logsFolder <- file.path(folder_param_value, "logs")
        unlink(logsFolder, recursive = TRUE);
        dir.create(logsFolder, showWarnings = FALSE)
        inputFolder <- file.path(folder_param_value, "input")
        unlink(inputFolder, recursive = TRUE);
        dir.create(inputFolder, showWarnings = FALSE)
        outputFolder <- file.path(folder_param_value, "output")
        unlink(outputFolder, recursive = TRUE);
        dir.create(outputFolder, showWarnings = FALSE)
        
        #Create the new yaml
        params_tunned <- params_default
        params_tunned[[paramName]] <- value_param
        write_yaml(params_tunned, file = file.path(inputFolder, paste0(algName, "params.yaml")))
        
        
      }
    
    
    }
  
  }
  
}
  
}



