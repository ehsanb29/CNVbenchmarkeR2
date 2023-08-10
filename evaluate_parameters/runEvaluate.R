#Rscript runEvaluate.R [-t tools_file] [-d datasets_file] [-f include_temp_files]

print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
library(methods)
library("optparse")
library("dplyr")
options(scipen = 999)

option_list <- list(
  make_option(c("-t", "--tools"), type="character", default="tools.yaml",
              help="Path to tools file (yaml)", metavar="character"),
  make_option(c("-d", "--datasets"), type="character", default="datasets.yaml",
              help="Path to datasets file (yaml)", metavar="character"),
  make_option(c("-f", "--include_temp_files"), type="character", default="true",
              help="Include temporary files in the output folders", metavar="character")
);
opt_parser <- OptionParser(option_list=option_list);

# Load params
args <- parse_args(opt_parser);
tools <- yaml.load_file(args$tools)
datasets <- yaml.load_file(args$datasets)
tools.run <- list.files(paste0(getwd(),"/tools"), pattern="run", recursive=T, full.names = T)

for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
  for (i in seq_len(length(tools))){

    if(isTRUE(tools[[i]])){
      algName <- names(tools[i])

      filesParams <- list.files(paste0(getwd(), "/evaluate_parameters/", algName, "/", name), pattern= "params.yaml", recursive=TRUE, full.names=TRUE)


      #Run benchmark for every param

      for(j in seq_len(length(filesParams))){

      num.tool <- which(stringr::str_detect(tools.run, algName))
      # query<- paste(tools.run[num.tool],
      #       filesParams[j],
      #       file.path(paste("",algName, name, sep="/"), "datasets.yaml"),
      #       args$include_temp_files,
      #       " > ", file.path(paste0(filesParams[j] %>% dirname() %>% dirname(),"/logs"),paste0(algName,".log 2>&1")))
      query <- paste(tools.run[num.tool],
                   filesParams[j],
                     file.path(paste(getwd(), "evaluate_parameters",algName, name, sep="/"), "datasets.yaml"),
                     args$include_temp_files,
                     file.path(paste0(filesParams[j] %>% dirname() %>% dirname(),"/logs"),paste0(algName,".log")))

      cmd <- paste("sbatch evaluate_parameters/job.sh", query)
      print(cmd)
      system(cmd)

        }
    }

  }

  }

}

