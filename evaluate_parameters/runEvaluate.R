#Rscript runEvaluate.R [-t tools_file] [-d datasets_file] [-f include_temp_files]

print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
library(methods)
library("optparse")
library("dplyr")
source("tools/germlineCNVcaller/germlineCNVcallerUtils.r") # Load GATK wrap functions
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

#loop through each dataset
for (name in names(datasets)) {
  
  dataset <- datasets[[name]]
  
  if (dataset$include){   
    
    # Loop through each tool
    for (i in seq_len(length(tools))){  
      
      if(isTRUE(tools[[i]])){
        
        algName <- names(tools[i])
        params <- yaml.load_file(paste0("tools/", algName, "/", algName, "Params.yaml"))

        # Run common part (precalc) for some tools  -----
        
        # Common part: Certain steps for Viscap, AtlasCNV, and GermlineCNVcaller have a "precalc" phase that do not depend on
        # on parameter values (.yaml) so it can be run only once for all parameter executions. This saves cpu time and disk space.
        
        ##Common part for Viscap and AtlasCNV:  Create Depth of Coverage files
        if(algName == "viscap" | algName == "atlasCNV"){
          depthCoverageFolder <- paste0("./evaluate_parameters/", algName, "/", name, "/DepthOfCoverage")

          if(!dir.exists(depthCoverageFolder)| dir.exists(depthCoverageFolder) &&  length(list.files(depthCoverageFolder)) < 7) {
            
            #get files and folders
            bamsDir <- file.path(dataset$bams_dir)
            bedFile <- file.path(dataset$bed_file)
            fastaFile <- file.path(dataset$fasta_file)
            bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)

            #create depthCoverageFolder
            dir.create(depthCoverageFolder, showWarnings = FALSE)
            bam <- basename(bamFiles) %>% tools::file_path_sans_ext()
            for (i in seq_len(length(bamFiles))){
              cmd <- paste("singularity exec ", 
                           " -B ", paste0(bamsDir), " ",
                           params$gatk,  "DepthOfCoverage",
                           "-R",  fastaFile,
                           "-I", bamFiles[i],
                           "-O",  ifelse(algName == "viscap",
                                         file.path(depthCoverageFolder,bam[i]),
                                         file.path(depthCoverageFolder,paste0(bam[i], ".DATA"))),
                           "--output-format", "TABLE",
                           "-L",  bedFile )
              paste(cmd);system(cmd)
            }
          }
        }

        ##Common part for GermlineCNVcaller: Interval_list, PreprocessIntervals, AnnotateIntervals, FilterIntervals, collectReadCounts and DetermineGermlineContigPloidy
        if(algName == "germlineCNVcaller"){
          
          #Create path to folders
          preCalcFolder <- paste0("./evaluate_parameters/", algName, "/", name, "/preCalculated")
          readCountsFolder <- paste0(preCalcFolder, "/ReadCounts")
          intervalListFolder <- paste0(preCalcFolder, "/IntervalList")
          contigPloidyFolder <- paste0(preCalcFolder, "/contigPloidy")
          contigPriorsFile <- file.path(getwd(), "tools/germlineCNVcaller/contig-ploidy-prior.tsv")

          # get dataset paths
          bamsDir <- file.path(dataset$bams_dir)
          bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
          bedFile <- file.path(dataset$bed_file)
          fastaFile <- file.path(dataset$fasta_file)
          fastaDict <- file.path(dataset$fasta_dict)

          
          # launch precalc part if contig DetermineContigPloidy did not finished
          if(!dir.exists(readCountsFolder) | 
             (dir.exists(readCountsFolder) && length(list.files(contigPloidyFolder, "*contig_ploidy.tsv", recursive = T)) != length(bamFiles)) ){
            
            #create directories
            dir.create(preCalcFolder, showWarnings = FALSE)
            dir.create(readCountsFolder, showWarnings = FALSE)
            dir.create(intervalListFolder, showWarnings = FALSE)
            outputFolder <- file.path(getwd(), preCalcFolder)

            # Get gatk and contig file for later use
            gatk <- params$gatk

            # create intervals file (.intervals_list) from .bed file
            interval_list <- paste0(intervalListFolder,"/list.interval_list")
            createIntervalList(params$picard, bedFile, interval_list, fastaDict)
            
            # preprocess intervals
            preprocessInterval_list <- paste0(intervalListFolder,"/preprocessed_intervals.interval_list")
            preprocessIntervals(gatk, fastaFile, outputFolder, interval_list, preprocessInterval_list)
            
            # annotate intervals
            annotatedIntervals_list <- paste0(intervalListFolder,"/annotated_intervals.interval_list")
            annotateIntervals(gatk, fastaFile, outputFolder, preprocessInterval_list, annotatedIntervals_list)

            # CollectReadCounts
            collectReadCounts(gatk, bamFiles, outputFolder, preprocessInterval_list, readCountsFolder)
            
            # merge all hdf5 files for later use
            hdf5s <- list.files(readCountsFolder, pattern = '*.hdf5$', full.names = TRUE)
            all_hdf5_names <- paste(hdf5s, collapse=" -I ")
            
            # FilterIntervals
            filteredIntervals_list <- paste0(intervalListFolder,"/filtered_intervals.interval_list")
            filterIntervals(gatk, outputFolder, preprocessInterval_list, all_hdf5_names, annotatedIntervals_list, filteredIntervals_list)

            # Determine contig ploidy
            determineContigPloidy(gatk, contigPloidyFolder, bamsDir, outputFolder, contigPriorsFile, all_hdf5_names, filteredIntervals_list)
            
          }
          
        } 
        # End of common part

        # Launch tool -----

        #Run evaluate: list all the files of the parameters that have to be evaluated
        filesParams <- list.files(paste0(getwd(), "/evaluate_parameters/", algName, "/", name), pattern= "params.yaml", recursive=TRUE, full.names=TRUE)


        #Run benchmark for every parameter file
        for(j in seq_len(length(filesParams))){
          
          #create the parameters that should be given to job.sh file
          num.tool <- which(stringr::str_detect(tools.run, algName))
          query <- paste(tools.run[num.tool],
                       filesParams[j],
                         file.path(paste(getwd(), "evaluate_parameters",algName, name, sep="/"), "datasets.yaml"),
                         args$include_temp_files,
                         "true",
                         file.path(paste0(filesParams[j] %>% dirname() %>% dirname(),"/logs"),paste0(algName,".log")))

          # Get job name
          parts <- strsplit(dirname(filesParams[j]), "/")[[1]]
          jobName <- paste0(parts[c(length(parts) - 2, length(parts) - 1)], collapse = "_")

          # Launch
          cmd <- paste("qsub -N", jobName, "evaluate_parameters/job.sh", query)
          print(cmd, quote=FALSE); system(cmd)
          
          #if(j == 1){
          #  break;
          #}
          
        }
      }
      
    }
  }
}

