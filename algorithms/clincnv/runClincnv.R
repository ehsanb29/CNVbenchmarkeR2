# Runs ClinCNV on the datasets configured at [datasets_params_file]
#USAGE: Rscript runClincnv.r [clincnv_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  clincnvParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  clincnvParamsFile <- "clincnvParams.yaml"
  datasetsParamsFile <- "../../datasets.yaml"
}


# translates DEL/DUP into a common format
auxCNname <- function(x) {
  if (x %in% c("CN0", "CN1")) return("deletion") 
  else if (x %in% c("CN3", "CN4")) return("duplication")
}

#Load the parameters file  
params <- yaml.load_file(clincnvParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
print(paste("Params for this execution:", list(params)))

# get clincnv folder
clincnvFolder <- file.path(params$clincnvFolder)


# run clincnv on selected datasets
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting ClinCNV for", name, "dataset", sep=" "))
    
    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)
    
    # Create output folder
    if (!is.null(params$outputFolder)) {
      outputFolder <- params$outputFolder
    } else{
      outputFolder <- file.path(getwd(), "output", paste0("clincnv-", name))  
    }
    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder)
    
    # create folders to be used later
    ontargetFolder <- file.path(outputFolder, "ontargetCov")
    dir.create(ontargetFolder)
    
    # 1: Calculate on-target coverage for each sample
    bamFiles <- list.files(dataset$bams_dir, "*.bam$")
    for (f in bamFiles){
      fullBamFile <- file.path(dataset$bams_dir, f)
      sampleName <- tools::file_path_sans_ext(f)
      fullOntargetFile <- file.path(ontargetFolder, paste0(sampleName, ".cov"))
      cmd <- paste(file.path(params$ngsbitsFolder, "BedCoverage"), '-bam', fullBamFile,
                   '-in', dataset$annotated_bed_file, '-decimals 4 -out ', fullOntargetFile)
      system(cmd)
    }
    
    # 2: Merge coverage files into one table
    covFiles <- list.files(ontargetFolder, "*.cov$", full.names = T)
    covData <- read.table(covFiles[1], stringsAsFactors = F)[,c(1:3)]
    names(covData) <- c("chr", "start", "end")
    coverageFile <- file.path(outputFolder, "coverage.cov")
    for (f in covFiles){
      sampleName <- tools::file_path_sans_ext(basename(f))
      covData[sampleName] <-  read.table(f, stringsAsFactors = F)[6]
    }
    write.table(covData, coverageFile, sep = "\t", row.names = F, quote = F)
    
    
    # 3: Call germline CNVs
    # --lengthG has to be 0 since we expect to find CNVs with at least one region
    # --maxNumGermCNVs: 20  # maximum number of allowed germline CNVs per sample. Default value is 10000, but we chose smaller value for gene panel following author's advice
    cmd <- paste("Rscript", file.path(params$clincnvFolder, "clinCNV.R"), 
                 "--normal", coverageFile, "--out", outputFolder, 
                 "--bed", dataset$annotated_bed_file, "--scoreG", params$scoreG,
                 "--minimumNumOfElemsInCluster", params$minimumNumOfElemsInCluster,
                 "--lengthG 0 --numberOfThreads 2 --reanalyseCohort")
    system(cmd)
    
    
    # 4: Summarize results into one file
    cnvFiles <- list.files(file.path(outputFolder, "normal"), "*_cnvs.tsv$", full.names = T, recursive = T)
    cnvData <- data.frame()
    for ( f in cnvFiles){
      nLines <- length(readLines(f))
      if (nLines > 9){
        sampleName <- tools::file_path_sans_ext(basename(f))
        sampleName <- strsplit(sampleName, "_cnvs")[[1]][1]
        fData <- read.table(f, stringsAsFactors = F, skip = 9)
        fData$Sample <- sampleName
        cnvData <- rbind(cnvData, fData)  
      }
    }
    names(cnvData) <- c("chr","start","end","CN_change","loglikelihood","no_of_regions","length_KB",
                        "potential_AF","genes","qvalue", "Sample")
    cnvData$CNV.type <- ifelse(cnvData$CN_change < 2, "deletion", "duplication") # assign common values
    finalSummaryFile <- file.path(outputFolder, "cnvs_summary.tsv")
    write.table(cnvData, finalSummaryFile, sep = "\t", quote = F, row.names = F)
    
    # Save results in GRanges format
    message("Saving CNV GenomicRanges results")
    saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "genes", 
                        sampleColumn = "Sample", chrColumn = "chr", startColumn = "start",
                        endColumn = "end", cnvTypeColumn = "CNV.type")
    
    print(paste("ClinCNV for", name, "dataset finished", sep=" "))
    cat("\n\n\n")
  }
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)

