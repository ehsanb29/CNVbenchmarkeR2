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
                 "--numberOfThreads", params$threads, "--lengthG 0 --reanalyseCohort")
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
    
    # assign common values for CNV type
    cnvData$CNV.type <- ifelse(cnvData$CN_change < 2, "deletion", "duplication") 
    
    # Split multi-gene CNVs in multiple lines
    bedGR <- regioneR::toGRanges(dataset$bed_file)
    cnvDataMultiGene <- setNames(data.frame(matrix(ncol = ncol(cnvData), nrow = 0)), names(cnvData))
    for (i in seq_len(nrow(cnvData))){
      
      # get each CNV and check which ROIs overlaps
      cnv <- cnvData[i,]
      cnvGR <- regioneR::toGRanges(cnv)
      GenomeInfoDb::seqlevelsStyle(cnvGR) <- "Ensembl"  # removes "chr" for chromosomes
      res <- as.data.frame(IRanges::subsetByOverlaps(bedGR, cnvGR))
      
      # add a line for each gene
      for (gene in unique(res$name)) {
        cnvToAdd <- cnvData[i,]
        cnvToAdd$genes <- gene
        cnvToAdd$start <- res[res$name == gene, "start"][1] - 1
        cnvToAdd$end <- tail(res[res$name == gene, "end"], n = 1)
        cnvDataMultiGene <- rbind(cnvDataMultiGene, cnvToAdd)
      }
    }
    
    # Write final CNV results file
    finalSummaryFile <- file.path(outputFolder, "cnvs_summary.tsv")
    write.table(cnvDataMultiGene, finalSummaryFile, sep = "\t", quote = F, row.names = F)
    
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

