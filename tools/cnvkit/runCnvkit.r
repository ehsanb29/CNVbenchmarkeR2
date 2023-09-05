# Runs cnvkit over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runCnvkit.R [clearCNV_params_file] [datasets_params_file] [include_temp_files]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)


# Read args
args <- commandArgs(TRUE)
if(length(args)>0) {
  cnvkitParamsFile <- args[1]
  datasetsParamsFile <- args[2]
  includeTempFiles <- args[3]
} else {
  cnvkitParamsFile <- "tools/cnvkit/cnvkitParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
  includeTempFiles <- "true"
}

#Load the parameters file
params <- yaml.load_file(cnvkitParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
print(paste("Params for this execution:", list(params)))

#get conda folder and the environment where cnvkit is insalled
condaFolder <-file.path(params$condaFolder)
environmentName <- file.path(params$environmentName)

#get cnvkit folder
cnvkitFolder <- file.path(params$cnvkitFolder)

# go over datasets and run cnvkit for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting cnvkit for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)
    currentFolder <- getwd()

    # Create output folder
    if (!is.null(params$outputFolder)) {
      if(stringr::str_detect(params$outputFolder, "^./")) params$outputFolder <- stringr::str_sub(params$outputFolder, 3, stringr::str_length(params$outputFolder))
      outputFolder <- file.path(currentFolder, params$outputFolder)
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("cnvkit-", name))
    }
    #unlink(outputFolder, recursive = TRUE);
    #dir.create(outputFolder, showWarnings = FALSE)

    # Get bam files
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)

    # Access: calculate accessible coordinates of the chromosome sequences
    # cmd <- paste(condaFolder,
    #              "run -n", environmentName,
    #              cnvkitFolder, "access",
    #              fastaFile,
    #              "-s", params$minGapSizeAccess,
    #              "-o", file.path(outputFolder, "access.bed"))

    cmd <- paste(".",
                 cnvkitFolder, "\n",
                 #"srun",
                 "cnvkit.py",  "access",
                 fastaFile,
                 "-s", params$minGapSizeAccess,
                 "-o", file.path(outputFolder, "access.bed"))

    print(cmd);system(cmd)

    #Autobin: Quickly estimate read counts or depths in a BAM file to estimate reasonable on- and (if relevant) off-target bin sizes.
    # cmd <- paste(condaFolder,
    #              "run -n", environmentName,
    #              cnvkitFolder, "autobin",
    #              paste0(bamsDir, "/*.bam"),
    #              "-m hybrid",
    #              "-t", bedFile,
    #              "-g", file.path(outputFolder, "access.bed"),
    #              "-b", params$bpPerBinAutobin,
    #              "--target-max-size", params$targetMaxSizeAutobin,
    #              "--target-min-size", params$targetMinSizeAutobin,
    #              "--antitarget-max-size", params$antitargetMaxSizeAutobin,
    #              "--antitarget-min-size", params$antitargetMinSizeAutobin,
    #              "--target-output-bed",  file.path(outputFolder, "ROIs-ICR96-panelcnDataset.target.bed"),
    #              "--antitarget-output-bed",  file.path(outputFolder, "ROIs-ICR96-panelcnDataset.antitarget.bed")
    # )

    cmd <- paste(".",
                 cnvkitFolder, "\n",
                 #"srun",
                 "cnvkit.py", "autobin",
                 paste0(bamsDir, "/*.bam"),
                 "-m hybrid",
                 "-t", bedFile,
                 "-g", file.path(outputFolder, "access.bed"),
                 "-b", params$bpPerBinAutobin,
                 "--target-max-size", params$targetMaxSizeAutobin,
                 "--target-min-size", params$targetMinSizeAutobin,
                 "--antitarget-max-size", params$antitargetMaxSizeAutobin,
                 "--antitarget-min-size", params$antitargetMinSizeAutobin,
                 "--target-output-bed",  file.path(outputFolder, "ROIs-ICR96-panelcnDataset.target.bed"),
                 "--antitarget-output-bed",  file.path(outputFolder, "ROIs-ICR96-panelcnDataset.antitarget.bed")
    )
    print(cmd);system(cmd)
    dir.create(file.path(outputFolder, "calls"))

    for(i in seq_len(length(bamFiles))){
      #For each sample
      sampleName <- basename(tools::file_path_sans_ext(bamFiles))

      #Get the coverage for target and antitarget
      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "coverage",
      #              bamFiles[i],
      #              "-q", params$minMapqCoverage,
      #              file.path(outputFolder, "ROIs-ICR96-panelcnDataset.target.bed"),
      #              "-o", file.path(outputFolder, paste0(sampleName[i], ".targetcoverage.cnn" )))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py",  "coverage",
                   bamFiles[i],
                   "-q", params$minMapqCoverage,
                   file.path(outputFolder, "ROIs-ICR96-panelcnDataset.target.bed"),
                   "-o", file.path(outputFolder, paste0(sampleName[i], ".targetcoverage.cnn" )))
      print(cmd);system(cmd)

      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "coverage",
      #              bamFiles[i],
      #              "-q", params$minMapqCoverage,
      #              file.path(outputFolder, "ROIs-ICR96-panelcnDataset.antitarget.bed"),
      #              "-o", file.path(outputFolder, paste0(sampleName[i], ".antitargetcoverage.cnn" )))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py","coverage",
                   bamFiles[i],
                   "-q", params$minMapqCoverage,
                   file.path(outputFolder, "ROIs-ICR96-panelcnDataset.antitarget.bed"),
                   "-o", file.path(outputFolder, paste0(sampleName[i], ".antitargetcoverage.cnn" )))
      #print(cmd); system(cmd)
    }

    for(i in seq_len(length(bamFiles))){
      #Build a reference for each sample, including all samples except the sample of interest
      #Get cnn Target and antitarget files
      cnnTarget <- paste0(file.path(outputFolder, paste0(sampleName[-i], ".targetcoverage.cnn" )), collapse = " ")
      cnnAntitarget <- paste0(file.path(outputFolder, paste0(sampleName[-i], ".antitargetcoverage.cnn" )), collapse=" ")
      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "reference",
      #              cnnTarget, cnnAntitarget,
      #              "-f", fastaFile,
      #              "--min-cluster-size", params$minClusterSizeReference,
      #              "-o", file.path(outputFolder, paste0(sampleName[i], "reference.cnn" )))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py", "reference",
                   cnnTarget, cnnAntitarget,
                   "-f", fastaFile,
                   "--min-cluster-size", params$minClusterSizeReference,
                   "-o", file.path(outputFolder, paste0(sampleName[i], "reference.cnn" )))
      paste(cmd);system(cmd)



      #Fix: uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference.
      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "fix",
      #              file.path(outputFolder, paste0(sampleName[i], ".targetcoverage.cnn" )),
      #              file.path(outputFolder, paste0(sampleName[i], ".antitargetcoverage.cnn" )),
      #              file.path(outputFolder, paste0(sampleName[i], "reference.cnn" )),
      #              "-o",  file.path(outputFolder, paste0(sampleName[i], ".cnr")))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py", "fix",
                   file.path(outputFolder, paste0(sampleName[i], ".targetcoverage.cnn" )),
                   file.path(outputFolder, paste0(sampleName[i], ".antitargetcoverage.cnn" )),
                   file.path(outputFolder, paste0(sampleName[i], "reference.cnn" )),
                   "-o",  file.path(outputFolder, paste0(sampleName[i], ".cnr")))
     # paste(cmd); system(cmd)


      #Segment infer discrete copy number segments from the given coverage tamble:
      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "segment",
      #              file.path(outputFolder, paste0(sampleName[i], ".cnr")),
      #              "-m", params$methodSegment,
      #              #"-t", params$thresehold,
      #              "--drop-outliers", params$DropOutliersSegment,
      #              "-o", file.path(outputFolder, paste0(sampleName[i], ".cns")))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py", "segment",
                   file.path(outputFolder, paste0(sampleName[i], ".cnr")),
                   "-m", params$methodSegment,
                   #"-t", params$thresehold,
                   "--drop-outliers", params$DropOutliersSegment,
                   "-o", file.path(outputFolder, paste0(sampleName[i], ".cns")))
      #paste(cmd);system(cmd)

      #Segmetrics: Calculate summary statistics of the residual bin-level log2 ratio estimates from the segment means, similar to the existing metrics command, but for each segment individually.
      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "segmetrics",
      #              file.path(outputFolder, paste0(sampleName[i], ".cnr")),
      #              "-s", file.path(outputFolder, paste0(sampleName[i], ".cns")),
      #              "--ci",
      #              "--alpha", params$alphaSegmetrics,
      #              "-b", params$bootstrapSegmetrics,
      #              "-o", file.path(outputFolder, paste0(sampleName[i],"segment.cns")))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py", "segmetrics",
                   file.path(outputFolder, paste0(sampleName[i], ".cnr")),
                   "-s", file.path(outputFolder, paste0(sampleName[i], ".cns")),
                   "--ci",
                   "--alpha", params$alphaSegmetrics,
                   "-b", params$bootstrapSegmetrics,
                   "-o", file.path(outputFolder, paste0(sampleName[i],"segment.cns")))

      #paste(cmd);system(cmd)

      #call cnvs
      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "call",
      #              file.path(outputFolder, paste0(sampleName[i], "segment.cns")),
      #              "--filter", "ci",
      #              "-m", params$methodCall,
      #              #"-t", params$thresholdsCall,
      #              "-o", file.path(outputFolder, paste0("calls/", sampleName[i], ".call.cns")))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py", "call",
                   file.path(outputFolder, paste0(sampleName[i], "segment.cns")),
                   "--filter", params$filterCall,
                   "-m", params$methodCall,
                   paste0("-t=", params$thresholdsCall_del, params$thresholdsCall_loss, params$thresholdsCall_gain, params$thresholdsCall_amp),
                   "-o", file.path(outputFolder, paste0("calls/", sampleName[i], ".call.cns")))
      #paste(cmd);system(cmd)

      #bintest

      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "bintest",
      #              file.path(outputFolder, paste0(sampleName[i], ".cnr")),
      #              "-s", file.path(outputFolder, paste0("calls/", sampleName[i], ".call.cns")),
      #              "-t",
      #              "-a", params$alphaBintest,
      #              "-o", file.path(outputFolder, paste0(sampleName[i], "bintest.cns"))
      #              )
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py", "bintest",
                   file.path(outputFolder, paste0(sampleName[i], ".cnr")),
                   "-s", file.path(outputFolder, paste0("calls/", sampleName[i], ".call.cns")),
                   "-t",
                   "-a", params$alphaBintest,
                   "-o", file.path(outputFolder, paste0(sampleName[i], "bintest.cns"))
      )
      #print(cmd); system(cmd)
      # cmd <- paste(condaFolder,
      #              "run -n", environmentName,
      #              cnvkitFolder, "call",
      #              file.path(outputFolder, paste0(sampleName[i], "bintest.cns")),
      #              "-m", params$methodCall,
      #              #"-t", params$thresholdsCall,
      #              "-o", file.path(outputFolder, paste0("calls/", sampleName[i], "bintest.call.cns")))
      cmd <- paste(".",
                   cnvkitFolder, "\n",
                   #"srun",
                   "cnvkit.py", "call",
                   file.path(outputFolder, paste0(sampleName[i], "bintest.cns")),
                   "-m", params$methodCall,
                   #"-t", params$thresholdsCall,
                   "-o", file.path(outputFolder, paste0("calls/", sampleName[i], "bintest.call.cns")))

      #print(cmd); system(cmd)


    }



#Put all results in the same file
cnvData  <- data.frame()

for(s in sampleName){
  #read the call files and merge them
  callFiles <- read.table(file.path(outputFolder, paste0("calls/", s, ".call.cns")), header=TRUE) %>% dplyr::mutate(p_bintest = NA) %>% dplyr::relocate(p_bintest, .after = depth)
  callBinFiles <- read.table(file.path(outputFolder, paste0("calls/", s, "bintest.call.cns")), header=TRUE)
  sam <- rbind(callFiles, callBinFiles)

  #Put sample names and the consequence of the cnv
  sam2 <- sam %>%
    dplyr::rowwise() %>%
    dplyr::filter(cn != 2) %>%
    dplyr::mutate(Sample= s,
                  CNV.type = ifelse(cn < 2 ,
                                    "deletion",
                                    ifelse(cn>2,
                                           "duplication",
                                            "")))
  #Merge all in the same file
  cnvData <- rbind(cnvData, sam2)
}

#Read the bed file and pass cnv results to a GRanges object
bedGR <- regioneR::toGRanges(bedFile)
resultsGR <- regioneR::toGRanges(cnvData %>% as.data.frame())

# Create a dataframe to store results with the associated gene
resGenes <- data.frame()
for (i in seq_len(length(resultsGR))){
  #get the CNV
  cnvGR <- regioneR::toGRanges(resultsGR[i,])
  # Merge the cnv with the bed file to know the "affected" genes.
  # Mutate columns to include information related to the cnv of interest (quality, the copy_number, Sample..)
  #Group by gene to only have one row per gene and get the smallest start coordinate and the biggest end coordinate
  res <- as.data.frame(IRanges::subsetByOverlaps(bedGR, cnvGR)) %>%
    dplyr::mutate(CNV.type = cnvGR$CNV.type,
                  log2 = cnvGR$log2,
                  #width = cnvGR$width,
                  cn = cnvGR$cn,
                  Sample = resultsGR$Sample[i]) %>%
    dplyr::group_by(name, Sample, seqnames, CNV.type) %>%
    dplyr::summarise(start_gene = min(start),
                     end_gene= max (end),
                     cn = min(cn),
                     log2 = min(log2))
  # Merge each cnv result per gene with the rest
  resGenes <- rbind(resGenes, res)
}

# Save results----
# Print results in a tsv file
finalSummaryFile <- file.path(outputFolder, "cnvs_summary.tsv")
write.table(resGenes, finalSummaryFile, sep = "\t", quote = F, row.names = F)

# Save results in a GenomicRanges object
message("Saving CNV GenomicRanges results")
saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "name",
                    sampleColumn = "Sample", chrColumn = "seqnames", startColumn = "start_gene",
                    endColumn = "end_gene", cnvTypeColumn = "CNV.type")

#Delete temporary files if specified
if(includeTempFiles == "false"){
  filesAll <- list.files(outputFolder, full.names = TRUE)
  filesToKeep <- c("failedRois.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
  filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse= "|"), filesAll, value=TRUE))])
  do.call(unlink, filesToRemove)
}

  }
}
