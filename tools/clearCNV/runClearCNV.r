# Runs clearCNV over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runclearCNV.R [clearCNV_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)


# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  clearCNVParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  clearCNVParamsFile <- "tools/clearCNV/clearCNVParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
}

#Load the parameters file
params <- yaml.load_file(clearCNVParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
print(paste("Params for this execution:", list(params)))

#get conda folder and the environment where clearCNV is insalled
condaFolder <-file.path(params$condaFolder)
environmentName <- file.path(params$environmentName)

#get clearCNV folder
clearCNVFolder <- file.path(params$clearCNVFolder)


# go over datasets and run clearCNV for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting clearCNV for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)


    # Create output folder
    if (!is.null(params$outputFolder)) {
      outputFolder <- params$outputFolder
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("clearCNV-", name))
    }
    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder, showWarnings = FALSE)

    # Get bam files
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)

    ## Create bam directory and merge all bams in a txt file
    dir.create(file.path(outputFolder, "bams"))
    bamTxt <- file.path(outputFolder, paste0("bams/all_bams.txt"))
    write.table(x = paste(bamFiles, sep="\n"), file = bamTxt , quote = FALSE, row.names = FALSE, col.names=FALSE)


    # Call germline CNV caller
    cmd <- paste( condaFolder,
                "run -n", environmentName,
                clearCNVFolder, "workflow_cnv_calling",
                "-w", outputFolder,
                "-p", name,
                "-r", fastaFile,
                "-b", bamTxt,
                "-d", bedFile,
                "-k", file.path("tools/clearCNV", "consensusBlacklist.bed"),
                "-c", params$cores,
                "--expected_artefacts", params$expected_artefacts,
                "--sample_score_factor", params$sample_score_factor,
                "--minimum_group_sizes", params$minimum_group_sizes,
                "--zscale", params$zscale,
                "--size", params$size,
                "--del_cutoff", params$del_cutoff,
                "--dup_cutoff", params$dup_cutoff,
                "--trans_prob", params$trans_prob
                #"--cluster_configfile", params$cluster_configfile,
                #"--drmaa_mem", params$drmaa_mem,
                #"--drmaa_time", params$drmaa_time,
                #"--drmaa_jobs", params$drmaa_jobs
                )
   paste(cmd);system(cmd);

   # Read output file, add CNV.type column and write the results in a tsv file
   resFile <- file.path(outputFolder, paste0(name, "/results/cnv_calls.tsv"))
   resDF <- read.table(resFile, header = TRUE) %>%
     dplyr::mutate(CNV.type = ifelse(aberration == "DEL",
                                     "deletion",
                                     "duplication"))
   write.table(resDF, file.path(outputFolder, "cnv_calls.tsv"), sep="\t", quote=F, row.names = FALSE, col.names = TRUE)

   # Save results----
   # Path to  tsv file
   finalSummaryFile <- file.path(outputFolder, "cnv_calls.tsv")
   # Save results in a GenomicRanges object
   message("Saving CNV GenomicRanges results")
   saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "gene",
                       sampleColumn = "sample", chrColumn = "chr", startColumn = "start",
                       endColumn = "end", cnvTypeColumn = "CNV.type")

  }
}
