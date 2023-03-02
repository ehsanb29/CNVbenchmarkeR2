# Runs Cobalt over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runCobalt.R [cobalt_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)
library(stringr)

# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  cobaltParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  cobaltParamsFile <- "tools/cobalt/cobaltParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
}

#Load the parameters file
params <- yaml.load_file(cobaltParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
print(paste("Params for this execution:", list(params)))

cobaltFolder <- file.path(params$cobaltFolder)


# go over datasets and run cobalt for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting cobalt for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)


    # Create output folder
    if (!is.null(params$outputFolder)) {
      outputFolder <- params$outputFolder
    } else {
      outputFolder <- file.path(getwd(), "output", paste0("cobalt-", name))
    }
    unlink(outputFolder, recursive = TRUE);
    dir.create(outputFolder, showWarnings = FALSE)


    # Read bamfiles and sample names
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    samplesName <- list.files(bamsDir, pattern = '*.bam$', full.names = FALSE) %>%
      stringr::str_replace_all( ".bam$", "")

    ## Create a coverage directory and coverage count for all samples----
    dir.create(file.path(outputFolder, "coverage"), showWarnings = FALSE)
    allCov <- file.path(outputFolder, "coverage/samples_coverages.bed")
    cmd <- paste( file.path(cobaltFolder, "covcounter"),
                  "--bed", bedFile,
                  "--bams",  paste(bamFiles, collapse = " "),
                  "--threads", params$threads,
                  ">", allCov
    )
    print(cmd); system(cmd)
    print("Covcounter finished");


    # Read all samples coverage file created in the last step
    allCovFile <- read.csv2(allCov, sep="\t")

    #Create a dataframe to store Cobalt results
    cnvData <- data.frame()

    # Go over each bam in bamFiles
    for (i in seq_len(length(bamFiles))){

    ## Training ----
      # Coverage for training model: Select params$training_samples (default: 50) random samples to compute the model, excluding the sample of interest from the training
      set.seed(123)
      randomSamples <- sample(c(1:length(bamFiles))[-i], params$training_samples, replace = FALSE)
      trainingSamples <- bamFiles[randomSamples]
      trainingSamplesName <-  samplesName[randomSamples]

      # Extract the coverages of the 50 samples randomly selected and store the results in a bed file
      covFile <- allCovFile %>%
        dplyr::select( "X.chrom" , "start", "end", contains(trainingSamplesName)) %>%
        dplyr::rename("#chrom" = X.chrom)
      trainingCov <- file.path(outputFolder, paste0("coverage/training_", samplesName[i], "_coverages.bed"))
      write.table(x = covFile,
                  file = trainingCov,
                  sep="\t",
                  quote=FALSE,
                  row.names = FALSE)


      # Create training directory
      dir.create(file.path(outputFolder, "training"), showWarnings = FALSE)

      # Build the model with the training samples (n=50)
      trainingModel <- file.path(outputFolder, paste0("training/training_model", samplesName[i], ".model"))
      cmd <- paste(file.path(cobaltFolder,"cobalt"), "train",
                   "--depths", trainingCov,
                   "-o", trainingModel,
                   "--chunk-size", params$chunk_size,
                   #"--no-mask", params$no_mask,
                   "--var-cutoff", params$var_cutoff,
                   "--min-depth", params$min_depth,
                   "--low-depth-trim-frac", params$low_depth_trim_frac,
                   "--high-depth-trim-frac", params$high_depth_trim_frac,
                   "--high-cv-trim-frac", params$high_cv_trim_frac,
                   "--cluster-width", params$cluster_width)
      print(cmd); system(cmd);


      ## Testing----

      # Create new test directories
      dir.create(file.path(outputFolder, "test"), showWarnings = FALSE)
      dir.create(file.path(outputFolder, paste0("test/coverages")), showWarnings = FALSE)
      dir.create(file.path(outputFolder, paste0("test/results")), showWarnings = FALSE)

      # Extract coverage of the sample of interest and store it in a bed file
      sampleCovFile <- allCovFile %>%
        dplyr::select( "X.chrom" , "start", "end", contains(samplesName[i])) %>%
        dplyr::rename("#chrom" = X.chrom)
      testCov <- file.path(outputFolder, paste0("test/coverages/", samplesName[i], "_coverages.bed"))
      write.table(x = sampleCovFile,
                  file = testCov,
                  sep="\t",
                  quote=FALSE,
                  row.names = FALSE)

      # Execute the model with a bed output format
      testResults <- file.path(outputFolder, paste0("test/results/test", samplesName[i], "_results.bed"))
      cmd <- paste(file.path(cobaltFolder, "cobalt"), "predict",
                   "-m", trainingModel,
                   "-d", testCov,
                   "-o", testResults)
      print(cmd); system(cmd);

      # Execute the model with a vcf file output
      #testResultsvcf <- file.path(outputFolder, paste0("test/results/test", samplesName[i], "_results.vcf"))
      # cmd <- paste(file.path(cobaltFolder, "cobalt"), "predict",
      #              "-m", trainingModel,
      #              "-d", testCov,
      #              "-o", testResultsvcf,
      #              "--vcf",
      #              "--reference", fastaFile)


      # Merge results of all samples in cnvData----
      testResultDf<- read.csv2(testResults, sep = "\t")
      if(nrow(testResultDf)>0){
        testResultDf <- testResultDf %>%
          dplyr::mutate(sample=samplesName[i])
        cnvData <- rbind(cnvData, testResultDf)
      }
    }

    #Rename cnvData columns
    names(cnvData) <- c("chr", "start", "end", "copy_number", "quality", "targets", "Sample")



    #Convert bed file and cnvData results to GRanges class object
    bedGR <- regioneR::toGRanges(bedFile)
    resultsGR <- regioneR::toGRanges(cnvData)

    # Create a dataframe to store results with the associated gene
    resGenes <- data.frame()

    # go over each CNV to identify genes
    for (i in seq_len(length(resultsGR))){

      #get the CNV
      cnvGR <- regioneR::toGRanges(resultsGR[i,])

      # Merge the cnv with the bed file to know the "affected" genes.
      # Mutate columns to include information related to the cnv of interest (quality, the copy_number, Sample..)
      #Group by gene to only have one row per gene and get the smallest start coordinate and the biggest end coordinate
      res <- as.data.frame(IRanges::subsetByOverlaps(bedGR, cnvGR)) %>%
        dplyr::mutate(copy_number = cnvGR$copy_number,
               quality = cnvGR$quality,
               #width = cnvGR$width,
               targets = cnvGR$targets,
               Sample = resultsGR$Sample[i]) %>%
        dplyr::group_by(name, Sample, seqnames) %>%
        dplyr::summarise(start_gene = min(start),
                         end_gene= max (end),
                         copy_number = min(copy_number),
                         quality = min(quality),
                         width = sum(width))

      # Merge each cnv result per gene with the rest
      resGenes <- rbind(resGenes, res)
    }

    # Convert copy_number nomenclature in CNV.type column. When copy_number <2 it is a deletion. Otherwise a duplication.
    resGenes <- resGenes %>%
      dplyr::rowwise() %>%
      dplyr::mutate(CNV.type = ifelse(copy_number < 2,
                                      "deletion",
                                      "duplication"))

    # Save results----
    # Print results in a tsv file
    finalSummaryFile <- file.path(outputFolder, "cnvs_summary.tsv")
    write.table(resGenes, finalSummaryFile, sep = "\t", quote = F, row.names = F)

    # Save results in a GenomicRanges object
    message("Saving CNV GenomicRanges results")
    saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "name",
                        sampleColumn = "Sample", chrColumn = "seqnames", startColumn = "start_gene",
                        endColumn = "end_gene", cnvTypeColumn = "CNV.type")

  }
}
