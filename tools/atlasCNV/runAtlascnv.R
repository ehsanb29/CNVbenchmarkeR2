# Runs Atlas-CNV over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runAtlasCNV.R [atlasCNVparams_file] [datasets_params_file] [include_temp_files]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

# Read args ----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  paramsFile <- args[1]
  datasetsParamsFile <- args[2]
  includeTempFiles <- args[3]
} else {
  paramsFile <- "params.yaml"
  datasetsParamsFile <- "../../datasets.yaml"
  includeTempFiles <- "true"
}

#Load the parameters file  ----
params <- yaml.load_file(paramsFile)
datasets <- yaml.load_file(datasetsParamsFile)


# print params
print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))


# Get AtlasCNV folder ----
atlascnvFolder <- file.path(params$atlascnvFolder)
#locate the reference and contig files
fastaFile <- file.path(datasets$ICR96$fasta_file)
contigFile <- file.path(params$contigFile)
gatkFolder <- file.path(params$gatkFolder)
picardJar <- file.path(params$picardJar)
print(picardJar)


#create input files required for running GATK
# Dictionary ----
#create dictionary file (.dict) from reference genome (.)
# if(!file.exists(paste0(atlascnvFolder, "/reference.dict"))){
#   cmd_dictionary<- paste0(" java -jar $PICARD",
#                           " CreateSequenceDictionary",
#                           " -R ", fastaFile,
#                           " -O ", atlascnvFolder, "/reference.dict")
#
#
#   paste(cmd_dictionary);system(cmd_dictionary);
# }
#create the variable for the reference dictionary
#dictionary<-paste0(atlascnvFolder,"/reference.dict")


fastaDict <- paste0(tools::file_path_sans_ext(fastaFile), ".dict")

if(!file.exists(fastaDict)){
  cmd <- paste0(" java -jar ", picardJar,
                " CreateSequenceDictionary",
                " -R ", fastaFile,
                " -O ", fastaDict)


  paste(cmd);system(cmd);
}


#create the folder for the coverage files

# coverageFiles <- file.path(paste0(atlascnvFolder, "/coverage_files"))
# if (!file.exists(coverageFiles)){
#   dir.create(coverageFiles)
# } else {
#   unlink(coverageFiles, recursive = TRUE)
#   dir.create(coverageFiles, recursive=TRUE)
# }


# Dataset iteration ----
# go over datasets and run cnvkit for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting Atlas-CNV for", name, "dataset", sep=" "))


    # Create output folder
    outputFolder <- file.path(getwd(), "output", paste0("atlasCNV-", name))
    print(outputFolder)
    dir.create(outputFolder, showWarnings = FALSE)

    # extract fields

    bamsDir <- file.path(dataset$bams_dir)
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    bedFile <- file.path(dataset$bed_file)
    #annotatedBedFile<-file.path(dataset$annotated_bed_file)

    # Input files for Atlas-CNV ----
    ## Interval list----
    #create intervals file (.intervals_list) from .bed file
    #if(!file.exists(paste0(atlascnvFolder, "/list.interval_list"))){
    # cmd_interval <- paste0( " java -jar  $PICARD",
    #                         " BedToIntervalList",
    #                         " -I ", bedFile,
    #                         " -O ", atlascnvFolder,"/list.interval_list",
    #                         " -SD ", dictionary)

    #paste(cmd_interval);system(cmd_interval); }


    #if(!file.exists(paste0(atlascnvFolder, "/list.interval_list"))){
    #dir.create(file.path(outputFolder,"list.interval_list"), showWarnings = FALSE)
    # print("file")
    cmd_interval <- paste0( " java -jar ", picardJar,
                            " BedToIntervalList",
                            " -I ", bedFile,
                            " -O ", outputFolder,"/list.interval_list",
                            " -SD ", fastaDict)
    paste(cmd_interval);system(cmd_interval);


    print("FIns aqui funciona")

    ## Panel file ----
    #create panel file to input in Atlas
    ##
    tempBed<-read.csv(bedFile, sep="\t", header = FALSE)

    #tempBed<-read.csv(bedFile, sep="\t", header = FALSE)

    panel<-tempBed %>%
      arrange(V1) %>%
      mutate(Exon_Target = paste0(V1,":", V2%>% as.numeric() +1,"-", V3))  %>%
      mutate(Gene_Exon=V4) %>%
      mutate(Call_CNV = ifelse(grepl("^Y", Exon_Target), "N", "Y"))%>%
      mutate(RefSeq=V4) %>%
      select(-starts_with("V"))
    #export panel file
    #write.table(panel, file=paste0(atlascnvFolder,'/panel',name,'.txt'),sep="\t", col.names=T, row.names=FALSE,quote=FALSE)
    write.table(panel, file=paste0(outputFolder,'/panel',name,'.txt'),sep="\t", col.names=T, row.names=FALSE,quote=FALSE)

    ## Sample file ----

    tempSample<-as.data.frame(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bamFiles)))
    sample<-tempSample %>%
      mutate(sex="F") %>%
      mutate(mp=name)
    #export sample file
    #write.table(sample, file=paste0(atlascnvFolder,'/',name,'.sample'),sep="\t", col.names=F, row.names=FALSE,quote=FALSE)
    write.table(sample, file=paste0(outputFolder,'/',name,'.sample'),sep="\t", col.names=F, row.names=FALSE,quote=FALSE)
    rm(tempSample, tempBed)

    #create folder for coverage files
    #coverageFiles <- file.path(paste0(outputFolder, "/coverage_files"))
    #create Depth of coverageFolder
    depthCoverageFolder <- paste0(outputFolder,"/DepthOfCoverage")
    dir.create(depthCoverageFolder, showWarnings = FALSE)

    #get depth of coverage for each BAM file, using GATK
    for (bam in bamFiles){
      bam <- basename(bamFiles) %>% tools::file_path_sans_ext()
      for (i in seq_len(length(bamFiles))){
        cmd <- paste(file.path(gatkFolder, "gatk"),  "DepthOfCoverage",
                     "-R",  fastaFile,
                     "-I", bamFiles[i],
                     "-O",  file.path(depthCoverageFolder,paste0(bam[i], ".DATA")),
                     "--output-format", "TABLE",
                     " --disable-sequence-dictionary-validation true ",
                     "-L", bedFile)
        #"-L", paste0(outputFolder,"/list.interval_list"))
        paste(cmd);system(cmd)
      }
      #   print(basename(bam))
      #   sample_name<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bam))
      #   cmd_read_counts <- paste0( " gatk DepthOfCoverage ",
      #                              " -R ", fastaFile,
      #                              " --disable-sequence-dictionary-validation true ",
      #                              " -I ", bamsDir,"/",sample_name, ".bam",
      #                              " -L ", bedFile,
      #                              " -O /", coverageFiles)
      #
      #   paste(cmd_read_counts);system(cmd_read_counts);
    }

    print("run atlasCNV")

    setwd(outputFolder)

    #read config file
    print(atlascnvFolder)
    config <- read.table(paste0(atlascnvFolder,"/config"))

    #Put the Depth of Coverage folder in config files
    config[1,1] <- paste0("GATKDIR=", depthCoverageFolder, "/[SAMPLE_FCLBC].DATA.sample_interval_summary")

    write.table(config, paste0(outputFolder,"/config"), col.names = FALSE, row.names = FALSE, quote = FALSE)


    # Run Atlas-CNV ----

    cmd_run_atlascnv <- paste0( "perl ",atlascnvFolder,"/atlas_cnv.pl",
                                " --config ", outputFolder,"/config",
                                " --panel ", outputFolder,'/panel',name,'.txt',
                                " --sample ", outputFolder,"/", name, ".sample")

    paste(cmd_run_atlascnv);system(cmd_run_atlascnv);

    #merge cnv output files and modify columns
    #resDF <- list.files(paste0(outputFolder,"/",name), pattern = ".cnv$", recursive = TRUE, full.names = TRUE)  %>%
    resDF <- list.files(name, pattern = ".cnv$", recursive = TRUE, full.names = TRUE)  %>%
      purrr::set_names() %>%
      purrr::map_dfr(read.csv, .id = "sample", sep='\t' ) %>%
      mutate(sample = str_replace_all(basename(sample),".cnv","")) %>%
      rename( gene= Gene_Exon) %>%
      tidyr::separate(Exon_Target, c("chr", "start","end"))%>%
      mutate(CNV.type = ifelse(cnv == "del",
                               "deletion",
                               "duplication"))


    # Read output file, add CNV.type column and write the results in a tsv file
    write.table(resDF, file.path(outputFolder, "cnv_calls.tsv"), sep="\t", quote=F, row.names = FALSE, col.names = TRUE)

    # Save results----
    # Path to  tsv file
    finalSummaryFile <- file.path(outputFolder, "cnv_calls.tsv")
    # Save results in a GenomicRanges object
    message("Saving CNV GenomicRanges results")
    saveResultsFileToGR(outputFolder, basename(finalSummaryFile), geneColumn = "gene",
                        sampleColumn = "sample", chrColumn = "chr", startColumn = "start",
                        endColumn = "end", cnvTypeColumn = "CNV.type")

    #Delete temporary files if specified
    if(includeTempFiles == "false"){
      filesAll <- list.files(outputFolder, full.names = TRUE)
      filesToKeep <- c("failedRois.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
      filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse= "|"), filesAll, value=TRUE))])
      do.call(unlink, filesToRemove)
    }


  }

  setwd("../..")
}

# Create output folder ----
# if (!is.null(params$outputFolder)) {
#   outputFolder <- params$outputFolder
# } else {
#   outputFolder <- file.path(getwd(), "output", paste0("atlasCNV-", name))
# }
# unlink(outputFolder, recursive = TRUE);
# dir.create(outputFolder, showWarnings = FALSE)
#


print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)

