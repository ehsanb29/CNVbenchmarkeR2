# Runs Atlas-CNV over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runAtlasCNV.R [atlasCNVparams_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions
library(dplyr)

# Read args ----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  paramsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  paramsFile <- "params.yaml"
  datasetsParamsFile <- "../../datasets.yaml"
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


#create input files required for running GATK
# Dictionary ----
#create dictionary file (.dict) from reference genome (.)
if(!file.exists(paste0(atlascnvFolder, "/reference.dict"))){
  cmd_dictionary<- paste0(" java -jar $PICARD",
                          " CreateSequenceDictionary",
                          " -R ", fastaFile,
                          " -O ", atlascnvFolder, "/reference.dict")
  
  
  paste(cmd_dictionary);system(cmd_dictionary);
}
#create the variable for the reference dictionary
dictionary<-paste0(atlascnvFolder,"/reference.dict")

#create the folder for the coverage files

coverageFiles <- file.path(paste0(atlascnvFolder, "/coverage_files"))
if (!file.exists(coverageFiles)){
  dir.create(coverageFiles)
} else {
  unlink(coverageFiles, recursive = TRUE)
  dir.create(coverageFiles, recursive=TRUE)
}


# Dataset iteration ----
# go over datasets and run cnvkit for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting Atlas-CNV for", name, "dataset", sep=" "))
    
    
    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
    bedFile <- file.path(dataset$bed_file)
    #annotatedBedFile<-file.path(dataset$annotated_bed_file)
    
    # Input files for Atlas-CNV ----
    ## Interval list----
    #create intervals file (.intervals_list) from .bed file
    if(!file.exists(paste0(atlascnvFolder, "/list.interval_list"))){
      cmd_interval <- paste0( " java -jar  $PICARD",
                              " BedToIntervalList",
                              " -I ", bedFile,
                              " -O ", atlascnvFolder,"/list.interval_list",
                              " -SD ", dictionary)
      paste(cmd_interval);system(cmd_interval); }
    

    ## Panel file ---- 
    #create panel file to input in Atlas
    ##
    tempBed<-read.csv(bedFile, sep="\t", header = FALSE)
    
    #tempBed<-read.csv(bedFile, sep="\t", header = FALSE)
    
    panel<-tempBed %>% 
      mutate(Exon_Target = paste0(V1,":", V2,"-", V3))  %>% 
      mutate(Gene_Exon=V4) %>%
      mutate(Call_CNV = ifelse(grepl("^Y", Exon_Target), "N", "Y"))%>%
      mutate(RefSeq=V4) %>%
      select(-starts_with("V")) 
    #export panel file
    write.table(panel, file=paste0(atlascnvFolder,'/panel',name,'.txt'),sep="\t", col.names=T, row.names=FALSE,quote=FALSE)
     
    ## Sample file ----  
    
    tempSample<-as.data.frame(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bamFiles)))
    sample<-tempSample %>%
      mutate(sex="F") %>%
      mutate(mp=name)
    #export sample file
    write.table(sample, file=paste0(atlascnvFolder,'/',name,'.sample'),sep="\t", col.names=F, row.names=FALSE,quote=FALSE)
    rm(tempSample, tempBed)
    
    #create folder for coverage files
    

    #get depth of coverage for each BAM file, using GATK
    for (bam in bamFiles){
      print(basename(bam))
      sample_name<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bam))
      cmd_read_counts <- paste0( " gatk DepthOfCoverage ",
                                 " -R ", fastaFile,
                                 " --disable-sequence-dictionary-validation true ",
                                 " -I ", bamsDir,"/",sample_name, ".bam",
                                 " -L ", bedFile,
                                 " -O /", coverageFiles)
      paste(cmd_read_counts);system(cmd_read_counts);
    }
    
    # Run Atlas-CNV ----

    cmd_run_atlascnv <- paste0( "perl ",atlascnvFolder,"/atlas_cnv.pl",
                               " --config ",atlascnvFolder,"/config",
                               " --panel ",atlascnvFolder,'/panel',name,'.txt',
                               " --sample ",atlascnvFolder,"/", name, ".sample")
    
    paste(cmd_run_atlascnv);system(cmd_run_atlascnv);

    #merge cnv output files and modify columns
    resDF <- list.files(paste0(atlascnvFolder,"/",name), pattern = ".cnv$", recursive = TRUE, full.names = TRUE)  %>% 
      set_names() %>% 
      map_dfr(read.csv, .id = "sample", sep='\t' ) %>%
      mutate(sample = str_replace_all(basename(sample),".cnv","")) %>%
      rename( gene= Gene_Exon) %>%
      separate(Exon_Target, c("chr", "start","end"))%>% 
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
    
    
    
    
  }  
  
}

# Create output folder ----
if (!is.null(params$outputFolder)) {
  outputFolder <- params$outputFolder
} else {
  outputFolder <- file.path(getwd(), "output", paste0("atlasCNV-", name))
}
unlink(outputFolder, recursive = TRUE);
dir.create(outputFolder, showWarnings = FALSE)



print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)
  
  