# 21 april

# Runs GermlineCNVcaller over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runGermlineCNVcaller.R [GermllineCNVcaller_params_file] [datasets_params_file]
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

#Load the parameters file  
params <- yaml.load_file(paramsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# print params
print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))

#locate the folder and the contig file
germlinecnvcallerFolder<- file.path(params$germlinecnvcallerFolder)
fastaFile <- file.path(datasets$ICR96$fasta_file)
contigFile <- file.path(params$contigFile)


#create input files required for running GATK
# Dictionary ----
#create dictionary file (.dict) from reference genome (.)
if(!file.exists(paste0(germlinecnvcallerFolder, "/reference.dict"))){
  cmd_dictionary<- paste0(" java -jar $PICARD",
                          " CreateSequenceDictionary",
                          " -R ", fastaFile,
                          " -O ", germlinecnvcallerFolder, "/reference.dict")
  
  
  paste(cmd_dictionary);system(cmd_dictionary);
}
#create the variable for the reference dictionary
dictionary<-paste0(germlinecnvcallerFolder,"/reference.dict")


# Dataset iteration ----

# go over datasets and run cnvkit for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting GermlineCNVcaller for", name, "dataset", sep=" "))
    
    
    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    annotatedBedFile<-file.path(dataset$annotated_bed_file)
    

    # Interval list----
    #create intervals file (.intervals_list) from .bed file
    if(!file.exists(paste0(germlinecnvcallerFolder, "/list.interval_list"))){
      cmd_interval <- paste0( " java -jar $PICARD",
                              " BedToIntervalList",
                              " -I ", bedFile,
                              " -O ", germlinecnvcallerFolder,"/list.interval_list",
                              " -SD ", dictionary)
      
      paste(cmd_interval);system(cmd_interval);
    }
    #create the variable for the interval list
    interval_list<-paste0(germlinecnvcallerFolder,"/list.interval_list")
    
    # Get bam files
    bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)

    


    # create GATKoutput folder if it doesnt exist
    gatkOutput <- file.path(paste0(germlinecnvcallerFolder, "/gatkOutput"))
    if (!file.exists(gatkOutput)){
      dir.create(gatkOutput)
      
      for (bam in bamFiles){
        print(basename(bam))
        sample_name<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bam))
        cmd_read_counts <- paste0( " gatk CollectReadCounts ",
                                   " -I ", outputFolder,"/",sample_name, ".bam",
                                   " -L ", interval_list,
                                   " --interval-merging-rule OVERLAPPING_ONLY",
                                   " -O /", gatkOutput,"/", sample_name, ".counts.hdf5")
        
        paste(cmd_read_counts);system(cmd_read_counts);
      }  
    }
    
    
    #GATK: Determine contig ploidy ----
    #merge all hdf5 files forrunning DetermineGermlineContigPloidy
    hdf5s <- list.files(gatkOutput, pattern = '*.hdf5$', full.names = TRUE)
    all_hdf5_names<- paste(basename(hdf5s), collapse=" -I ")
    
    contigPloidy <- file.path(paste0(gatkOutput, "/contigPloidy"))
    if (!file.exists(contigPloidy)){
      dir.create(contigPloidy)
      #contigPloidy<-paste0(gatkOutput,"/contigPloidy")
      cmd_det_contig_ploidy <- paste0( " gatk DetermineGermlineContigPloidy",
                                       " -I ",all_hdf5_names, 
                                       " --contig-ploidy-priors ", contigFile,
                                       " -O ", contigPloidy,
                                       " --output-prefix contig")
      
      paste(cmd_det_contig_ploidy);system(cmd_det_contig_ploidy);
    }
    
    
    #GATK: Run GermlineCNVcaller ----
    cnvRaw <- file.path(paste0(gatkOutput, "/cnvRaw"))
    if (!file.exists(cnvRaw)){
      dir.create(cnvRaw)
    } else {
      unlink(cnvRaw, recursive = TRUE)
      dir.create(cnvRaw, recursive=TRUE)
    }
    
    cnvRaw<-paste0(gatkOutput,"/cnvRaw")
    cmd_GermlineCNVCaller<- paste0( " srun --mem 10G gatk GermlineCNVCaller ",
                                    " --run-mode COHORT",
                                    " -L ", interval_list,
                                    " --interval-merging-rule OVERLAPPING_ONLY",
                                    "--contig-ploidy-calls",
                                    " -I ",all_hdf5_names, 
                                    " --contig-ploidy-priors ", contigPloidy, "/contig-calls",
                                    " -O ", cnvRaw,
                                    " --output-prefix cohort_run")
    
    paste(cmd_GermlineCNVCaller);system(cmd_GermlineCNVCaller);


#GATK: Run PostprocessingGermlineCNVcaller
    
    cnvProcessed <- file.path(paste0(gatkOutput, "/cnvProcessed"))
    if (!file.exists(cnvProcessed)){
      dir.create(cnvProcessed)
    } else {
      unlink(cnvProcessed, recursive = TRUE)
      dir.create(cnvProcessed, recursive=TRUE)
    }
    
    samples<-list.dirs(cnvRaw, "/cohort_run-calls")
    callsFolders<- paste0(samples, collapse=" -contig-ploidy-calls ")
    name<- paste0(basename(samples), collapse=" -I ")
    
    cmd_GermlineCNVCaller<- paste0( " gatk PostprocessGermlineCNVCalls",
                                    " --calls-shard-path ", callsFolders,
                                    " --contig-ploidy-calls ", contigPloidy,
                                    " --model-shard-path ", cnvRaw, "/cohort_run-model",

                                    " --output-genotyped-intervals postprocessing_GCC/sample_0_genotyped_intervals.vcf ",

                                    " --output-genotyped-segments postprocessing_GCC/sample_0_genotyped_segments.vcf ",
                                    " --output-denoised-copy-ratios postprocessing_GCC/sample_0_denoised_copy_ratios.tsv",

                                    " -O ", cnvProcessed,
                                    " --output-prefix cohort_run")
    
    
    #for each folder in cnvProcessed
    

    #format them as benchmark desired output
    
    
    
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
    
  }
}
  
  



# Create output folder
if (!is.null(params$outputFolder)) {
  outputFolder <- params$outputFolder
} else {
  outputFolder <- file.path(getwd(), "output", paste0("GermlineCNVcaller-", name))
}
unlink(outputFolder, recursive = TRUE);
dir.create(outputFolder, showWarnings = FALSE)



print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)