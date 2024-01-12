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
      params <- yaml.load_file(paste0("tools/", algName, "/", algName, "Params.yaml"))

      #Create Depth of Coverage files for Viscap and AtlasCNV
      if(algName == "viscap" | algName == "atlasCNV"){
        depthCoverageFolder <- paste0("./evaluate_parameters/", algName, "/", name, "/DepthOfCoverage")

        if(!dir.exists(depthCoverageFolder)| dir.exists(depthCoverageFolder) &&  length(list.files(depthCoverageFolder))<7){
          bamsDir <- file.path(dataset$bams_dir)
          bedFile <- file.path(dataset$bed_file)
          fastaFile <- file.path(dataset$fasta_file)
          gatkFolder<- params$gatkFolder

          bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
          dir.create(depthCoverageFolder, showWarnings = FALSE)
          bam <- basename(bamFiles) %>% tools::file_path_sans_ext()
          for (i in seq_len(length(bamFiles))){
            cmd <- paste(file.path(gatkFolder, "gatk"),  "DepthOfCoverage",
                         "-R",  fastaFile,
                         "-I", bamFiles[i],
                         "-O",  ifelse(algName == "viscap",
                                       file.path(depthCoverageFolder,bam[i]),
                                       file.path(depthCoverageFolder,paste0(bam[i], ".DATA"))),
                         "--output-format", "TABLE",
                         #" --disable-sequence-dictionary-validation true ",
                         "-L",  bedFile )
            paste(cmd);system(cmd)
          }
          }
      }

      #Common part for GermlineCNV
      if(algName == "germlineCNVcaller"){
        preCal <- paste0("./evaluate_parameters/", algName, "/", name, "/preCalculated")
        readCountsFolder <- paste0(preCal, "/ReadCounts")
        intervalListFolder <- paste0(preCal, "/IntervalList")
        contigPloidy <- paste0(preCal, "/contigPloidy")

        #get bams
        bamsDir <- file.path(dataset$bams_dir)
        bamFiles <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)

        if(!dir.exists(readCountsFolder)| dir.exists(readCountsFolder) &&  length(list.files(readCountsFolder))!=length(bamFiles)){

          #create directories
          dir.create(preCal, showWarnings = FALSE)
          dir.create(readCountsFolder, showWarnings = FALSE)
          dir.create(intervalListFolder, showWarnings = FALSE)
          dir.create(contigPloidy, showWarnings = FALSE)



          #get singularity folder
          singularityFolder <- params$singularityFolder

          #get bed file
          bedFile <- file.path(dataset$bed_file)

          #Create fasta dict
          fastaFile <- file.path(dataset$fasta_file)
          fastaDict <- paste0(tools::file_path_sans_ext(fastaFile), ".dict")
          picardFolder <- params$picardFolder

          if(!file.exists(fastaDict)){
            cmd <- paste0(" java -jar ", picardFolder,
                          " CreateSequenceDictionary",
                          " -R ", fastaFile,
                          " -O ", fastaDict)


            print(cmd);system(cmd);
          }

          #create intervals file (.intervals_list) from .bed file
          interval_list<-paste0(intervalListFolder,"/list.interval_list")
          if(!file.exists(interval_list)){
            cmd_interval <- paste0( " java -jar ", picardFolder,
                                    " BedToIntervalList",
                                    " -I ", bedFile,
                                    " -O ", interval_list,
                                    " -SD ", fastaDict)

            print(cmd_interval);system(cmd_interval);
          }
          for (bam in bamFiles){
            print(basename(bam))
            sample_name<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bam))
            cmd_read_counts <- paste0( "apptainer run ",
                                       " -B ", bamsDir, " ",
                                       singularityFolder,
                                       ' gatk --java-options "-Xmx30G" CollectReadCounts',
                                       " -I ", bamsDir,"/", sample_name, ".bam",
                                       " -L ", interval_list,
                                       " --interval-merging-rule OVERLAPPING_ONLY",
                                       " -O ", readCountsFolder,"/", sample_name, ".counts.hdf5")
            print(cmd_read_counts);system(cmd_read_counts);
          }

          #GATK: Determine contig ploidy ----
          #merge all hdf5 files forrunning DetermineGermlineContigPloidy
          hdf5s <- list.files(readCountsFolder, pattern = '*.hdf5$', full.names = TRUE)
          all_hdf5_names<- paste(hdf5s, collapse=" -I ")
          contigFile <- file.path(params$contigFile)
          cmd_det_contig_ploidy <- paste0( "apptainer run ",
                                           #"-B ", currentFolder, " ",
                                           "-B  /data ",
                                           singularityFolder,
                                           ' gatk --java-options "-Xmx30G" DetermineGermlineContigPloidy',
                                           " -I ", all_hdf5_names,
                                           #" -L ", interval_list,
                                           " --contig-ploidy-priors ", contigFile,
                                           " -O ", contigPloidy,
                                           " --output-prefix contig")

          print(cmd_det_contig_ploidy);system(cmd_det_contig_ploidy);
        }
      }

      filesParams <- list.files(paste0(getwd(), "/evaluate_parameters/", algName, "/", name), pattern= "params.yaml", recursive=TRUE, full.names=TRUE)


      #Run benchmark for every param

      for(j in seq_len(length(filesParams))){

      print(j)
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
                     "true",
                     file.path(paste0(filesParams[j] %>% dirname() %>% dirname(),"/logs"),paste0(algName,".log")))

      cmd <- paste("sbatch evaluate_parameters/job.sh", query)
      print(cmd, quote=FALSE)
      system(cmd)

        }
    }
}
  }



}

