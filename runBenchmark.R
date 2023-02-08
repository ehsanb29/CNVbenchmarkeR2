# Description: Run the benchmark analysis
# USAGE: Rscript runBenchmark.R [-t tools_file] [-d datasets_file]

# libs
library("optparse")
library(yaml)

# Check runBenchmark.R is being called from CNVbenchmarkeR2 folder
if (length(list.files(pattern = "runBenchmark.R"))== 0){
  cat("Sorry, runBenchmark.R should be called from CNVbenchmarkeR2 folder\n")
  quit()
}

# Build options list
option_list <- list(
  make_option(c("-t", "--tools"), type="character", default="tools.yaml", 
              help="Path to tools file (yaml)", metavar="character"),
  make_option(c("-d", "--datasets"), type="character", default="datasets.yaml", 
              help="Path to datasets file (yaml)", metavar="character")
); 
opt_parser <- OptionParser(option_list=option_list);

# Load params
args <- parse_args(opt_parser);
tools <- yaml.load_file(args$tools)


# create logs/output folder if not exists
dir.create("logs", showWarnings = F)
dir.create("output", showWarnings = F)


## Execute tools on selected datasets ##

# Panelcn.mops
if (tools$panelcn == T){
  cat(as.character(Sys.time()), " - Executing panelcn.MOPS\n")
  cmd <- paste("Rscript tools/panelcnmops/runPanelcnmops.r tools/panelcnmops/panelcnmopsParams.yaml", args$datasets, " > logs/panelcnmops.log 2>&1")
  system(cmd)
}

# Decon
if (tools$decon == T){
  cat(as.character(Sys.time()), " - Executing DECoN\n")
  cmd <- paste("Rscript tools/decon/runDecon.r tools/decon/deconParams.yaml", args$datasets, " > logs/decon.log 2>&1")
  system(cmd)
}

# ExomeDepth
if (tools$panelcn == T){
  cat(as.character(Sys.time()), " - Executing ExomeDepth\n")
  cmd <- paste("Rscript tools/exomedepth/runExomedepth.r tools/exomedepth/exomedepthParams.yaml", args$datasets, " > logs/exomedepth.log 2>&1")
  system(cmd)
}

# CODEX2
if (tools$panelcn == T){
  cat(as.character(Sys.time()), " - Executing CODEX2\n")
  cmd <- paste("Rscript tools/codex2/runCodex2.r tools/codex2/codex2Params.yaml", args$datasets, " > logs/codex2.log 2>&1")
  system(cmd)
}

# ClinCNV
if (tools$clincnv == T){
  cat(as.character(Sys.time()), " - Executing ClinCNV\n")
  cmd <- paste("Rscript tools/clincnv/runClincnv.R tools/clincnv/clincnvParams.yaml", args$datasets, " > logs/clincnv.log 2>&1")
  system(cmd)
}

# CoNVaDING
if (tools$clincnv == T){
  cat(as.character(Sys.time()), " - Executing CoNVaDING\n")
  cmd <- paste("Rscript tools/convading/runConvading.r tools/convading/convadingParams.yaml", args$datasets, " > logs/convading.log 2>&1")
  system(cmd)
}

##  Generate summary file  ##

cat(as.character(Sys.time()), " - Generating summary file")
cmd <- paste("Rscript utils/summary.r", args$tools, args$datasets, "> logs/summary.log 2>&1")
system(cmd)
