# Description: Run the benchmark analysis
# USAGE: Rscript runBenchmark.R [-t tools_file] [-d datasets_file] [-tf include_temp_files]

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
              help="Path to datasets file (yaml)", metavar="character"),
  make_option(c("-tf", "--include_temp_files"), type="character", default="true",
              help="Include temporary files in the output folders", metavar="character")
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
  cmd <- paste("Rscript tools/panelcnmops/runPanelcnmops.r tools/panelcnmops/panelcnmopsParams.yaml", args$datasets, args$include_temp_files, " > logs/panelcnmops.log 2>&1")
  system(cmd)
}

# Decon
if (tools$decon == T){
  cat(as.character(Sys.time()), " - Executing DECoN\n")
  cmd <- paste("Rscript tools/decon/runDecon.r tools/decon/deconParams.yaml", args$datasets, args$include_temp_files, " > logs/decon.log 2>&1")
  system(cmd)
}

# ExomeDepth
if (tools$exomedepth == T){
  cat(as.character(Sys.time()), " - Executing ExomeDepth\n")
  cmd <- paste("Rscript tools/exomedepth/runExomedepth.r tools/exomedepth/exomedepthParams.yaml", args$datasets, args$include_temp_files," > logs/exomedepth.log 2>&1")
  system(cmd)
}

# CODEX2
if (tools$codex2 == T){
  cat(as.character(Sys.time()), " - Executing CODEX2\n")
  cmd <- paste("Rscript tools/codex2/runCodex2.r tools/codex2/codex2Params.yaml", args$datasets, args$include_temp_files," > logs/codex2.log 2>&1")
  system(cmd)
}

# ClinCNV
if (tools$clincnv == T){
  cat(as.character(Sys.time()), " - Executing ClinCNV\n")
  cmd <- paste("Rscript tools/clincnv/runClincnv.R tools/clincnv/clincnvParams.yaml", args$datasets, args$include_temp_files, " > logs/clincnv.log 2>&1")
  system(cmd)
}

# CoNVaDING
if (tools$clincnv == T){
  cat(as.character(Sys.time()), " - Executing CoNVaDING\n")
  cmd <- paste("Rscript tools/convading/runConvading.r tools/convading/convadingParams.yaml", args$datasets, args$include_temp_files," > logs/convading.log 2>&1")
  system(cmd)
}
# Cobalt
if (tools$cobalt == T){
  cat(as.character(Sys.time()), " - Executing Cobalt\n")
  cmd <- paste("Rscript tools/cobalt/runCobalt.r tools/cobalt/cobaltParams.yaml", args$datasets, args$include_temp_files," > logs/cobalt.log 2>&1")
  system(cmd)
}

# ClearCNV
if (tools$clearCNV == T){
  cat(as.character(Sys.time()), " - Executing clearCNV\n")
  cmd <- paste("Rscript tools/clearCNV/runClearCNV.r tools/clearCNV/clearCNVParams.yaml", args$datasets, args$include_temp_files," > logs/clearCNV.log 2>&1")
  system(cmd)
}

# Cnvkit
if (tools$cnvkit == T){
  cat(as.character(Sys.time()), " - Executing cnvkit\n")
  cmd <- paste("Rscript tools/cnvkit/runCnvkit.r tools/cnvkit/cnvkitParams.yaml", args$datasets, args$include_temp_files," > logs/cnvkit.log 2>&1")
  system(cmd)
}

# GermlineCNVcaller
if (tools$germlineCNVcaller == T){
  cat(as.character(Sys.time()), " - Executing germlineCNVcaller\n")
  cmd <- paste("Rscript tools/germlineCNVcaller/runGermlinecnvcaller.r tools/germlineCNVcaller/germlinecnvcallerParams.yaml", args$datasets, args$include_temp_files," > logs/germlineCNVcaller.log 2>&1")
  system(cmd)
}

# AtlasCNV
if (tools$atlasCNV == T){
  cat(as.character(Sys.time()), " - Executing Atlas-CNV\n")
  cmd <- paste("Rscript tools/atlasCNV/runAtlascnv.r tools/atlasCNV/atlasCNVParams.yaml", args$datasets, args$include_temp_files,">  logs/atlasCNV.log 2>&1")
  system(cmd)
}




##  Generate summary file  ##

cat(as.character(Sys.time()), " - Generating summary file")
cmd <- paste("Rscript utils/summary.r", args$tools, args$datasets, "> logs/summary.log 2>&1")
system(cmd)
