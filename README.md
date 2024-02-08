# CNVbenchmarkeR2 #

CNVbenchmarkeR2 is a framework to benchmark germline copy number variant (CNV) calling tools on different NGS datasets. Current version supports DECoN, CoNVaDING, panelcn.MOPS, ExomeDepth, CODEX2, ClinCNV, clearCNV, GATK-gCNV, Atlas-CNV, Cobalt, CNVkit and VisCap tools

Previous version, CNVbenchmarkeR, is available [here](https://github.com/TranslationalBioinformaticsIGTP/CNVbenchmarkeR).


### Prerequisites ###

Tools should be properly installed. Links for tools installation:

- [Panelcn.mops](https://github.com/bioinf-jku/panelcn.mops)
- [CoNVaDING](https://molgenis.gitbooks.io/convading)
- [DECoN](https://github.com/RahmanTeam/DECoN)
- [CODEX2](https://github.com/yuchaojiang/CODEX2)
- [ExomeDepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html)
- [ClinCNV](https://github.com/imgag/ClinCNV)
- [clearCNV](https://github.com/bihealth/clear-cnv)
- [Atlas-CNV](https://github.com/theodorc/Atlas-CNV)
- [Cobalt](https://github.com/ARUP-NGS/cobalt)
- [CNVkit](https://github.com/etal/cnvkit)
- [VisCap](https://github.com/pughlab/VisCap)
- [GATK](https://hub.docker.com/r/broadinstitute/gatk) (At the moment, it only works using [singularity]https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

Also, R/Bioconductor should be installed with at least this packages: GenomicRanges, biomaRt, regioneR.

### How to use
1. Get Code
```
git clone https://github.com/jpuntomarcos/CNVbenchmarkeR2 
```

2. **Configure tools.yaml** to set which tools will be benchmarked. For DECoN, Convading, ClinCNV, clearCNV, GATK-gCNV, Atlas-CNV, Cobalt, CNVkit and Viscap modify the corresponding params.yaml file stored in:
tools/[name_of_the_tool]/[name_of_the_tool]Params.yaml 
to include the tool foler installation. 

Additional specification to consider: 

- For ClinCNV, set ngsbits folder
- For Atlas-CNV, set the contigFile folder, the gatkFolder and the Rscript and rpath folder. 
- For VisCap, set the gatkFolder and picardJar folder




3. **Configure datasets.yaml** to define against which datasets the tools will be executed. Within this file, it is important to provide files with the exact expected format (**special attention** to `validated_results_file` and `bed_file` that are **tab-delimited** files). To do so, please **check the [examples](https://github.com/jpuntomarcos/CNVbenchmarkeR2/tree/master/examples) folder**.


4. Launch CNVbenchmarkeR2
```
cd CNVbenchmarkerR2
Rscript runBenchmark.R [-t tools_yaml] [-d datasets_yaml] [-f include_temp_files]
```


### Output ###

A summary file and a .csv results file will be generated at output/summary folder. Stats include sensitivity, specificity, no-call rate, precision (PPV), NPV, F1, MCC and kappa coefficient.

Statistics are calculated per ROI, per gene and at whole strategy level.

Logs files will be generated in the logs folder. Output for each tool and dataset will be generated at output folder.


### Troubleshooting  ###

Two important checks to ensure that metrics are computed correctly:

- The **sample names in the `validated_results_file` should match the file names of your bam files** (excluding the .bam extension). For example, if the `validated_results_file` contains sample names like mySample2312, your bam files should have file names like mySample2312.bam .
- Provide and use chromosomes names with the same format, that is, do not use "chr5" and "5" in you bed and `validated_results_file` files, for example.


## Extra feature: evaluate parameters ##
A parameter evaluator is also attached in the framework. It executes each tool parameter over a broad range of values to assess its impact on tool performance. Up to 15 values are evaluated for each numerical param and all the available options for categorical ones. 
The parameter evaluator empowers users to gain deeper insights into how individual parameters influence the overall performance of CNV calling tools.


### Prerequisites ###

An SGE cluster system has to be available.

### Run evaluate parameters
1. Configure all the steps needed to run the benckmark (expained in section how to use)

2. Run the evaluate script to create the necessary files for 'runEvaluate.R' . This step must be executed within the evaluate_parameters folder.
```
cd evaluate_parameters
Rscript evaluate.R [-t tools_yaml] [-d datasets_yaml]
```
3. Execute the runEvaluate script. This step should be performed in the CNVbenchmarkeR2 folder.
```
cd ..
Rscript evaluate_parameters/runEvaluate.R 

```
For space optimization, it's recommended to set the -f parameter to false, which deletes all intermediate files.




