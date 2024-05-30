# nf-core/evoverse: Output
## Introduction

This document describes the output produced by the pipeline. All plots generated in each step are summarised into the final report.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to both published and unpublished results. 



## Pipeline overview

<!-- The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Variant Annotation](#variant-annotation) - annotation of variants and cohort summary visualization
* [Formatter](#formatter) - coversion of files to different formats (unpubblished)
* [Lifter](#lifter) - pileup of private mutations of the other samples in multi-sample setting
* [Driver Annotation](#driver-annotation) - _add description_ (unpubblished)
* [QC](#qc) - quality control of copy-number and somatic mutation calling and creation of multi-CNAqc object
* [Subclonal Deconvolution](#subclonal-deconvolution) - 
* [Signature Deconvolution](#signature-deconvolution) - 

Intermediate steps of the pipeline will output unpublished results which will be available for the user in the working directory of the pipeline. -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and consists in five main subworkflows:

* [Variant Annotation](#variant-annotation)
* [QC](#qc)
* [Lifter](#lifter)
* [Signature Deconvolution](#signature-deconvolution)
* [Subclonal Deconvolution](#subclonal-deconvolution)

Intermediate steps connetting the main subworkflows will output [unpublished results](#unpublished-results) which will be available in the working directory of the pipeline. These steps consist in:

* [Formatter](#formatter)
* [Lifter](#lifter)
* [Driver Annotation](#driver-annotation)

## Variant Annotation

This directory contains results from the variant annotation subworkflow. At the level of individual samples, genomic variants present in the input VCF files are annotated using [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and further converted into MAF format using [vcf2maf](https://github.com/mskcc/vcf2maf). Genomic information of samples from the same cohort are summarized into a unique MAF object with [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html).

### VEP
[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html) is a `Ensembl` tool that determines the effect of all variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence.
This step starts from VCF files.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/variant_annotation/VEP/<dataset>/<patient>/<sample>/`**

- `<sample>.vep.summary.html`
  - Summary of the VEP run to be visualised with a web browser
- `<sample>.vep.vcf.gz`
  - annotated VCF file

</details>


### vcf2maf

[vcf2maf](https://github.com/mskcc/vcf2maf) convert a VCF file into a MAF (Mutation Annotation Format) file, where each variant is mapped to only one of all possible gene transcripts/isoforms that it might affect. vcf2maf is designed to work with VEP annotation output. 

> **NB:** While VEP is tolerant of chromosome format mismatches (when the input .VCF file uses the UCSC format chrN and the reference fasta uses Ensembl/NCBI format N), vcf2maf is not.
>  Make sure the reference fasta chromosome format matches that of your input

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/variant_annotation/vcf2maf/<dataset>/<patient>/<sample>/`**
- `<sample>.vcf2maf.maf`
  - annotated MAF file

</details>


### maftools

[maftools]( https://bioconductor.org/packages/release/bioc/html/maftools.html) is an R-based tool that provides a comprehensive set of functions for processing MAF files and to perform most commonly used analyses in cancer genomics. In particular, maftools summarize, analyze and visualize MAF files from the same cohort providing summary statistics of the top mutated genes and aggregating metadata into a single oncoplot visualization.
This step started from annotated MAF files.

MAF fields requirements:

- Mandatory fields: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type` and `Tumor_Sample_Barcode`.

- Optional fields: `VAF` (Variant Allele Frequency), `amino acid change` information.

<details markdown="1">
<summary>Output files for the dataset</summary>

**Output directory: `{publish_dir}/variant_annotation/maftools/<dataset>/`**
- `<dataset>.maftools.rds`
  - summarized MAF object
- `<dataset>.maftools.pdf`
  - summary plots (oncoplot, summary stats)

</details>

## Lifter

The Lifter subworkflow is an optional step and it is run when `--mode multisample` is used. When multiple samples from the same patient are provided, the user can specify either a single joint VCF file, containing variant calls from all tumor samples of the patient, or individual sample specific VCF files. In the latter case, path to tumor BAM files must be provided in order to collect all mutations from the samples and perform pile-up of sample's private mutations in all the other samples. Two intermediate steps, [get_positions](#get_positions) and [mpileup](#mpileup), are performed to identify private mutations in all the samples and retrieve their variant allele frequency. Once private mutations are properly defined, they are merged back into the original VCF file during the [join_positions](#join_positions) step. The updated VCF file is then converted into a `vcfR` RDS object.

The output files of [get_positions](#get_positions) and [mpileup](#mpileup) are intermediate and by default not kept in the result directory. _add the possibility to save also this files in the results dir??_ 

### join_positions
In this step, all retrieved mutations are joined with original mutations present in input VCF, which is in turn converted into an RDS object using [vcfR](https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/lifter/mpileup/<dataset>/<patient>/<sample>/`**
- `pileup_VCF.rds`
  - RDS containing shared and private mutations

</details>

## QC
The QC subworkflows requires in input a segmentation file from allele-specific copy number callers (either [Sequenza](https://sequenzatools.bitbucket.io/#/home), [ASCAT](https://github.com/VanLoo-lab/ascat)) and the joint VCF file from [join_positions](#join_positions) subworkflow. The QC sub-workflows first conduct quality control on CNV and somatic mutation data for individual samples in [CNAqc](#cnaqc) step, and subsequently summarize validated information at patient level in [join_CNAqc](#join_cnaqc) step.
The QC subworkflow is a crucial step of the pipeline as it ensures high confidence in identifying clonal and subclonal events while accounting for variations in tumor purity.

### CNAqc
[CNAqc](https://caravagnalab.github.io/CNAqc/) is a package to quality control (QC) bulk cancer sequencing data for validating copy number segmentations against variant allele frequencies of somatic mutations.

<details markdown="1">
<summary>Output files for all samples</summary>

<!-- **Output directory: `{outdir}/subclonal_deconvolution/mobster/<dataset>/<patient>/<sample>/`** -->
**Output directory: `{publish_dir}/QC/CNAqc/<dataset>/<patient>/<sample>/`**
- `data.pdf`
  - CNAqc report with genome wide mutation and allele specific copy number plots
- `qc.pdf`
  - QC step report resulting from peak analysis
- `qc.rds`
  -  CNAqc object
</details>

### join_CNAqc
This module creates a multi-CNAqc object for patient by summarizing the quality check performed at the single sample level. For more information about the strucutre of multi-CNAqc object see [CNAqc documentation]((https://caravagnalab.github.io/CNAqc/)).

<!-- **Output directory: `results/datasetID/patientID/join_CNAqc`**
* `multi_cnaqc.rds`
  * multi-CNAqc object  -->

<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{publish_dir}/QC/join_CNAqc/<dataset>/<patient>/`**
- `multi_cnaqc.rds`
  - multi-CNAqc object
</details>

## Signature Deconvolution

<!-- Mutational signatures represent characteristic patterns of somatic mutations in cancer genomes, reflecting the underlying mutational processes at the basis of tumor evolution and progression. Mutational signatures are discovered by analyzing ensemble point-mutation counts from a set of individual samples. Validated mutations from [join_CNAqc](#join_cnaqc) step are converted into a TSV joint table in (see [tsvparse](#tsvparse) module), subsequently given as input to signature deconvolution subworkflow, which performs de novo extraction, inference, deciphering or deconvolution of mutational counts.  -->

Mutational signatures are distinctive patterns of somatic mutations in cancer genomes that reveal the underlying mutational processes driving tumor evolution and progression. These signatures are identified by analyzing aggregated point-mutation counts from multiple samples. Validated mutations from the [join_CNAqc](#join_cnaqc) step are converted into a joint TSV table (see [tsvparse](#tsvparse)) and then input into the signature deconvolution subworkflow, which performs de novo extraction, inference, interpretation, or deconvolution of mutational counts.

The results of this step are collected in `{pubslish_dir}/signature_deconvolution/`. Two tools can be specified by using `--tools` parameter: [SparseSignatures](#sparsesignatures) and [SigProfiler](#sigprofiler).

### SparseSignatures

[SparseSignatures](https://www.bioconductor.org/packages/release/bioc/html/SparseSignatures.html) is an R-based computational framework that provides a set of functions to extract and visualize the mutational signatures that best explain the mutation counts of a large number of patients.

<details markdown="1">
<summary>Output files for dataset</summary>

**Output directory: `{publish_dir}/signatures_deconvolution/SparseSignatures/<dataset>/`**
- `<dataset>.SparseSig.rds`
  - signatures best configiration object
- `<dataset>.SparseSig.pdf`
  - signatures plot

</details>
  

### SigProfiler

[SigProfiler](https://osf.io/t6j7u/wiki/home/) is a python framework that allows _de novo_ extraction of mutational signatures from data generated in a matrix format. The tool identifies the number of operative mutational signatures, their activities in each sample, and the probability for each signature to cause a specific mutation type in a cancer sample. The tool makes use of `SigProfilerMatrixGenerator` and `SigProfilerPlotting`, seamlessly integrating with other `SigProfiler` tools.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/signatures_deconvolution/SigProfiler/<dataset>/`**
- `<dataset>.` _missing_
  - _missing_
</details>


## Subclonal Deconvolution

The subclonal deconvolution subworkflow requires in input a joint `mCNAqc` object resulting from the [join_CNAqc](#join_cnaqc) step. Two different modalities can be specified by the user using `--mode` parameter: single sample and multi sample mode.
<!-- If `--mode singlesample` is provided, each sample is analysed individually providing a snapshot of clonal and subclonal diversity starting from allele frequency of detected somatic variants. When multiple samples from the same patient are provided, the user may take advantage of the multisample modality, by setting `--mode multisample`. -->
<!-- This approach allows for a more detailed and accurate identification of subclonal populations, as it can capture spatial and temporal heterogeneity within the tumor. By integrating data across multiple samples, it improves the resolution of subclonal structures and provides insights into the evolutionary dynamics and progression of the tumor. -->

<!-- Various tools can be specified using the `--tools` parameter, leading to different methods for performing subclonal deconvolution analysis. Among the available tools, [MOBSTER](https://caravagnalab.github.io/mobster/) operates only on single samples but can still be specified for use in multi-sample mode, while [VIBER](https://caravagnalab.github.io/VIBER/index.html) and [PyClone-VI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03919-2) can operates in both modalities. In case `--mode multisample` and `--tools mobster,pyclone-vi` are specified, first MOBSTER is run on individual samples to remove tail mutations and then PyClone-VI operates a multivariate subclonal deconvolution on the preprocessed MOBSTER mutations. A similar procedure is perfomed when  `--mode multisample` and `--tools mobster,viber` are defined. More detailed explanation is provided in the following sections.  -->

The results of subclonal decovnultion step are collected in `{publish_dir}/subclonal_deconvolution/` directory.

### Single sample

If `--mode singlesample` is provided, each sample is analysed individually providing a snapshot of clonal and subclonal diversity starting from allele frequency of detected somatic variants. Available tools provides different ways of performing subclonal deconvolution. Both MOBSTER and VIBER model mutation counts as mixture of binomial distribution; however, MOBSTER includes a pareto Type-I power law distribution to model within-clone neutral dynamics. Instead, PyClone-VI models clonal architecture taking into account variant allele frequencies corrected for coincident copy number variation.

#### MOBSTER

[MOBSTER](https://caravagnalab.github.io/mobster/) processes mutant allelic frequencies to identify and remove neutral tails from the input data, so that subclonal reconstruction algorithms can be applied downstream to find subclones from the processed read counts.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/subclonal_deconvolution/mobster/<dataset>/<patient>/<sample>/`**

- `<sample>_fit.rds`
  - `.rds` object contains best fit of subclonal deconvolution
- `<sample>.pdf`

</details>

#### PyClone-VI

[PyClone-VI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03919-2) is a computationally efficient Bayesian statistical method for inferring the clonal population structure of cancers, by considering allele fractions and coincident copy number variation using a variational inference approach.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/subclonal_deconvolution/pyclone_vi/<dataset>/<patient>/<sample>/`**

- `<sample>_all_fits.h5`
  - HDF5 file for all possible fit
- `<sample>_best_fit.txt`
  - TSV file for the best fit
</details>

#### VIBER

[VIBER](https://caravagnalab.github.io/VIBER/index.html) is an R package that implements a variational Bayesian model to fit multi-variate Binomial mixtures. In the context of subclonal deconvolution in singlesample modality, VIBER models read counts that are associated with the most represented karyotype.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/subclonal_deconvolution/viber/<dataset>/<patient>/<sample>/`**

- `<sample>_best_st_fit.rds`
  - RDS file for best fit

</details>

### Multi samples

When multiple samples for the same patient are available (e.g. multi-regional or longitudinal experiment), we can be intersted in visualizing the tumor evolution at different level. In multi sample mode, a patient-specific subclonal deconvolution is perfomed.

> **NB:** At this stage, it is strongly recommended that the user conducts single-sample subclonal deconvolution using `mosbter` to eliminate mutations assigned to the neutral tail before proceeding with multi-variate clone identification.
> This process ensures a cleaned signal for downstream analyses aimed at focusing on functional intratumor heterogeneity.

#### Pyclone-VI

This folder contains the results of multivariate analysis using Pyclone-VI, which can be run prior to removal of mutations assigned to tail in all the samples.

<details markdown="1">
<summary>Output files for all patients with MOBSTER</summary>

**Output directory: `{publish_dir}/subclonal_deconvolution/pyclone_vi/<dataset>/<patient>/`**

- `<patient>_with_mobster_all_fits.h5`
  - HDF5 file for all possible fit
- `<patient>_with_mobster_best_fit.txt`
  - TSV file for the best fit
</details>

<details markdown="1">
<summary>Output files for all patients without MOBSTER</summary>

**Output directory: `{publish_dir}/subclonal_deconvolution/pyclone_vi/<dataset>/<patient>/`**

- `<patient>_all_fits.h5`
  - HDF5 file for all possible fit
- `<patient>_best_fit.txt`
  - TSV file for the best fit
  
</details>

#### VIBER

This folder contains the results of multivariate analysis using VIBER, which can be run prior to removal of mutations assigned to tail in all the samples.

<details markdown="1">
<summary>Output files for all patients with MOBSTER</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>/`**

- `<patient>_with_mobster_best_fit.rds`
  - RDS file for best fit
</details>


<details markdown="1">
<summary>Output files for all patients without MOBSTER</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>/`**

- `<patient>_best_fit.rds`
  - RDS file for best fit
  
</details>


## Clone Tree Inference

Subclonal deconvolution results are used to build clone tree from both single samples and multple samples using [ctree](https://caravagnalab.github.io/ctree/index.html). ctree is a R-based package which implements basic functions to create, manipulate and visualize clone trees by modelling Cancer Cell Fractions (CCF) clusters. Annotated driver genes must be provided in the input data.

> **NB:** When `-- tools pyclone-vi` is used, the output of PyClone-VI subclonal deconvolution is preprocessed prior to clone tree inference. Since ctree requires labeling one of the clusters as "clonal," the one with the highest CCF across all samples is choosen.

VIBER and MOBSTER fits are already suitable for ctree analysis.

### Single sample
<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{publish_dir}/subclonal_deconvolution/ctree/<patient>/<sample>/`**

- `<sample>_<tool>_ctree.rds`
  - RDS file for single sample clone tree

</details>

### Multi sample
<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{publish_dir}/subclonal_deconvolution/ctree/<patient>/`**

- `<patient>_<tool>_ctree.rds`
  - RDS file for single patient clone tree

</details>

## Unpublished results

### Formatter
The Formatter subworkflow is used to convert file to other formats and to standardize the output files resulting from different mutation (Mutect2, Strelka) and cna callers (ASCAT,Sequenza). Output files from this step are not published.

#### cnaparse

This parser aims at standardize into a unique format copy number calls and purity estimate from different callers.
<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{work_dir}/formatter/cna2CNAqc/<dataset>/<patient>/<sample>/`**
* `CNA.rds`
  * `rds` file containing parsed cna output in table format 
</details>

#### vcfparse

This parser aims at standardize into a unique format single nucleotide variants from different callers.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{work_dir}/formatter/vcf2CNAqc/<dataset>/<patient>/<sample>/`**
* `VCF.rds`
  * `rds` file containing parsed vcf in table format
</details>

#### tsvparse

This parser aims at converting mutations data of joint CNAqc analysis from CNAqc format (RDS file) into a tabular format (TSV file). This step is mandatory for running python-based tools (e.g. PyClone-VI, SigProfiler).

<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{work_dir}/formatter/CNAqc2tsv/<dataset>/<patient>/`**
* `joint_table.tsv`
  * TSV file containing cna and variants joint information.
</details>

### Lifter
The Lifter subworkflow is optional in multi-sample mode, when for a patient more samples are provided. The sub-workflow collect all mutations from the samples and perform pile-up of sample's private mutations in all the other samples.

#### get_positions
This intermediate step allows to retrieve private and shared mutations across samples originated from the same patient. 

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{work_dir}/lifter/mpileup/<dataset>/<patient>/<sample>/`**

* `positions_missing.bed`
  * BED file containing mutations to be retrieved for a given sample
</details>

<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{work_dir}/lifter/mpileup/<dataset>/<patient>/`**
* `all_positions.bed`
  * BED file containing mutation of all the samples of the same patient
</details>

#### mpileup

At this stage, [bcftools](https://samtools.github.io/bcftools/bcftools.html) is used to perform the pileup in order to retrieve frequency information of private mutations across all samples.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{work_dir}/lifter/mpileup/<dataset>/<patient>/<sample>/`**
- `pileup.vcf`
  - VCF file with called mutations
</details>


### Driver Annotation

_Add small description of this step_

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{work_dir}/Driver_Annotation/<dataset>/<patient>/<sample>/`**
* `driver_vcf.rds`
  * `rds` file containing variants with annotated drivers.
</details>