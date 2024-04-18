# nf-core/evoverse: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images


## Variant annotation
This step perform annotation of variants and mutations using genomic data.

### VEP
[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) (Variant Effect Predictor) is a `Ensembl` tool that determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence.

This step should be started from `vcf` files.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/results/VEP/{dataset,patient,sample}/`**

- `<sample>.VEP.summary.html`
  - Summary of the VEP run to be visualised with a web browser
- `<sample>.vep.vcf.gz`
  - annotated vcf file
  </details>


### vcf2maf

[vcf2maf](https://github.com/mskcc/vcf2maf) convert a VCF file into a MAF (Mutation annotation Format), where each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. vcf2maf is designed to work with VEP. 

**NB: While VEP is tolerant of chromosome format mismatches (when the input .vcf file uses the UCSC format chrN and the reference fasta uses Ensembl/NCBI format N), vcf2maf is not. Make sure the reference fasta chromosome format matches that of your input.**

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/results/vcf2maf/{dataset,patient,sample}/`**
- `<sample>.vcf2maf.maf`
  - annotated maf file
  </details>


### maftools

[maftools]( https://bioconductor.org/packages/release/bioc/html/maftools.html) R-based tool that provides a comprehensive set of functions for processing MAF files and to perform most commonly used analyses in cancer genomics. In particular, maftools summarize, analyze and visualize MAF (Mutation Annotation Format) files. 
It requires somatic variants in a MAF file which must be gz compressed.

MAF fields requirements:

- Mandatory fields: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type` and `Tumor_Sample_Barcode`.

- Optional fields: `VAF` (Variant Allele Frequency), `amino acid change` information.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/results/maftools/{dataset}/`**
- `<dataset>.maftools.rds`
  - summarized MAF object
- `<dataset>.maftools.pdf`
  - summary plots
  </details>


## Signature Deconvolution

Mutational signatures are characteristic patterns of somatic mutations in cancer genomes, reflecting the underlying mutational processes. This step performs de novo extraction, inference, deciphering or deconvolution of mutational counts.

The available tools for this step are:
- SparseSignatures (Bioconductor R package)
- SigProfilerMatrixGenerator (Python framework)
- SigProfilerExtractor (Python framework)
- SigProfilerPlotting
- CNAqc (R package)

### SparseSignatures

[SparseSignatures](https://www.bioconductor.org/packages/release/bioc/html/SparseSignatures.html) R-based computational framework that provides a set of functions to extract and visualize the mutational signatures that best explain the mutation counts of a large number of patients. In particular:
- reliably extracts mutational signatures and quantifies their activity;
- incorporates an explicit background model to improve the inference;
- exploits LASSO regularization to reduce the impact of overfitting;
- implements bi-cross-validation to select the best number of signatures

The following parameters can be tuned for this step:

- `K` - the candidate numbers of signatures (min.value = 2) to be fit to the dataset;
- `lambda_values_beta` - the range of values of the signature sparsity parameter;
- `cross_validation_repetitions` - the number of repetitions of the cross-validation procedure.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/results/SparseSignatures/{dataset}/`**
- `<dataset>.SparseSig.rds`
  - signatures best configiration object
- `<dataset>.SparseSig.pdf`
  - signatures plot
  </details>

### SigProfiler

[SigProfiler](https://osf.io/t6j7u/wiki/home/) is a python framework that allows de novo extraction of mutational signatures from data generated in a matrix format. The tool identifies the number of operative mutational signatures, their activities in each sample, and the probability for each signature to cause a specific mutation type in a cancer sample. The tool makes use of `SigProfilerMatrixGenerator` and `SigProfilerPlotting`, seamlessly integrating with other `SigProfiler` tools.

The following parameters can be tuned for this step:

- `minimum_signatures` - the minimum number of signatures to be extracted (default = 1); 
- `maximum_signatures` - the maximum number of signatures to be extracted (default = 25). 

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/results/SigProfiler/results/{dataset}/`**


## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)



