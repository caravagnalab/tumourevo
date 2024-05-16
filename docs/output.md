# nf-core/evoverse: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

<!-- ## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info) -->


## Variant annotation
This step perform annotation of variants and mutations using genomic data.

### VEP
[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) (Variant Effect Predictor) is a `Ensembl` tool that determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence.

This step should be started from `vcf` files.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/variant_annotation/VEP/<dataset>/<patient>/<sample>/`**

- `<sample>.vep.summary.html`
  - Summary of the VEP run to be visualised with a web browser
- `<sample>.vep.vcf.gz`
  - annotated vcf file

</details>


### vcf2maf

[vcf2maf](https://github.com/mskcc/vcf2maf) convert a VCF file into a MAF (Mutation annotation Format), where each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. vcf2maf is designed to work with VEP. 

**NB: While VEP is tolerant of chromosome format mismatches (when the input .vcf file uses the UCSC format chrN and the reference fasta uses Ensembl/NCBI format N), vcf2maf is not. Make sure the reference fasta chromosome format matches that of your input.**

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/vcf2maf/<dataset>/<patient>/<sample>/`**
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

**Output directory: `{outdir}/maftools/<dataset>/`**
- `<dataset>.maftools.rds`
  - summarized MAF object
- `<dataset>.maftools.pdf`
  - summary plots

</details>


## Signature Deconvolution

Mutational signatures are characteristic patterns of somatic mutations in cancer genomes, reflecting the underlying mutational processes. This step performs de novo extraction, inference, deciphering or deconvolution of mutational counts.

The available tools for this step are:
- SparseSignatures (Bioconductor R package)
- SigProfiler (Python framework)
- CNAqc (R package)


### SparseSignatures

[SparseSignatures](https://www.bioconductor.org/packages/release/bioc/html/SparseSignatures.html) R-based computational framework that provides a set of functions to extract and visualize the mutational signatures that best explain the mutation counts of a large number of patients. In particular:
- reliably extracts mutational signatures and quantifies their activity;
- incorporates an explicit background model to improve the inference;
- exploits LASSO regularization to reduce the impact of overfitting;
- implements bi-cross-validation to select the best number of signatures

`mCNAqc` object is used as input for this step.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/results/SparseSignatures/<dataset>/`**
- `<dataset>.SparseSig.rds`
  - signatures best configiration object
- `<dataset>.SparseSig.pdf`
  - signatures plot

</details>
  

### SigProfiler

[SigProfiler](https://osf.io/t6j7u/wiki/home/) is a python framework that allows de novo extraction of mutational signatures from data generated in a matrix format. The tool identifies the number of operative mutational signatures, their activities in each sample, and the probability for each signature to cause a specific mutation type in a cancer sample. The tool makes use of `SigProfilerMatrixGenerator` and `SigProfilerPlotting`, seamlessly integrating with other `SigProfiler` tools.

`mCNAqc` object is used as input for this step.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/results/SigProfiler/<dataset>/results/`**

</details>


## Subclonal deconvolution

In the context of tumor evolution, subclonal reconstruction consists in the identification of cancer clones by leveraging  variant read counts and the associated variant allele frequency (VAF) of somatic mutations, adjusted for copy-number status and tumor purity. 

The results of subclonal decovnultion step are collected in `{outdir}/subclonal_deconvolution/` directory.

`mCNAqc` object is used as input for this step.

### Single sample

For single sample mode, clones clustering and neutral tail mutations are detected seperately for each sample.

#### MOBSTER

[MOBSTER](https://caravagnalab.github.io/mobster/) processes mutant allelic frequencies to identify and remove neutral tails from the input data, so that subclonal reconstruction algorithms can be applied downstream to find subclones from the processed read counts.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/mobster/<dataset>/<patient>/<sample>/`**

- `<sample>_fit.rds`
  - `.rds` object contains the fit of subclonal deconvolution
- `<sample>.pdf`

</details>

#### PyClone-VI

[PyClone-VI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03919-2) is a computationally efficient Bayesian statistical method for inferring the clonal population structure of cancers, by considering allele fractions and coincident copy number variation using a variational inference approach. 


<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/pyclone_vi/<dataset>/<patient>/<sample>/`**

- `<sample>_all_fits.h5`
  - HDF5 file for all fit
- `<sample>_best_fit.txt`
  - tsv file for the best fit
  
</details>

#### VIBER

[VIBER](https://caravagnalab.github.io/VIBER/index.html) is a package that implements a variational Bayesian model to fit multi-variate Binomial mixtures. In the context of subclonal deconvolution, VIBER models read counts. 

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>/<sample>/`**

- `<sample>_best_st_fit.rds`
  - rds file for best fit

</details>

### Multi samples

When multiple samples for the same patient are available (e.g. multi-regional or longitudinal experiment), we can be intersted in visualizing the tumor evolution at different level. In multi sample mode, a patient-specific subclonal deconvolution is perfomed. In this step, the user has the possibility to use the result of single sample subclonal deconvolution using `MOBSTER` to remove mutations assigned to the tail before proceeding with multi-variate clone identification. This option is controlled by the `--tail` parameter.


#### Pyclone-VI

##### Without MOBSTER

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/pyclone_vi/<dataset>/<patient>/`**

- `<patient>_all_fits.h5`
  - HDF5 file for all fit
- `<patient>_best_fit.txt`
  - tsv file for the best fit
  
</details>

##### With MOBSTER

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/pyclone_vi/<dataset>/<patient>/`**

- `<patient>_with_mobster_all_fits.h5`
  - HDF5 file for all fit
- `<patient>_with_mobster_best_fit.txt`
  - tsv file for the best fit

</details>

#### VIBER


##### Without MOBSTER

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>/`**

- `<patient>_best_fit.rds`
  - rds file for best fit

</details>

##### With MOBSTER

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>/`**

- `<patient>_with_mobster_best_fit.rds`
  - rds file for best fit
  
</details>

<!-- ## Clone Tree Inference

### Single sample
<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/ctree/<patient>/<sample>/`**

- `<sample>_ctree.rds`
  - rds file for best fit

</details>

### Multi sample
<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/ctree/<patient>/`**

- `<patient>_ctree.rds`
  - rds file for best fit

</details> -->
