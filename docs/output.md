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
