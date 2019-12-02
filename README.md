# RNA-Seq based pipeline

This repository contains the backbone RNA-Seq pipeline following the Snakemake gold standard
implementation. Currently, the following steps are implemented:

* QC with fastqc, multiQC, and Picard tool
* alignment against the genome with STAR
* expression quantification using Htseq2 and RSEM at both gene and transcript levels
* construction of splicing graphs per input sample using SplAdder
* compute the hypoxia score using R package sspath
* differential expression analysis using DeSeq2 and LIMMA-Voom
* generate TCGA plots for HIF genes using R-scripts


## Conda environment

The tools used in the pipeline are defined in a conda environment file in the `envs` directory.
This environment can be built once before the pipeline is run on many samples and then re-used for
any further run. However, if the definition of the environment (`envs/wf_rna_seq.yml`) changes, the
environment will be recreated.

## Preliminary steps

### Conda

To generate the conda environment upfront, one can run the following on a server with internet
access:

```bash
snakemake --cores 1 -p --use-conda --conda-prefix $CONDA_DIR --create-envs-only
```
Here, `$CONDA_DIR` should contain (or be replaced with) a parent directory for the respective conda
environment to be created.

### Samples

To specify any samples to be run on or specific parameters to be changed `samples.tsv` and
`config.yaml` need to be adapted.

## Executing the pipeline

To use the `workflow-rna-seq-star` pipeline on a single compute node, you will need to follow the command:

```
bash
snakemake --cores 10 -p --use-conda --conda-prefix $CONDA_DIR
```
Again, `$CONDA_DIR` contains the (previously used) path to the conda environments.

## Examples of outputs from each of the rules

Successful execution of the snakemake workflow will generate multiple files in different formats for each input sample. In addition to that, each rule generates a set of output files in different directories. The following table briefly defines the outputs files from each step.

| Rules | Description | Directory | Tool |
| ---------- | --- | --------------------------------------------------------------------------------------- |------ |
| QC_FASTQC | generates quality control reports of the raw fastq files and returns `html`and `zip` files  | qc/fastqc | fastqc =0.11.8 |
| QC_PICARD | runs Picard CollectRnaSeqMetrics on the BAm files | qc/picard | picard =2.18.4|
| QC_QORTS | runs QoRTs and generates geneBodyCoverage plot and QoRTs multiplot for each sample | qc/qorts | qorts =1.3.0 | 
| QC_RSEQC_BAMSTAT | runs RSeQC tool to generate BAMSTAT on BAMs | qc/rseqc | rseqc =3.0.0|
| QC_RSEQC_READGC | runs RSeQC tool to generate GC content | qc/rseqc | rseqc =3.0.0 |
| QC_MULTIQC | aggregates results from tools like FASTQC, RSEM, Picard, HTSeq to a html file  | qc/multiqc | multiqc =1.7|
| RSEQC_PLOT | generates multisample RSeQC GC content plot | qc/rseqc | rseqc =3.0.0 |
| QORTS_PLOT | generates multisample QoRTs geneBodyCoverage plot | qc/qorts | qorts =1.3.0 & R QoRTs library | 
| ALIGN_GENOME| aligns the raw fastq files, returns aligned files in `*.BAM` format | alignments | star =2.7.0 |
| ALIGN_TRANSCRIPT| aligns the raw fastq files returns aligned files at transcript level in `.bam` format | alignments_transcript | star =2.7.0 |
| ALIGN_RSEM EXPR_RSEM| returns the aligned `.bam` file for each input genes, and quantification results of those gene on isoform and gene level as `.genes.results` and `.isoforms.results` files | expression/rsem| rsem =1.3.1|
| EXPR_HTSEQ| returns gene-level quantified files in `*.txt` format | expression/htseq | htseq =0.9.1|
| EXPR_SIMPLE EXPR_SIMPLE_NA| returns expression quantified files in `tsv` format for each input sample. NA stands for expression quantification for non-alternative, reutrns `non_alt.tsv` file, and a combined `hdf5` file for all samples| expression/simple_counting | count_hdf2tsv |
| DE_GENES_DESEQ DE_GENES_LIMMA| returns `*.tsv` files each from two tools with differential expression values for each genes with estimated log2 fold changes , p-values, and adjusted p-values (FDR)  | diffexp/deseq diffexp/limma_voom | deseq2 =1.22.1 limma =3.28.2|
| HYPO_SCORE | returns hypoxia score for selected HIF genes (n=30) as a `tsv` file| hypoxia | hypoxia | r-sspaths =0.1.1 |
| HYPO_PLOT_ALL HYPO_PLOT_SKIN | generates hypoxia boxplots for a sample with all TCGA cancer types, GTEX samples and boxplot with TCGA Melanoma and GTEX skin | hypoxia | r-sspaths =0.1.1 & rscripts |
| TCGA_BOXPLOT | generates boxplots for patient samples in the context of TCGA for the genes in [`this file`](https://github.com/ratschlab/workflows-rna_seq_star/blob/master/resources/tcga_boxplot/genes_of_interest.txt). To plot more genes, add them to the above file in the format `HUGO_symbol,Ensembl_ID`| tcga_boxplot/png_images | r-scripts|
| SPLICE_EVENTS_S SPLICE_COUNTS_S SPLICE_COUNTS_M SPLICE_BURDEN  | returns graphs of `BAM` files in `pickle` format, quantified graphs in `.count.hdf5` formats, merged graphs and splicing burden in same formats| splicing/spladder | spladder==2.2.0 |

## Customization of the workflow for different projects

The current workflow can be composed and customized according to your requirements by modifying the configuration settings and the snakefile. For the sake of the documentation, let us take the example of the in-house Fly-Base (FB) project. The output files for FB project are hypoxia score, splicing burden and TCGA plots. To generate only those outputs you will need to modify and comment some of the global variables and rules within some files. By doing so, you will see the execution of other rules will be skipped, and it results in only the required outputs. The following are the files you will need to modify.

* **sample**: Add the new TP samples IDs, lane information and fastq file names to [this tsv file](https://github.com/ratschlab/workflows-rna_seq_star/blob/master/samples.tsv) in the corresponding columns.

* **snakefile**: Comment the following target rules, in the [`snakefile`](https://github.com/ratschlab/workflows-rna_seq_star/blob/master/Snakefile).

    - `DE_GENES_DESEQ`: compute differential expression analysis using deseq2
    - `DE_GENES_LIMMA`: compute differential expression analysis using LIMMA -voom
    - `EXPR_HTSEQ`    : quantify the expression of genes using HTSeq2
    -  `include: "rules/diffexp.smk"` : Rules for differential expression analysis using deseq2 and LIMMA-voom
    

* **config**: Within the [configuration file](https://github.com/ratschlab/workflows-rna_seq_star/blob/master/config.yaml) file please change the paths for the following entries

    - `fastqdir`: Replace with the absolute path where your `fastq` files are located.
    - `outdir`: Replace with the path where you want your `results`
    - `logdir`: Replace with the path where you want your `log files`
    - `workdir`: Replace with the path where you want your `workdir`
    - `samples_list`: The list of samples you plan to run

Leave all other entries in above-mentioned files as default. By this, you are now ready to execute the workflow for the TP project. To run the workflow in the cluster please follow the following instructions.

## Executing the pipeline on the cluster

To run the pipeline using the LSF scheduler, the following command can be run:

```bash
snakemake --use-conda --conda-prefix $CONDA_DIR --configfile config.yaml --cluster 'bsub -M {params.res_mem} -W {params.res_time} -n {threads} -R "rusage[mem={params.res_mem}]" -o {params.lsf_log}' -j 20 -k -p --latency-wait 120
```
Specifically, the following items are:
* `$CONDA_DIR` contains the path to the conda environments
* `config.yaml` contains the specific run-configurations for the project
* `{params.*}` and `{threads}` contain the tool-specific resources defined in the rules / config
* `-j 20` allows a max of 20 jobs to be run at the same time
* `-k` will continue the workflow even if a single jkob crashed but independent jobs can continue
* `-latency-wait 120` increases the wait-time for result files to 2 minutes (useful in a distributed system)

## Snakemake Reports

Snakemake reports `report.html` can be generated once the pipeline run is over.

The current report have the hypoxia plots, qorts geneBody Coverage plot, RSeQC GC content plot and MultiQC html file.

To generate the snakemake report
```bash
snakemake --report report.html
```

To manage the order of tools in MultiQC html file, please edit [`this file`](https://github.com/ratschlab/workflows-rna_seq_star/blob/master/multiqc_config.html)

### ToDo

- scripts/tcga_boxplot/*, move to central gromics repo
- Make the rule to generate rsem_reference
