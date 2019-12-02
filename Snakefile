import pandas as pd
from snakemake.utils import min_version

### set min snakemake requirement
min_version("5.4.1")

### load config
configfile: "config.yaml"

### load sample data
samples          = pd.read_csv(config['general']['samples'], sep='\t').set_index(['sample', 'lane'], drop=False)
samples.index    = samples.index.set_levels([i.astype(str) for i in samples.index.levels])

### define commonly used variables
outdir = config['general']['paths']['outdir']
logdir = config['general']['paths']['logdir']

ALIGN_GENOME     = expand(os.path.join(outdir, 'alignments', '{sample}.all.bam.bai'), sample=samples['sample'].unique())
ALIGN_GENOME_CRM = expand(os.path.join(outdir, 'alignments', '{sample}.all.crai'), sample=samples['sample'].unique())
ALIGN_TRANSCRIPT = expand(os.path.join(outdir, 'alignments_transcript', '{sample}.all.transcript.sort.bam.bai'), sample=samples['sample'].unique()) 
ALIGN_RSEM       = expand(os.path.join(outdir, 'expression', 'rsem', '{sample}.all.transcript.sort.rsem.bam'), sample=samples['sample'].unique()) 
DE_GENES_DESEQ   = os.path.join(outdir, 'diffexp','deseq','DE_deseq.tsv')
DE_GENES_LIMMA   = os.path.join(outdir, 'diffexp','limma_voom', 'DE_limma.tsv')
EXPR_RSEM        = expand(os.path.join(outdir, 'expression', 'rsem', '{sample}.genes.results'), sample=samples['sample'].unique()) 
EXPR_SIMPLE      = os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.tsv')
EXPR_SIMPLE_NA   = os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.tsv')
EXPR_HTSEQ       = expand(os.path.join(outdir, 'expression', 'htseq', '{sample}.htseq_counts.txt'), sample=samples['sample'].unique())
QC_FASTQC        = expand(os.path.join(outdir, 'qc', 'fastqc', '{i.sample}_{i.lane}_fastqc.done'), i=samples.itertuples())
QC_MULTIQC       = os.path.join(outdir, 'qc', 'multiqc', 'multiqc_report.html')
QC_RNASEQC        = expand(directory(os.path.join(config['general']['paths']['outdir'], 'qc', 'rnaseqc', '{sample}')),sample=samples['sample'].unique())
QC_PICARD         = expand(os.path.join(outdir, 'qc', 'picard', '{sample}.RNAMetrics.txt'), sample=samples['sample'].unique())
QC_RSEQC_BAMSTAT  = expand(os.path.join(outdir, 'qc', 'rseqc', '{sample}.rseqc.bam_stat.txt'), sample=samples['sample'].unique())
QC_RSEQC_READGC   = expand(os.path.join(outdir, 'qc', 'rseqc', '{sample}.rseqc.read_gc.GC_plot.pdf'), sample=samples['sample'].unique()) 
RSEQC_PLOT        = expand(os.path.join(outdir, 'qc', 'rseqc', 'gc_content.pdf'))
QC_QORTS          = expand(os.path.join(outdir, 'qc', 'qorts', '{sample}.qorts', '{sample}_QC.multiPlot.png'), sample=samples['sample'].unique())
QORTS_PLOTS       = os.path.join(outdir, 'qc', 'qorts', 'qorts.genebody.pdf')
TCGA_BOXPLOT     = expand(os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.boxplot.done'), sample=samples['sample'].unique()) 
SPLICE_EVENTS_S  = expand(os.path.join(outdir, 'splicing', 'spladder', 'genes_graph_conf' + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.pickle'), sample=samples['sample'].unique())
SPLICE_COUNTS_S  = expand(os.path.join(outdir, 'splicing', 'spladder', 'genes_graph_conf' + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.hdf5'), sample=samples['sample'].unique())
SPLICE_COUNTS_M  = os.path.join(outdir, 'splicing', 'spladder', 'genes_graph_conf' + str(config['tools']['spladder']['confidence'])  + '.merge_graphs.count.hdf5')
SPLICE_BURDEN    = expand(os.path.join(outdir, 'splicing', 'spladder', 'genes_graph_conf' + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.G0.01.globsum{globsum}.conf3_neojunctions.pdf'), sample=samples['sample'].unique(), globsum=config['tools']['burden']['globsum'])
HYPO_SCORE       = expand(os.path.join(outdir, 'hypoxia', '{sample}.hypoxia_score.tsv'), sample=samples['sample'].unique())
HYPO_PLOT_ALL    = expand(os.path.join(outdir, 'hypoxia', '{sample}.hypoxia_allcancer_boxplot.pdf'), sample=samples['sample'].unique())
HYPO_PLOT_SKIN   = expand(os.path.join(outdir, 'hypoxia', '{sample}.hypoxia_skin_boxplot.pdf'), sample=samples['sample'].unique())

### define rules not to be executed on the cluster
localrules: all

### define target rules
rule all:
    input:
        ALIGN_GENOME,
        ALIGN_GENOME_CRM,
        ALIGN_TRANSCRIPT,
        ALIGN_RSEM,
        QC_FASTQC,
        QC_MULTIQC,
        QC_PICARD,
        QC_RSEQC_BAMSTAT,
        QC_RSEQC_READGC,
        RSEQC_PLOT,
        QC_QORTS,
        QORTS_PLOTS,
        EXPR_SIMPLE,
        EXPR_SIMPLE_NA,
        EXPR_HTSEQ,
        EXPR_RSEM,
        TCGA_BOXPLOT,
        HYPO_SCORE,
        HYPO_PLOT_ALL,
        HYPO_PLOT_SKIN,
        SPLICE_EVENTS_S,
        SPLICE_COUNTS_S,
        SPLICE_COUNTS_M,
        SPLICE_BURDEN
      
### include relevant rules
include: "rules/common.smk"
include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/splicing.smk"
include: "rules/expression.smk"
include: "rules/hypo_score.smk"
include: "rules/diffexp.smk"
include: "rules/tcga_boxplot.smk"
