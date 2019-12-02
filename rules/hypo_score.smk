"""
This rule generates hypoxia scores for each sample
"""
rule hypoxia_score:
    input:
        count=os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.tsv')
    output:
        fname=os.path.join(outdir, 'hypoxia', '{sample}.hypoxia_score.tsv'),
    params:
        result_dir=os.path.join(outdir, 'hypoxia'),
        tcga_expression=config['general']['paths']['tcga_expression'],
        cancertype=config['general']['cancertype'],
        res_mem=config['tools']['hypoxia']['mem'],
        res_time=config['tools']['hypoxia']['time'],
        lsf_log=os.path.join(logdir, 'hypoxia', '{sample}.hypoxia.lsf.log')
    conda:
        "../envs/wf_rna_seq_R_hypoxia.yml"
    threads: 
        config['tools']['hypoxia']['threads']
    log:
       os.path.join(logdir, 'hypoxia', '{sample}.log')
    benchmark:
       os.path.join(logdir, 'hypoxia', '{sample}.benchmark.log')
    shell:
        """
        tupro_pipeline_hypoxia.R {params.cancertype} {params.tcga_expression} {input.count} {output.fname}
        """

"""
This rule makes hypoxia boxplots for a sample with all TCGA cancer types, GTEX samples
Also another boxplot with TCGA Melanoma and GTEX skin
"""
rule hypoxia_boxplots:
    input:
        patient_hypoScoreFile=os.path.join(outdir, 'hypoxia', '{sample}.hypoxia_score.tsv')
    params:
        res_mem=config['tools']['hypoxia_plot']['mem'],
        res_time=config['tools']['hypoxia_plot']['time'],
        sampleName='{sample}',
        hypoxia_plot=config['tools']['hypoxia_plot']['call'],
        tcga_gtex_hypo_scores=config['general']['paths']['tcga_gtex_hypo_scores'],
        lsf_log=os.path.join(logdir, 'hypoxia', '{sample}.hypoxia_boxplot.lsf.log'),
        outdir=os.path.join(outdir, 'hypoxia')
    output:
        out_fig1=report(os.path.join(outdir, 'hypoxia', '{sample}.hypoxia_allcancer_boxplot.pdf'), caption="../report/hypoxia_allcancer.rst", category="Hypoxia Boxplots"),
        out_fig2=report(os.path.join(outdir, 'hypoxia', '{sample}.hypoxia_skin_boxplot.pdf'), caption="../report/hypoxia.rst", category="Hypoxia Boxplots")
    log:
        os.path.join(logdir, 'hypoxia', '{sample}.hypoxia_boxplot.log')
    benchmark:
        os.path.join(logdir, 'hypoxia', '{sample}.hypoxia_boxplot.benchmark.log')
    conda:
        "../envs/wf_rna_seq_R_hypoxia.yml"  
    threads:
        config['tools']['hypoxia_plot']['threads']
    shell:
        'Rscript --vanilla {params.hypoxia_plot} {params.sampleName} {params.tcga_gtex_hypo_scores} {input.patient_hypoScoreFile} {params.outdir} {output.out_fig1} {output.out_fig2} ' 
