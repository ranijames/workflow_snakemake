"""
FastQC
"""
rule fastqc:
    input:
        fastq=get_fastqs
    params:
        outdir=os.path.join(outdir, 'qc', 'fastqc'),
        res_mem=config['tools']['fastqc']['mem'],
        res_time=config['tools']['fastqc']['time'],
        lsf_log=os.path.join(logdir, 'qc', '{sample}_{lane}_fastqc.lsf.log')
    output:
        outfile=os.path.join(outdir, 'qc', 'fastqc', '{sample}_{lane}_fastqc.done')
    log:
        os.path.join(logdir, 'qc', '{sample}_{lane}_fastqc.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['fastqc']['threads']
    shell:
        'mkdir -p {params.outdir}; fastqc -o {params.outdir} -t {threads} {input.fastq} > {log} 2>&1 && touch {output.outfile}'

"""
Multi QC
"""
rule multiqc:
    input:
        expand(os.path.join(outdir, 'qc', 'fastqc', '{i.sample}_{i.lane}_fastqc.done'), i=samples.itertuples()),
        expand(os.path.join(outdir, 'expression', 'htseq', '{sample}.htseq_counts.txt'),sample=samples['sample'].unique()),
        expand(os.path.join(outdir, 'qc', 'picard', '{sample}.RNAMetrics.txt'),sample=samples['sample'].unique()),
        expand(os.path.join(outdir, 'expression', 'rsem', '{sample}.genes.results'),sample=samples['sample'].unique()),
        expand(os.path.join(outdir, 'qc', 'qorts', '{sample}.qorts', '{sample}_QC.multiPlot.png'), sample=samples['sample'].unique()),
        expand(os.path.join(outdir, 'qc', 'rseqc', '{sample}.rseqc.bam_stat.txt'), sample=samples['sample'].unique()), 
        expand(os.path.join(outdir, 'qc', 'rseqc', '{sample}.rseqc.read_gc.GC_plot.pdf'), sample=samples['sample'].unique()) 
    params:
        fastqc_dir=os.path.join(outdir, 'qc', 'fastqc'),
        htseq_dir=os.path.join(outdir, 'expression', 'htseq'),
        picard_dir=os.path.join(outdir, 'qc', 'picard'),
        rsem_dir=os.path.join(outdir, 'expression', 'rsem'),
        rseqc_dir=os.path.join(outdir, 'qc', 'rseqc'),
        config=config['tools']['multiqc']['config'],
        out_dir=os.path.join(outdir, 'qc', 'multiqc'),
        res_mem=config['tools']['multiqc']['mem'],
        res_time=config['tools']['multiqc']['time'],
        lsf_log=os.path.join(logdir, 'qc', 'multiqc.lsf.log')
    output:
        multiqcoutfile=report(os.path.join(outdir, 'qc', 'multiqc', 'multiqc_report.html'), caption="../report/multiqc.rst", category="Quality Control")
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        'mkdir -p {params.out_dir}; multiqc -d {params.fastqc_dir} -c {params.config} {params.htseq_dir} {params.picard_dir} {params.rsem_dir} {params.rseqc_dir} -f -o {params.out_dir}'

"""
rseqc bamstat
"""
rule rseqc_bam_stat:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
        bai=os.path.join(outdir, 'alignments', '{sample}.all.bam.bai')
    params:
        res_mem=config['tools']['rseqc']['bamstat']['mem'],
        res_time=config['tools']['rseqc']['bamstat']['time'],
        lsf_log=os.path.join(logdir, 'qc', 'rseqc', '{sample}_rseqc_bamstat.lsf.log')
    output:
        outfile=os.path.join(outdir, 'qc', 'rseqc', '{sample}.rseqc.bam_stat.txt')
    log:
        os.path.join(logdir, 'qc', 'rseqc', '{sample}_rseqc_bamstat.log')
    benchmark:
        os.path.join(logdir, 'qc', 'rseqc', '{sample}_rseqc_bamstat.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['rseqc']['bamstat']['threads']
    shell:
        ('{config[tools][rseqc][bamstat][call]} -i {input.bam} > {output.outfile}')

"""
rseqc readGC
"""
rule rseqc_readGC:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
        bai=os.path.join(outdir, 'alignments', '{sample}.all.bam.bai')
    params:
        res_mem=config['tools']['rseqc']['readGC']['mem'],
        res_time=config['tools']['rseqc']['readGC']['time'],
        lsf_log=os.path.join(logdir, 'qc', 'rseqc', '{sample}_rseqc_readGC.lsf.log'),
        txt=os.path.join(outdir, 'qc', 'rseqc','{sample}.rseqc.read_gc'),
        #donefile=os.path.join(outdir, 'qc', 'rseqc', 'rseqc_readGC.done')
    output:
        outfile=os.path.join(outdir, 'qc', 'rseqc', '{sample}.rseqc.read_gc.GC_plot.pdf')
    log:
        os.path.join(logdir, 'qc', 'rseqc', '{sample}_rseqc_readGC.log')
    benchmark:
        os.path.join(logdir, 'qc', 'rseqc', '{sample}_rseqc_readGC.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['rseqc']['readGC']['threads']
    shell:
        '{config[tools][rseqc][readGC][call]} -i {input.bam} -o {params.txt} '

"""
QoRTs
"""

rule qorts:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.readname_sort.all.bam')
    output:
        out=os.path.join(outdir, 'qc', 'qorts', '{sample}.qorts', '{sample}_QC.multiPlot.png')
    params:
        res_mem=config['tools']['qorts']['mem'],
        res_time=config['tools']['qorts']['time'],
        lsf_log=os.path.join(logdir, 'qc', 'qorts', '{sample}.qorts.lsf.log'),
        gtf = config['general']['paths']['annotation'],
        outdir = os.path.join(outdir, 'qc', 'qorts', '{sample}.qorts'),
        name = '{sample}'
    log:
        os.path.join(logdir, 'qc', 'qorts', '{sample}.qorts.log')
    benchmark:
        os.path.join(logdir, 'qc', 'qorts', '{sample}.qorts.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['qorts']['threads']
    shell:
        ('mkdir -p {params.outdir}; samtools view -b -F4 -F8 {input.bam} | {config[tools][qorts][call]} QC --runFunctions writeGeneBody,writeGenewiseGeneBody --generatePdfReport --generatePlots - {params.gtf} {params.outdir}; cp {params.outdir}/QC.multiPlot.png {params.outdir}/{params.name}_QC.multiPlot.png')


"""
Qorts geneBodyCoverage plot
"""
rule qorts_plot:
    input:
        expand(os.path.join(outdir, 'qc', 'qorts', '{sample}.qorts', '{sample}_QC.multiPlot.png'), sample=samples['sample'].unique())
    params:
        out_dir=os.path.join(outdir, 'qc', 'qorts'),
        qorts = config['tools']['qorts']['plot'],
        res_mem=config['tools']['plots_stats']['mem'],
        res_time=config['tools']['plots_stats']['time'],
        lsf_log=os.path.join(logdir, 'qc', 'qorts', 'qorts_plot.lsf.log'),
        htseq_dir=os.path.join(outdir, 'quantification', 'htseq'),
        qorts_dir=os.path.join(outdir, 'qc', 'qorts'),
        samples_list=config['general']['paths']['samples_list'],
    output:
        out=report(os.path.join(outdir, 'qc', 'qorts', 'qorts.genebody.pdf'), caption="../report/qorts.rst", category="Qorts Plots")
    log:
        os.path.join(logdir, 'qc', 'qorts', 'qorts_MultiPlot.log')
    benchmark:
        os.path.join(logdir, 'qc', 'qorts', 'qorts_MultiPlot.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['plots_stats']['threads']
    shell:
        'echo -e "unique.ID\tsample.ID\tgroup.ID\tqc.data.dir" >{params.out_dir}/decoder.txt;'
        'while read samp; do echo -e "$samp.fastq.gz\t$samp.fastq.gz\t$samp\t{params.out_dir}/$samp.qorts">>{params.out_dir}/decoder.txt; done <{params.samples_list};'
        'Rscript {params.qorts} {params.out_dir}/decoder.txt {output.out}'

"""
RSeQC GC content plot
"""
rule rseqc_plot:
    input:
        expand(os.path.join(outdir, 'qc', 'rseqc', '{sample}.rseqc.read_gc.GC_plot.pdf'), sample=samples['sample'].unique())
    params:
        out_dir=os.path.join(outdir, 'qc', 'rseqc'),
        plot_gc=config['tools']['plots_stats']['gc_plot'],
        plot_genebodycoverage=config['tools']['plots_stats']['genebody_coverage'],
        generate_readStats = config['tools']['plots_stats']['generate_readStats'],
        qorts = config['tools']['qorts']['plot'],
        res_mem=config['tools']['plots_stats']['mem'],
        res_time=config['tools']['plots_stats']['time'],
        lsf_log=os.path.join(logdir, 'qc', 'rseqc', 'rseqc_plot.lsf.log'),
        htseq_dir=os.path.join(outdir, 'quantification', 'htseq'),
        qorts_dir=os.path.join(outdir, 'qc', 'qorts'),
    output:
        rseqc_gcplot=report(os.path.join(outdir, 'qc', 'rseqc', 'gc_content.pdf'), caption="../report/rseqc_multisample.rst", category="RSeQC Plots")
    log:
        os.path.join(logdir, 'qc', 'rseqc', 'rseqc_plot.log')
    benchmark:
        os.path.join(logdir, 'qc', 'rseqc', 'rseqc_plot.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['plots_stats']['threads']
    shell:
       'R --vanilla --args {params.out_dir} {output.rseqc_gcplot} < {params.plot_gc}' 


"""
rule rnaseqc:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}_all.Aligned.sortedByCoord.out.bam')
    params:
        out_dir=os.path.join(outdir, 'qc', 'rnaseqc')
    output:
        outfile=directory(os.path.join(outdir, 'qc', 'rnaseqc', '{sample}'))
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'mkdir -p {params.out_dir}; rnaseqc {config[general][paths][annotation_collapsed]} {input.bam} {output} && touch {output.outfile}'
"""

"""
Picard CollectRnaSeqMetrics -  produces metrics describing the distribution of the bases within the transcripts
"""
rule picardCollectRnaSeqMetrics:
    input:
         bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
         bai=os.path.join(outdir, 'alignments', '{sample}.all.bam.bai')
    params:
        ref_flat=os.path.join(config['general']['paths']['ref_flat']),
        strand=os.path.join(config['general']['picard_strand']),
        res_mem=config['tools']['picard']['CollectRnaSeqMetrics']['mem'],
        res_time=config['tools']['picard']['CollectRnaSeqMetrics']['time'],
        lsf_log=os.path.join(logdir, 'qc', 'picard', '{sample}_picard.lsf.log')
    output:
        outfile=os.path.join(outdir, 'qc', 'picard', '{sample}.RNAMetrics.txt'),
    log:
        os.path.join(logdir, 'qc', 'picard', '{sample}_picard.log')
    benchmark:
        os.path.join(logdir, 'qc', 'picard', '{sample}_picard.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['picard']['CollectRnaSeqMetrics']['threads']
    shell:
        ('picard -Xmx$(({params.res_mem} * {threads}))m CollectRnaSeqMetrics ' +
        'INPUT={input.bam} ' +
        'REF_FLAT={params.ref_flat} ' +
        'OUTPUT={output.outfile} ' +
        'VALIDATION_STRINGENCY=LENIENT ' + 
        'STRAND_SPECIFICITY={params.strand} > {log} 2>&1 ')
