def deseq2_threads(wildcards=None):
    # suggested by the author of DESEQ2 #MikeLove
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(samples) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule limma:
    input:
        counts=os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.tsv'),
        conditions=config['general']['conditions']
    output:
        diff_counts=os.path.join(outdir, 'diffexp','limma_voom', 'DE_limma.tsv')
    threads: 
        deseq2_threads
    params:
        res_mem=config['tools']['diff']['mem'],
        res_time=config['tools']['diff']['time'],
        lsf_log=os.path.join(logdir, 'diffexp', 'limma.lsf.log')
    log:
        os.path.join(logdir, 'diffexp', 'limma.log')
    benchmark:
        os.path.join(logdir, 'diffexp', 'limma.benchmark.log')
    conda:
        "../envs/wf_rna_seq_R_limma.yml"
    shell:
        """
        Limma.R {input.counts} {input.conditions} {output.diff_counts}
        """

rule deseq2:
    input:
        counts=os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.tsv'),
        conditions=config['general']['conditions']
    output:
        diff_counts=os.path.join(outdir, 'diffexp','deseq','DE_deseq.tsv')
    threads: 
        deseq2_threads
    log:
        os.path.join(logdir, 'diffexp', 'deseq2.log')
    benchmark:
        os.path.join(logdir, 'diffexp', 'deseq2.benchmark.log')
    params:
        groupA=config['general']['diffexp']['contrasts']['groupA'],
        groupB=config['general']['diffexp']['contrasts']['groupB'],
        res_mem=config['tools']['deseq']['mem'],
        res_time=config['tools']['deseq']['time'],
        lsf_log=os.path.join(logdir, 'diffexp', 'deseq2.lsf.log')
    conda:
        "../envs/wf_rna_seq_R.yml"
    shell:
        """
        deseq2.R {input.counts} {input.conditions} {params.groupA} {params.groupB} {output.diff_counts}
        """
