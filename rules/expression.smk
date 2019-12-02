rule simple_prep_anno:
    input:
        anno=config['general']['paths']['annotation'],
        genome=config['general']['paths']['genome']
    output:
        config['general']['paths']['annotation'] + '.mask_go.exons.hdf5'
    params:
        res_mem=config['tools']['prep_anno']['mem'],
        res_time=config['tools']['prep_anno']['time'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', 'prep_annonation.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', 'prep_annonation.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "count_expression_prep_anno -a {input.anno} -g {input.genome} -m"  

rule simple_prep_anno_non_alt:
    input:
        anno=config['general']['paths']['annotation'],
        genome=config['general']['paths']['genome']
    output:
        config['general']['paths']['annotation'] + '.mask_go.mask_ao.exons.hdf5'
    params:
        res_mem=config['tools']['prep_anno']['mem'],
        res_time=config['tools']['prep_anno']['time'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', 'prep_annonation.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', 'prep_annonation.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "count_expression_prep_anno -a {input.anno} -g {input.genome} -m -M"


rule simple_counting:
    input: 
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
        bam_index=os.path.join(outdir, 'alignments', '{sample}.all.bam.bai'),
        anno=config['general']['paths']['annotation'] + '.mask_go.exons.hdf5'
    output:
        os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.tsv')
    params:
        annotation=config['general']['paths']['annotation'],
        sample='{sample}',
        res_mem=config['tools']['simple_counting']['mem'],
        res_time=config['tools']['simple_counting']['time'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.expression.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.expression.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "count_expression -a {params.annotation} -A {input.bam} -o {output} -H {params.sample} -v -B -m"


rule simple_counting_non_alt:
    input:  
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
        bam_index=os.path.join(outdir, 'alignments', '{sample}.all.bam.bai'),
        anno=config['general']['paths']['annotation'] + '.mask_go.mask_ao.exons.hdf5'
    output:
        os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.non_alt.tsv')
    params:
        annotation=config['general']['paths']['annotation'],
        sample='{sample}',
        res_mem=config['tools']['simple_counting']['mem'],
        res_time=config['tools']['simple_counting']['time'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.expression.non_alt.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.expression.non_alt.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "count_expression -a {params.annotation} -A {input.bam} -o {output} -H {params.sample} -v -B -m -M"


rule compute_libsize:
    input:
        os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.tsv')
    output:
        os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.libsize.tsv')
    params:
        annotation=config['general']['paths']['annotation'],
        res_mem=config['tools']['simple_counting']['mem_merge'],
        res_time=config['tools']['simple_counting']['time_merge'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.libsize.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.libsize.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "compute_lib_size -i {input} -a {params.annotation} --coding > {log} 2>&1"
        

rule compute_libsize_non_alt:
    input:
        os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.non_alt.tsv')
    output:
        os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.non_alt.libsize.tsv')
    params:
        annotation=config['general']['paths']['annotation'],
        res_mem=config['tools']['simple_counting']['mem_merge'],
        res_time=config['tools']['simple_counting']['time_merge'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.libsize.non_alt.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', '{sample}.all.libsize.non_alt.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "compute_lib_size -i {input} -a {params.annotation} --coding > {log} 2>&1"


rule merge_expression_simple:
    input:
        expand(os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.tsv'), sample=samples['sample'].unique()) 
    output:
        os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.hdf5')
    params:
        res_mem=config['tools']['simple_counting']['mem_merge'],
        res_time=config['tools']['simple_counting']['time_merge'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', 'all_samples.expression.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', 'all_samples.expression.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "collect_counts -i {input} -o {output} -v > {log} 2>&1"


rule merge_expression_simple_non_alt:
    input:
        expand(os.path.join(outdir, 'expression', 'simple_counting', '{sample}.all.expression.non_alt.tsv'), sample=samples['sample'].unique()) 
    output:
        os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.hdf5')
    params:
        res_mem=config['tools']['simple_counting']['mem_merge'],
        res_time=config['tools']['simple_counting']['time_merge'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.lsf.log')
    log:
        os.path.join(logdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "collect_counts -i {input} -o {output} -v > {log} 2>&1"
 

rule convert_hdf5_to_tsv:
    input:
        os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.hdf5')
    output:
        os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.tsv')
    params:
        annotation=config['general']['paths']['annotation'],
        res_mem=config['tools']['simple_counting']['mem_merge'],
        res_time=config['tools']['simple_counting']['time_merge'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', 'all_samples.expression.hdf2tsv.lsf.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "count_hdf2tsv -i {input} -o {output} -v"
 

rule convert_hdf5_to_tsv_non_alt:
    input:
        os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.hdf5')
    output:
        os.path.join(outdir, 'expression', 'simple_counting', 'all_samples.expression.non_alt.tsv')
    params:
        annotation=config['general']['paths']['annotation'],
        res_mem=config['tools']['simple_counting']['mem_merge'],
        res_time=config['tools']['simple_counting']['time_merge'],
        lsf_log=os.path.join(logdir, 'expression', 'simple_counting', 'all_samples.expression.hdf2tsv.non_alt.lsf.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        1
    shell:
        "count_hdf2tsv -i {input} -o {output} -v"


"""
RSEM BAM prep to run RSEM_calculate_expression
"""
rule rsem_bam_prep:
    input:
        bam=os.path.join(outdir, 'alignments_transcript', '{sample}.all.transcript.sort.bam')
    params:
        prefix=os.path.join(outdir, 'expression', 'rsem', '{sample}.all.transcript.sort.rsem'),
        scratch = config['tools']['RSEM_calc_expr']['scratch'],
        res_mem = config['tools']['RSEM_calc_expr']['mem'],
        res_time = config['tools']['RSEM_calc_expr']['time'],
        lsf_log=os.path.join(logdir, 'expression', 'rsem', '{sample}_rsem-prep.lsf.log')
    output:
        outfile=os.path.join(outdir, 'expression', 'rsem', '{sample}.all.transcript.sort.rsem.bam'),
    log:
        os.path.join(logdir, 'expression', 'rsem', '{sample}_rsem-prep.log')
    benchmark:
        os.path.join(logdir, 'expression', 'rsem', '{sample}_rsem-prep.log')
    conda:
        "../envs/rsem.yml"
    threads: config['tools']['RSEM_calc_expr']['threads']
    shell:
        '{config[tools][convert-sam-for-rsem][call]} -p {threads} {input.bam} {params.prefix} > {log}'

"""
Calculate expression using the tool RSEM_calc_expression
"""
rule rsem_calc_expression:
    input:
        bam=os.path.join(outdir, 'expression', 'rsem', '{sample}.all.transcript.sort.rsem.bam')
    params:
        prefix=os.path.join(outdir, 'expression', 'rsem', '{sample}'),
        reference=config['general']['paths']['rsem_reference'],
        variousParams = config['tools']['RSEM_calc_expr']['variousParams'],
        scratch = config['tools']['RSEM_calc_expr']['scratch'],
        res_mem = config['tools']['RSEM_calc_expr']['mem'],
        res_time = config['tools']['RSEM_calc_expr']['time'],
        lsf_log=os.path.join(logdir, 'expression', 'rsem', '{sample}_rsem.lsf.log')
    output:
        iso=os.path.join(outdir, 'expression', 'rsem', '{sample}.isoforms.results'),
        gene=os.path.join(outdir, 'expression', 'rsem', '{sample}.genes.results'),
    log:
        os.path.join(logdir, 'expression', 'rsem', '{sample}_rsem.log')
    benchmark:
        os.path.join(logdir, 'expression', 'rsem', '{sample}_rsem.benchmark.log')
    conda:
        "../envs/rsem.yml"
    threads: config['tools']['RSEM_calc_expr']['threads']
    shell:
        '{config[tools][RSEM_calc_expr][call]} {params.variousParams} -p {threads} --bam {input.bam} {params.reference} {params.prefix}'
        

rule htseq:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam')
    params:
        outdir = os.path.join(outdir, 'expression', 'htseq'),
        annotation = config['general']['paths']['annotation'],
        strandedness = config['general']['strandedness'],
        scratch = config['tools']['htseq']['scratch'],
        res_mem = config['tools']['htseq']['mem'],
        res_time = config['tools']['htseq']['time'],
        lsf_log=os.path.join(logdir, 'expression', 'htseq', '{sample}_htseq.lsf.log')
    output:
        outfile=os.path.join(outdir, 'expression', 'htseq', '{sample}.htseq_counts.txt'),
    log:
        os.path.join(logdir, 'expression', 'htseq', '{sample}_htseq.log')
    benchmark:
        os.path.join(logdir, 'expression', 'htseq', '{sample}_htseq.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads: config['tools']['htseq']['threads']
    shell:  
        ('{config[tools][htseq][call]} -f bam -m intersection-nonempty --stranded={params.strandedness} \
        --idattr gene_id {input.bam} {params.annotation} > {output.outfile} 2> {log} ')
