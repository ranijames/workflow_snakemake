def get_all_lanes(wildcards):
    """Assemble bam file names integrating sample and lane"""
    return expand(os.path.join(outdir, 'alignments', '{sample}_{lane}', '{sample}_{lane}_Aligned.sortedByCoord.out.bam'), lane=samples.loc[wildcards.sample].lane.tolist(), **wildcards)

def get_outdir_star(wildcards):
    """Assemble output directory for STAR"""
    return os.path.join(outdir, 'alignments', '{sample}_{lane}'.format(sample=wildcards.sample, lane=wildcards.lane))

def get_outprefix_star(wildcards):
    """Assemble output prefix for STAR"""
    return os.path.join(get_outdir_star(wildcards), '{sample}_{lane}_'.format(sample=wildcards.sample, lane=wildcards.lane))

def get_all_lanes_transcript(wildcards):
    """Assemble bam file names integrating sample and lane"""
    return expand(os.path.join(outdir, 'alignments_transcript', '{sample}_{lane}', '{sample}_{lane}_Aligned.toTranscriptome.out.bam'), lane=samples.loc[wildcards.sample].lane.tolist(), **wildcards)

def get_outdir_star_transcript(wildcards):
    """Assemble output directory for STAR transcript BAMs"""
    return os.path.join(outdir, 'alignments_transcript', '{sample}_{lane}'.format(sample=wildcards.sample, lane=wildcards.lane))

def get_outprefix_star_transcript(wildcards):
    """Assemble output prefix for STAR transcript BAMs"""
    return os.path.join(get_outdir_star_transcript(wildcards), '{sample}_{lane}_'.format(sample=wildcards.sample, lane=wildcards.lane))


rule index_genome:
    input:
        genome=config['general']['paths']['genome']
    params:
        outdir=config['general']['paths']['star_genome'],
        annotation=config['general']['paths']['annotation'],
        overhang=config['tools']['STAR_index']['sjdbOverhang'],
        res_mem=config['tools']['STAR_index']['mem'],
        res_time=config['tools']['STAR_index']['time'],
        lsf_log=os.path.join(config['general']['paths']['star_genome'], 'generate_index.lsf.log')
    output:
        os.path.join(config['general']['paths']['star_genome'], 'gen_index.done')
    log:
        os.path.join(config['general']['paths']['star_genome'], 'generate_index.log')
    benchmark:
        os.path.join(config['general']['paths']['star_genome'], 'generate_index.benchmark.log')
    threads:
        config['tools']['STAR_index']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['STAR_index']['threads']
    shell:
        ('mkdir -p {params.outdir}; \
STAR \
--genomeDir {params.outdir} \
--genomeFastaFiles {input.genome} \
--sjdbOverhang {params.overhang} \
--sjdbGTFfile {params.annotation} \
--runThreadN {threads} \
{config[tools][STAR_index][otherParams]} > {log} 2>&1 && touch {output}')

rule align:
    input:
        index_done=os.path.join(config['general']['paths']['star_genome'], 'gen_index.done'),
        fastq=get_fastqs
    params:
        outdir=get_outdir_star, 
        outprefix=get_outprefix_star,
        genome=config['general']['paths']['star_genome'],
        overhang=config['tools']['STAR_index']['sjdbOverhang'],
        res_mem=config['tools']['STAR']['mem'],
        res_time=config['tools']['STAR']['time'],
        lsf_log=os.path.join(config['general']['paths']['logdir'], 'alignments', '{sample}_{lane}_align.lsf.log')
    output:
        temp(os.path.join(outdir, 'alignments', '{sample}_{lane}', '{sample}_{lane}_Aligned.out.bam'))
    log:
        os.path.join(logdir, 'alignments', '{sample}_{lane}_align.log')
    benchmark:
        os.path.join(logdir, 'alignments', '{sample}_{lane}_align.benchmark.log')
    threads:
        config['tools']['STAR']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['STAR']['threads']
    shell:
        ('mkdir -p {params.outdir}; \
{config[tools][STAR][call]} \
--genomeDir {params.genome} \
--readFilesIn {input.fastq} \
--runThreadN {threads} \
--sjdbOverhang {params.overhang} \
--outFileNamePrefix {params.outprefix} \
{config[tools][STAR][otherParams]} {config[tools][STAR][passmodeParams]} > {log} 2>&1')
        
rule sort_bam:
    input:
        os.path.join(outdir, 'alignments', '{sample}_{lane}', '{sample}_{lane}_Aligned.out.bam')
    output:
        temp(os.path.join(outdir, 'alignments', '{sample}_{lane}', '{sample}_{lane}_Aligned.sortedByCoord.out.bam'))
    log:
        os.path.join(logdir, 'alignments', '{sample}_{lane}_align.sort.log')
    benchmark:
        os.path.join(logdir, 'alignments', '{sample}_{lane}_align.sort.benchmark.log')
    threads:
        config['tools']['samtools']['threads']
    params:
        res_time=config['tools']['samtools']['time'],
        res_mem=config['tools']['samtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments', '{sample}_{lane}_align.sort.lsf.log'),
        tmp_dir=os.path.join(outdir, 'alignments', '{sample}_{lane}', 'sorting_tmp')
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'mkdir -p {params.tmp_dir} && samtools sort -o {output} {input} > {log} 2>&1 && rm -r {params.tmp_dir}'   

rule merge_lanes:
    input:
        genome_lane_bams=get_all_lanes,
    output:
        genome_bam=temp(os.path.join(outdir, 'alignments', '{sample}.all.bam')) if config['general']['delete_bam'] else os.path.join(outdir, 'alignments', '{sample}.all.bam') 
    log:
        genome_merge_log=os.path.join(logdir, 'alignments', '{sample}_merge.log'),
    params:
        res_time=config['tools']['samtools']['time'],
        res_mem=config['tools']['samtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments', '{sample}_merge.lsf.log')
    threads:
        config['tools']['samtools']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'samtools merge {output.genome_bam} {input.genome_lane_bams} > {log.genome_merge_log} 2>&1'

rule index_bam:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
    output:
        index=temp(os.path.join(outdir, 'alignments', '{sample}.all.bam.bai')) if config['general']['delete_bam'] else os.path.join(outdir, 'alignments', '{sample}.all.bam.bai')
    log:
        index_log=os.path.join(logdir, 'alignments', '{sample}_index.log')
    params:
        res_time=config['tools']['samtools']['time'],
        res_mem=config['tools']['samtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments', '{sample}_index.lsf.log')
    threads:
        config['tools']['samtools']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'samtools index {input.bam} > {log.index_log} 2>&1' 

rule sort_readname_bam:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam')
    output:
        out_bam=os.path.join(outdir, 'alignments', '{sample}.readname_sort.all.bam')
    log:
        sort_readname_log=os.path.join(logdir, 'alignments', '{sample}_readname_sort.log')
    params:
        res_time=config['tools']['samtools']['time'],
        res_mem=config['tools']['samtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments', '{sample}_readname_sort.lsf.log'),
        tmp_dir=os.path.join(outdir, 'alignments', 'sort_tmp')
    threads:
        config['tools']['samtools']['threads_sort']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'mkdir -p {params.tmp_dir};'
        'samtools sort -n -T {params.tmp_dir} --output-fmt BAM -o {output.out_bam} -@ {threads} {input.bam}'

rule compress_bam:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
        genome=config['general']['paths']['genome']
    output:
        cram=os.path.join(outdir, 'alignments', '{sample}.all.cram'),
    params:
        res_time=config['tools']['cramtools']['time'],
        res_mem=config['tools']['cramtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments', '{sample}_compress.lsf.log')
    threads:
        config['tools']['cramtools']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        "samtools view -h -F4 {input.bam} | cramtools cram -O {output.cram} -n --capture-tags 'NM' -L '*8' -R {input.genome}"

rule index_cram:
    input:
        cram=os.path.join(outdir, 'alignments', '{sample}.all.cram')
    output:
        crai=os.path.join(outdir, 'alignments', '{sample}.all.crai')
    params:
        res_time=config['tools']['cramtools']['time'],
        res_mem=config['tools']['cramtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments', '{sample}_index_cram.lsf.log')
    threads:
        config['tools']['cramtools']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        "cramtools index -I {input.cram} -O {output.crai}"


rule align_transcript:
    input:
        index_done=os.path.join(config['general']['paths']['star_genome'], 'gen_index.done'),
        fastq=get_fastqs
    params:
        outdir=get_outdir_star_transcript, 
        outprefix=get_outprefix_star_transcript,
        genome=config['general']['paths']['star_genome'],
        overhang=config['tools']['STAR_index']['sjdbOverhang'],
        res_mem=config['tools']['STAR']['mem'],
        res_time=config['tools']['STAR']['time'],
        lsf_log=os.path.join(logdir, 'alignments_transcript', '{sample}_{lane}_align-transcript.lsf.log')
    output:
        temp(os.path.join(outdir, 'alignments_transcript', '{sample}_{lane}', '{sample}_{lane}_Aligned.toTranscriptome.out.bam'))
    log:
        os.path.join(logdir, 'alignments_transcript', '{sample}_{lane}_align-transcript.log')
    benchmark:
        os.path.join(logdir, 'alignments_transcript', '{sample}_{lane}_align-transcript.benchmark.log')
    threads:
        config['tools']['STAR']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        ('mkdir -p {params.outdir}; \
{config[tools][STAR][call]} \
--genomeDir {params.genome} \
--readFilesIn {input.fastq} \
--runThreadN {threads} \
--sjdbOverhang {params.overhang} \
--outFileNamePrefix {params.outprefix} \
{config[tools][STAR][otherParams]} {config[tools][STAR][quantmodeParams]} > {log} 2>&1')
        
rule merge_lanes_transcript:
    input:
        transcript_lane_bams=get_all_lanes_transcript
    output:
        transcript_bam=temp(os.path.join(outdir, 'alignments_transcript', '{sample}.all.transcript.bam'))
    log:
        os.path.join(logdir, 'alignments_transcript', '{sample}_merge-transcript.log')
    params:
        res_time=config['tools']['samtools']['time'],
        res_mem=config['tools']['samtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments_transcript', '{sample}_merge-transcript.lsf.log')
    threads:
        config['tools']['samtools']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'samtools merge {output.transcript_bam} {input.transcript_lane_bams} > {log} 2>&1'


rule sort_transcript_bam:
    input:
        bam=os.path.join(outdir, 'alignments_transcript', '{sample}.all.transcript.bam'),
    output:
        os.path.join(outdir, 'alignments_transcript', '{sample}.all.transcript.sort.bam')
    log:
        os.path.join(logdir, 'alignments_transcript', '{sample}_sort-transcript.log')
    params:
        res_time=config['tools']['samtools']['time'],
        res_mem=config['tools']['samtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments_transcript', '{sample}_sort-trancript.lsf.log')
    threads:
        config['tools']['samtools']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'samtools sort {input.bam} -o {output} > {log} 2>&1' 

rule index_transcript_bam:
    input:
        bam=os.path.join(outdir, 'alignments_transcript', '{sample}.all.transcript.sort.bam'),
    output:
        index=os.path.join(outdir, 'alignments_transcript', '{sample}.all.transcript.sort.bam.bai')
    log:
        os.path.join(logdir, 'alignments_transcript', '{sample}_index-transcript.log')
    params:
        res_time=config['tools']['samtools']['time'],
        res_mem=config['tools']['samtools']['mem'],
        lsf_log=os.path.join(logdir, 'alignments_transcript', '{sample}_index-trancript.lsf.log')
    threads:
        config['tools']['samtools']['threads']
    conda:
        "../envs/wf_rna_seq.yml"
    shell:
        'samtools index {input.bam} > {log} 2>&1' 

