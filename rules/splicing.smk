### general variables
genes_graph_conf_path = os.path.join(outdir, 'splicing', 'spladder', 'genes_graph_conf')

rule spladder_single_graphs:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'),
        bam_index=os.path.join(outdir, 'alignments', '{sample}.all.bam.bai') 
    output:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.pickle'
    params:
        outdir=os.path.join(outdir, 'splicing'),
        annotation=config['general']['paths']['annotation'],
        confidence=config['tools']['spladder']['confidence'],
        readlen=config['tools']['spladder']['readlen'],
        res_mem=config['tools']['spladder']['mem'],
        res_time=config['tools']['spladder']['time'],
        lsf_log=os.path.join(logdir, 'splicing', 'single_graphs_build.' + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.lsf.log')
    log:
       os.path.join(logdir, 'splicing', 'single_graphs_build.' + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads: 
        config['tools']['spladder']['threads']
    shell:
        ('{config[tools][spladder][call]} build -v -b {input.bam} \
-o {params.outdir} \
-a {params.annotation} \
-c {params.confidence} \
--merge-strat single \
--parallel {threads} \
--readlen {params.readlen} \
{config[tools][spladder][otherParams]} > {log} 2>&1')

rule spladder_single_quantify:
    input:
        bam=os.path.join(outdir, 'alignments', '{sample}.all.bam'), 
        bam_index=os.path.join(outdir, 'alignments', '{sample}.all.bam.bai'), 
        pickle=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.pickle'
    output:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.hdf5'
    params:
        outdir=os.path.join(outdir, 'splicing'),
        annotation=config['general']['paths']['annotation'],
        confidence=config['tools']['spladder']['confidence'],
        readlen=config['tools']['spladder']['readlen'],
        res_mem=config['tools']['spladder']['mem'],
        res_time=config['tools']['spladder']['time'],
        lsf_log=os.path.join(logdir, 'splicing', 'single_graphs_quantify.' + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.lsf.log')
    log:
       os.path.join(logdir, 'splicing', 'single_graphs_quantify.' + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads: 
        config['tools']['spladder']['threads']
    shell:
        ('{config[tools][spladder][call]} build -v -b {input.bam} \
-o {params.outdir} \
-a {params.annotation} \
-c {params.confidence} \
--merge-strat single \
--parallel {threads} \
--readlen {params.readlen} \
--quantify-graph \
{config[tools][spladder][otherParams]} > {log} 2>&1')


rule spladder_merge_graphs:
    input:
        bam=expand(os.path.join(outdir, 'alignments', '{sample}.all.bam'), sample=samples['sample'].unique()),
        bam_index=expand(os.path.join(outdir, 'alignments', '{sample}.all.bam.bai'), sample=samples['sample'].unique()),
        pickle=expand(genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.pickle', sample=samples['sample'].unique())
    output:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.merge_graphs.pickle',
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.merge_graphs.count.hdf5'
    params:
        outdir=os.path.join(outdir, 'splicing'),
        annotation=config['general']['paths']['annotation'],
        confidence=config['tools']['spladder']['confidence'],
        readlen=config['tools']['spladder']['readlen'],
        res_mem=config['tools']['spladder']['mem'],
        res_time=config['tools']['spladder']['time'],
        lsf_log=os.path.join(logdir, 'splicing', 'merge_graphs.conf' + str(config['tools']['spladder']['confidence'])  + '.lsf.log')
    log:
       os.path.join(logdir, 'splicing', 'merge_graphs.conf' + str(config['tools']['spladder']['confidence'])  + '.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads: 
        config['tools']['spladder']['threads']
    shell:
        ('bams=$(echo {input.bam} | tr " " ","); {config[tools][spladder][call]} build -v -b $bams \
-o {params.outdir} \
-a {params.annotation} \
-c {params.confidence} \
--merge-strat merge_graphs \
--parallel {threads} \
--readlen {params.readlen} \
 --quantify-graph \
{config[tools][spladder][otherParams]} > {log} 2>&1')

rule splice_burden_projection_single:
    input:
        junc_db=config['general']['paths']['burden_outgroup_junction_db'],
        sample_pickle=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.pickle',
        sample_count=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.hdf5'
    output:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.junc_db_proj.hdf5'
    params:
        res_mem=config['tools']['burden']['mem'],
        res_time=config['tools']['burden']['time'],
        lsf_log=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.junc_db_proj.lsf.log'
    log:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.junc_db_proj.log'
    conda:
        "../envs/wf_rna_seq.yml"
    threads: 
        1
    shell:
        'splice_burden_project {input.junc_db} {input.sample_pickle} {input.sample_count} {output}'

rule splice_burden_compute_single:
    input:
        anno_junc=config['general']['paths']['burden_anno_junction_db'],
        sample_count=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.hdf5',
        sample_pickle=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.pickle',
        sample_project=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.junc_db_proj.hdf5',
        tcga_junc=config['general']['paths']['tcga_neojunctions']
    params:
        res_mem=config['tools']['burden']['mem'],
        res_time=config['tools']['burden']['time'],
        lsf_log=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.burden.glob{globsum}.lsf.log',
        globsum='{globsum}'
    output:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.G0.01.globsum{globsum}.conf3_neojunctions.tsv.gz'
    log:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.burden.glob{globsum}.log'
    conda:
        "../envs/wf_rna_seq.yml"
    threads: 
        1
    shell:
        'splice_burden_compute {input.anno_junc} {input.sample_count} {input.sample_pickle} {input.sample_project} {input.tcga_junc} {params.globsum}'

rule plot_splice_burden_single:
    input:
        juncs=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.G0.01.globsum{globsum}.conf3_neojunctions.tsv.gz',
        tcga_count=config['general']['paths']['tcga_neojunction_counts'],
        tcga_meta=config['general']['paths']['tcga_metadata']
    output:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.count.G0.01.globsum{globsum}.conf3_neojunctions.pdf'
    params:
        res_mem=config['tools']['burden_plot']['mem'],
        res_time=config['tools']['burden_plot']['time'],
        lsf_log=genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.burden.glob{globsum}.plotting.lsf.log',
        sample='{sample}'
    log:
        genes_graph_conf_path + str(config['tools']['spladder']['confidence'])  + '.{sample}.all.burden.glob{globsum}.plotting.log'
    conda:
        "../envs/wf_rna_seq.yml"
    threads: 
        1
    shell:
        'splice_burden_plot_tcga {input.juncs} {input.tcga_count} {input.tcga_meta} {params.sample} {output}'

