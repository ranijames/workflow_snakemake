# This rule normalizes the gene expression quantification data of files of type {sample}.genes.results that are given out by rsem_calculate_expression
rule normalize_gene_quant:
    input:
        ingene=os.path.join(outdir, 'expression', 'rsem', '{sample}.genes.results')
    params:
        res_mem=config['tools']['normalize_gene_quant']['mem'],
        res_time=config['tools']['normalize_gene_quant']['time'],
        scratch=config['tools']['normalize_gene_quant']['scratch'],
        variousParams = config['tools']['normalize_gene_quant']['variousParams'],
        lsf_log=os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.results.normalize_gene_quant.lsf.log')
    output:
        outgene=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.results.normalized.txt')
    log:
        os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.results.normalize_gene_quant.log')
    benchmark:
        os.path.join(logdir, 'tcga_boxplot', '{sample}_normalize_gene_quant.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['normalize_gene_quant']['threads']
    shell:
        '{config[tools][normalize_gene_quant][call]} {params.variousParams} -t 1000 -o {output.outgene} {input.ingene} '

# This rule changes the header of the file header of normalized gene output file
rule change_header:
    input:
        infile=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.results.normalized.txt')
    params:
        res_mem=config['tools']['change_header']['mem'],
        res_time=config['tools']['change_header']['time'],
        scratch=config['tools']['change_header']['scratch'],
        lsf_log=os.path.join(logdir, 'tcga-boxplot', '{sample}.genes.results.normalize_change_header.lsf.log'),
        sampleName='{sample}'
    output:
        out=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.results.normalized.header.txt')
    log:
        os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.results.normalize_change_header.log')
    benchmark:
        os.path.join(logdir, 'tcga_boxplot', '{sample}_normalize_change_header.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['change_header']['threads']
    shell:
        '{config[tools][change_header][call]} "1 s/0.0000/{params.sampleName}/" {input.infile} > {output.out}'

# This rule filters genes of interest from the normalized output file 
rule filter_genes:
    input:
        infile=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.results.normalized.header.txt')
    params:
        res_mem=config['tools']['filter_genes']['mem'],
        res_time=config['tools']['filter_genes']['time'],
        scratch=config['tools']['filter_genes']['scratch'],
        lsf_log=os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.results.normalized.genes_filtered.lsf.log'),
    output:
        outfile=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.results.normalized.genes_filtered.txt')
    log:
        os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.results.normalized.genes_filtered.log')
    benchmark:
        os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.results.normalized.genes_filtered.benchmark.log')
    conda:
        "../envs/wf_rna_seq.yml"
    threads:
        config['tools']['filter_genes']['threads']
    shell:
        '{config[tools][filter_genes][call]} resources/tcga_boxplot/genes_of_interest.txt {input.infile} > {output.outfile}'

# This rule parses the tcga pipeline results and cohort values to a useful table format
rule parse_cohort_comparison:
    input:
        infile=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.results.normalized.genes_filtered.txt')
    params:
        res_mem=config['tools']['parse_cohort_comparison']['mem'],
        res_time=config['tools']['parse_cohort_comparison']['time'],
        scratch=config['tools']['parse_cohort_comparison']['scratch'],
        cohort=config['tools']['parse_cohort_comparison']['cohort'],
        lsf_log=os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.parse_cohort.lsf.log')
    output:
        out=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.final_rna_table.txt')
    log:
        os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.parse_cohort.log')
    benchmark:
        os.path.join(logdir, 'tcga_boxplot', '{sample}_parse_cohort_comparison.benchmark.log')
    conda:
        "../envs/tcga_boxplot.yml"
    threads:
        config['tools']['parse_cohort_comparison']['threads']
    shell:
        '{config[tools][parse_cohort_comparison][call]} {params.cohort} {input.infile} {output.out}'

# This rule plots a boxplot of the expression values of all cohort samples and marks the position of the patient sample in the boxplot. It does so for the subset of genes in the list that it is given.
rule boxplot_expression_value:
    input:
        patient=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.results.normalized.genes_filtered.txt'),
        filtered=os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.final_rna_table.txt')
    params:
        res_mem=config['tools']['boxplot_expression_value']['mem'],
        res_time=config['tools']['boxplot_expression_value']['time'],
        scratch=config['tools']['boxplot_expression_value']['scratch'],
        cohort=config['tools']['boxplot_expression_value']['cohort'],
        sampleName='{sample}',
        out_dir=os.path.join(outdir, 'tcga_boxplot', 'png_images'),
        lsf_log=os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.boxplot.lsf.log')
    output:
        out=touch(os.path.join(outdir, 'tcga_boxplot', '{sample}.genes.boxplot.done'))
    log:
        os.path.join(logdir, 'tcga_boxplot', '{sample}.genes.boxplot.log')
    benchmark:
        os.path.join(logdir, 'tcga_boxplot', '{sample}_boxplot.benchmark.log')
    conda:
        "../envs/tcga_boxplot.yml"
    threads:
        config['tools']['boxplot_expression_value']['threads']
    shell:
        'if [ ! -d {params.out_dir} ]; then mkdir {params.out_dir}; fi; Rscript {config[tools][boxplot_expression_value][call]} {params.out_dir} {params.cohort} {input.patient} {input.filtered} {params.sampleName}'
        #; touch {output.out}'


