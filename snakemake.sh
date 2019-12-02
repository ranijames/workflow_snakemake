CONDA_DIR=/cluster/home/thomasti/software/anaconda/envs


snakemake --use-conda --conda-prefix ${CONDA_DIR} --configfile config.yaml --cluster 'bsub -M $(({params.res_mem} * {threads})) -W {params.res_time} -n {threads} -R "rusage[mem={params.res_mem}]" -o {params.lsf_log}' -j 20 -k -p --notemp --latency-wait 120 --rerun-incomplete

