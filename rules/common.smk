def is_single(sample, lane):
    """Return true if given lane of sample has only one fastq file."""
    return pd.isnull(samples.loc[(sample, lane), 'fq2'])


def get_fastqs(wildcards):
    """Return original fastq file name(s) for the current sample and lane"""
    if is_single(**wildcards):
        return os.path.join(config['general']['paths']['fastqdir'], samples.loc[(wildcards.sample, wilacards.lane), 'fq1'])
    else:
        return [os.path.join(config['general']['paths']['fastqdir'], samples.loc[(wildcards.sample, wildcards.lane), 'fq1']),
                os.path.join(config['general']['paths']['fastqdir'], samples.loc[(wildcards.sample, wildcards.lane), 'fq2'])]

