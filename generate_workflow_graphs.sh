# The commands I am trying on snakemake file
# For plotting the rule graph for the given snakefile
snakemake -s Snakefile --rulegraph | dot -Tpng > rulegraph.png
