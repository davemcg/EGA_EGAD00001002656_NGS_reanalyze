#!/bin/bash

# to run snakemake as batch job
# run in the data folder for this project

module load snakemake || exit 1

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s /home/$USER/git/EGA_EGAD00001002656_NGS_reanalyze/Snakefile_variant_prioritization \
-pr --local-cores 2 --jobs 1999 \
--configfile $1 \
--cluster-config /home/$USER/git/EGA_EGAD00001002656_NGS_reanalyze/cluster_variant_prioritization.json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k 

