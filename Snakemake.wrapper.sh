#!/bin/bash

# to run snakemake as batch job
# run in the data folder for this project
# on biowulf2:
# /data/mcgaugheyd/projects/nei/mcgaughey/EGA_EGAD00001002656

cd /data/mcgaugheyd/projects/nei/mcgaughey/EGA_EGAD00001002656

module load smakemake || exit 1

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s ~/git/EGA_EGAD00001002656_NGS_reanalyze/Snakemake \
-pr --local-cores 2 --jobs 999 \
--cluster-config ~/git/EGA_EGAD00001002656_NGS_reanalyze/cluster.json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete
