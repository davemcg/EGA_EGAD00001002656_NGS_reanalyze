# EGA_EGAD00001002656 Reanalysis Project
Snakemake pipeline to (re)call genotypes from EGA EGAD00001002656 using a bwa mem / GATK based approach

Working dir on biowulf2:
`/data/mcgaugheyd/projects/nei/mcgaughey/EGA_EGAD00001002656`

Two notes:
1. The faux_cram directory is populated with the headers of the cram/bam files (that way you can transfer the crucial metadata to the server to determine the processing steps without having to transfer the entire file first). You generate the faux_cram files by running `for i in *.bam; do samtools view -h $i | head -n 1000 | samtools view -b  > faux_cram/$i; done` and `for i in *.cram; do samtools view -h $i | head -n 1000 | samtools view -b  > faux_cram/$i; done`
2. Make sure you have created a `00log folder` to hold the log files. Snakemake won't do this.
