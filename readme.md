# EGA_EGAD00001002656 Reanalysis Project
Snakemake pipeline to (re)call genotypes from EGA EGAD00001002656 using a bwa mem / GATK based approach

Two notes:
1. The faux_cram directory is populated with the headers of the cram/bam files (that way you can transfer the crucial metadata to the server to determine the processing steps without having to transfer the entire file first). You generate the faux_cram files by running `samtools view -bH original.cram > faux_original.cram`
2. Make sure you have created a 00log folder to hold the log files. Snakemake won't do this.
