#!/bin/bash

module load gemini
gemini query -q 'select * from samples' EGAD00001002656.GATK.PED_EGAD00001002656.gemini.db | cut -f2 > ~/git/EGA_EGAD00001002656_NGS_reanalyze/data/present_samples.txt
