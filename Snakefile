from os.path import join

configfile: "/data/mcgaugheyd/projects/nei/mcgaughey/EGA_EGAD00001002656_7n/xab/config.yaml"

def return_ID(wildcards):
    # returns the ID in the read group from the header
    import subprocess
    import glob
    cram_file = glob.glob('faux_cram/' + wildcards + '.*')[0]
    command = 'module load samtools; samtools view -H ' + cram_file + ' | grep ^@RG'
    RG_info = subprocess.check_output(command, shell = True)
    RG_info = RG_info.decode().split('\n')[:-1]
    rg_id = []
    for RG in RG_info:
        RG = RG.split('\t')
        the_id = [x for x in RG if 'ID' in x][0]
        the_id = the_id.replace('ID:','')
        rg_id.append(the_id)
    return(rg_id)
    
def lane_bam_names(wildcards):
    # creates the filenames for each lane bam possibility for a sample
    rg_ids = return_ID(str(wildcards))
    lane_bam_files = []
    for rg_id in rg_ids:
        lane_bam_files.append('bam/realigned/{}.ID{}.realigned.bam'.format(wildcards, rg_id))
    return(lane_bam_files)

def build_RG(wildcards):
    # builds a sam header for the lane bam
    import subprocess
    lane_file = wildcards
    command = 'samtools view -H ' + lane_file + ' | grep ^@RG'
    RG_info = subprocess.check_output(command, shell = True)
    RG_info = RG_info.decode().split('\t')
    pl = [x for x in RG_info if 'PL' in x][0].rstrip()
    sm = [x for x in RG_info if 'SM' in x][0].rstrip()
    rg_id = [x for x in RG_info if 'ID' in x][0].rstrip()
    pu = [x for x in RG_info if 'PU' in x][0].rstrip()
    new_RG = '\\@RG\\\\t' + str(rg_id) + '\\\\t' + str(pu) + '\\\\t' + str(sm) + '\\\\t' + str(pl)
    return(new_RG)

def chr_GVCF_to_single_GVCF(wildcards):
	# creates the filenames for the chr level GVCFs to use to concatenate to a single file
	# ensures that input GVCF chrs are provided in order (same as CHRS) below
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('GVCFs/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.g.vcf.gz')
	return(sample_by_chr)

def return_correct_chr_set(wildcards):
	# replaces 'MT_contigs' with actual values	
	MT_CONTIGS="MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"
	chr = str(wildcards.chr)
	if chr == 'MT_contigs':
		chr = MT_CONTIGS
	return(chr)

 	
(SAMPLES, FILE_ENDINGS) = glob_wildcards(join('faux_cram/', '{sample}.ba{file_ending}'))
CHRS=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT_contigs"]
#CHRS=["1","2","3"]

MT_CONTIGS="MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"

wildcard_constraints:
	sample="_EGAR.*_[A-Z]\d{6}_[A|B]"

rule all:
	input:
		expand('GVCFs/{sample}.g.vcf.gz', sample=SAMPLES),
	#	expand('GATK_metrics/{sample}__{chr}.BQSRplots.pdf', sample=SAMPLES, chr=CHRS),
		'GATK_metrics/multiqc_report',
		'fastqc/multiqc_report'

rule globus_cram_transfer_from_Arges:
	input:
		'faux_cram/{sample}.bam.cram'
	output:
		temp('cram/{sample}.bam.cram')
	resources:
		parallel=1
	shell:
		"""
		globus_out=$(globus transfer --sync-level size \
			d0960b02-c5f7-11e5-9a36-22000b96db58:{config[orig_data_path]}/{wildcards.sample}.bam.cram \
			e2620047-6d04-11e5-ba46-22000b92c6ec:{config[project_path]}/{output})
		id=$(echo $globus_out | rev | cut -f1 -d ' ' | rev)
		globus task wait $id
		"""

rule globus_bam_transfer_from_Arges:
    input:
        'faux_cram/{sample}.bam'
    output:
        temp('cram/{sample}.bam')
    resources:
        parallel=1
    shell:
        """
        globus_out=$(globus transfer --sync-level size \
            d0960b02-c5f7-11e5-9a36-22000b96db58:{config[orig_data_path]}/{wildcards.sample}.bam \
            e2620047-6d04-11e5-ba46-22000b92c6ec:{config[project_path]}/{output})
        id=$(echo $globus_out | rev | cut -f1 -d ' ' | rev)
        globus task wait $id
        """
		
rule split_original_cram_by_rg:
	input:
		'cram/{sample}.bam.cram'
	output:
		temp('bam/lane_bam/{sample}/{sample}.ID{rg_id}.bam')
	shell:
		"""
		# ref cache goes to ~/ by default. TERRIBLE
		export REF_CACHE=/lscratch/$SLURM_JOB_ID/hts-refcache
		module load {config[samtools_version]}
		samtools view -b -r {wildcards.rg_id} {input} > {output}
		"""

rule split_original_bam_by_rg:
    input:
        'cram/{sample}.bam'
    output:
        temp('bam/lane_bam/{sample}/{sample}.ID{rg_id}.bam')
    shell:
        """
        # ref cache goes to ~/ by default. TERRIBLE
        export REF_CACHE=/lscratch/$SLURM_JOB_ID/hts-refcache
        module load {config[samtools_version]}
        samtools view -b -r {wildcards.rg_id} {input} > {output}
        """

rule align:
    input:
        'bam/lane_bam/{sample}/{sample}.ID{rg_id}.bam'
    output:
        temp('bam/realigned/{sample}.ID{rg_id}.realigned.bam')
    threads: 16 
    run:
        import subprocess
        rg = build_RG(str(input))
        call = 'module load ' + str(config["samtools_version"]) + '; ' \
               'module load ' + str(config["bwa_version"]) + '; \
				mkdir -p /scratch/mcgaugheyd/$SLURM_JOB_ID/; \
                export REF_CACHE=/scratch/mcgaugheyd/$SLURM_JOB_ID/hts-refcache;  \
                samtools collate -uOn 128 ' + str(input) + ' /scratch/mcgaugheyd/$SLURM_JOB_ID/TMP_' + str(wildcards.sample) + '_ID' + str(wildcards.rg_id) + ' | \
                samtools fastq - | \
                bwa mem -M -t ' + str(threads) + ' -B 4 -O 6 -E 1 -M -p -R ' + str(rg) + \
                    ' ' + str(config["bwa_genome"]) + ' - | \
                    samtools view -1 - > ' + str(output)
        print(call)
        subprocess.check_call(call, shell=True)

rule merge_RG_bams_back_together:
    input:
        lane_bam_names
    output:
        temp('bam/{sample}.realigned.bam')
    threads: 2
    shell:
        """
        module load {config[picard_version]}
        picard_i=""
        for bam in {input}; do
            picard_i+=" I=$bam"
        done
        java -Xmx20g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
                MergeSamFiles \
                TMP_DIR=/lscratch/$SLURM_JOB_ID \
                $picard_i \
                O={output}
        """

rule build_index:
	input:
		'bam/{sample}.realigned.bam'
	output:
		temp('bam/{sample}.realigned.bam.bai')
	threads: 2
	shell:
		"""
		export REF_CACHE=/lscratch/$SLURM_JOB_ID/hts-refcache
		module load {config[samtools_version]}
		samtools index {input} {output}
		"""

rule fastqc:
	input:
		'bam/{sample}.realigned.bam'
	output:
		'fastqc/{sample}'
	threads: 8
	shell:
		"""
		module load fastqc
		mkdir -p fastqc 
		mkdir fastqc/{wildcards.sample}
		fastqc -t {threads} -o {output} {input}
		"""
rule split_bam_by_chr:
	input:
		bam = 'bam/{sample}.realigned.bam',
		bai = 'bam/{sample}.realigned.bam.bai'
	output:
		temp('bam/chr_split/{sample}/{sample}__{chr}.realigned.bam')
	threads: 2
	params:
		chromosome = return_correct_chr_set
	shell:
		"""
		module load {config[samtools_version]}
		samtools view -bh {input.bam} {params.chromosome}  > {output}
		"""
rule picard_clean_sam:
# "Soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads"
	input:
		'bam/chr_split/{sample}/{sample}__{chr}.realigned.bam'
	output:
		temp('bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.bam')
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CleanSam \
			TMP_DIR=/lscratch/$SLURM_JOB_ID \
			INPUT={input} \
			OUTPUT={output}
		"""

rule picard_fix_mate_information:
# "Verify mate-pair information between mates and fix if needed."
# also coord sorts
	input:
		'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.bam'
	output:
		temp('bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.bam')
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
		FixMateInformation \
			SORT_ORDER=coordinate \
			INPUT={input} \
			OUTPUT={output}
		"""

rule picard_mark_dups:
# Mark duplicate reads
	input:
		'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.bam'
	output:
		bam = temp('bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.bam'),
		metrics = 'GATK_metrics/{sample}__{chr}.markDup.metrics'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MarkDuplicates \
			INPUT={input} \
			OUTPUT={output.bam} \
			METRICS_FILE={output.metrics}
		"""

rule picard_bam_index:
# Build bam index
	input:
		'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.bam'
	output:
		temp('bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.bam.bai')
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
		BuildBamIndex \
			INPUT={input} \
			OUTPUT={output}
		"""

rule gatk_realigner_target:
# identify regions which need realignment
	input:
		bam = 'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.bam',
		bai = 'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.bam.bai'
	output:
		temp('bam/chr_split/{sample}/{sample}__{chr}.forIndexRealigner.intervals')
	threads: 2
	params:
		chromosome = return_correct_chr_set
	shell:
		"""
		chromosome_interval=$(sed 's/^\|\s/ -L /g' <(echo {params.chromosome}))
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g RealignerTargetCreator  \
			-R {config[ref_genome]}  \
			$chromosome_interval \
			-I {input.bam} \
			--known {config[1000g_indels]} \
			--known {config[mills_gold_indels]} \
			-o {output}
		"""

rule gatk_indel_realigner:
# realigns indels to improve quality
	input:
		bam = 'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.bam',
		bai = 'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.bam.bai',
		targets = 'bam/chr_split/{sample}/{sample}__{chr}.forIndexRealigner.intervals'
	output:
		temp('bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam')
	threads: 2
	params:
		chromosome = return_correct_chr_set
	shell:
		"""
		chromosome_interval=$(sed 's/^\|\s/ -L /g' <(echo {params.chromosome}))
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g IndelRealigner \
			-R {config[ref_genome]} \
			$chromosome_interval \
			-I {input.bam} \
			--knownAlleles {config[1000g_indels]} \
			--knownAlleles {config[mills_gold_indels]} \
			-targetIntervals {input.targets} \
			-o {output} 
		"""

rule gatk_base_recalibrator:
# recalculate base quality scores
	input:
		'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam'
	output:
		'GATK_metrics/{sample}__{chr}.recal_data.table1'
	threads: 2
	params:
		chromosome = return_correct_chr_set
	shell:
		"""
		chromosome_interval=$(sed 's/^\|\s/ -L /g' <(echo {params.chromosome}))
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g BaseRecalibrator  \
			-R {config[ref_genome]} \
			$chromosome_interval \
			-I {input} \
			--knownSites {config[1000g_indels]} \
			--knownSites {config[mills_gold_indels]} \
			--knownSites {config[dbsnp_var]} \
			-o {output}
		"""

rule gatk_print_reads:
# print out new bam with recalibrated scoring
	input:
		bam = 'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam',
		bqsr = 'GATK_metrics/{sample}__{chr}.recal_data.table1'
	output:
		temp('bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.gatk_realigner.recalibrated.bam')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g PrintReads \
			-R {config[ref_genome]} \
			-I {input.bam} \
			-BQSR {input.bqsr} \
			-o {output}
		"""

rule gatk_base_recalibrator2:
# recalibrate again
	input:
	    bam = 'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam',
		bqsr = 'GATK_metrics/{sample}__{chr}.recal_data.table1'
	output:
		'GATK_metrics/{sample}__{chr}.recal_data.table2'
	threads: 2
	params:
		chromosome = return_correct_chr_set
	shell:
		"""
		chromosome_interval=$(sed 's/^\|\s/ -L /g' <(echo {params.chromosome}))
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g BaseRecalibrator  \
			-R {config[ref_genome]} \
			$chromosome_interval \
			-I {input.bam} \
			--knownSites {config[1000g_indels]} \
			--knownSites {config[mills_gold_indels]} \
			--knownSites {config[dbsnp_var]} \
			-BQSR {input.bqsr} \
			-o {output}
			"""

rule gatk_analyze_covariates:
	input:
		one = 'GATK_metrics/{sample}__{chr}.recal_data.table1',
		two = 'GATK_metrics/{sample}__{chr}.recal_data.table2'
	output:
		'GATK_metrics/{sample}__{chr}.BQSRplots.pdf'
	threads: 2
	params:
		chromosome = return_correct_chr_set
	shell:
		"""
		chromosome_interval=$(sed 's/^\|\s/ -L /g' <(echo {params.chromosome}))
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g AnalyzeCovariates \
			-R {config[ref_genome]} \
			$chromosome_interval \
			-before {input.one} \
			-after {input.two} \
			-plots {output}
		"""

rule gatk_haplotype_caller:
# call gvcf
	input:
		bam = 'bam/chr_split/{sample}/{sample}__{chr}.realigned.CleanSam.sorted.markDup.gatk_realigner.recalibrated.bam',
		bqsr = 'GATK_metrics/{sample}__{chr}.recal_data.table1' 
	output:
		temp('GVCFs/chr_split/{sample}/{sample}__{chr}.g.vcf.gz')
	threads: 2
	params:
		chromosome = return_correct_chr_set
	shell:
		"""
		chromosome_interval=$(sed 's/^\|\s/ -L /g' <(echo {params.chromosome}))
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g HaplotypeCaller  \
			-R {config[ref_genome]} \
			$chromosome_interval \
			-I {input.bam} \
			--emitRefConfidence GVCF \
			-BQSR {input.bqsr} \
			-o {output}
		"""

rule picard_merge_gvcfs:
# merge chr split gvcf back into one gvcf per sample
	input:
		chr_GVCF_to_single_GVCF
	output:
		'GVCFs/{sample}.g.vcf.gz'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		cat_inputs_i=""
		for gvcf in {input}; do
			cat_inputs_i+="I=$gvcf "; done
		java -Xmx15g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeVcfs \
			$cat_inputs_i \
			O={output}
		"""

rule multiqc_gatk:
# run multiqc on recalibrator metrics
	input:
		expand('GATK_metrics/{sample}__{chr}.recal_data.table1',sample=SAMPLES, chr=CHRS),
		expand('GATK_metrics/{sample}__{chr}.recal_data.table2', sample=SAMPLES, chr=CHRS)
	output:
		'GATK_metrics/multiqc_report'
	shell:
		"""
		module load multiqc
		multiqc -f -o {output} GATK_metrics
		"""

rule multiqc_fastqc:
	input:
		expand('fastqc/{sample}', sample=SAMPLES)
	output:
		fastqc = 'fastqc/multiqc_report'
	shell:
		"""
		module load multiqc
		multiqc -f -o {output} fastqc/
		"""
