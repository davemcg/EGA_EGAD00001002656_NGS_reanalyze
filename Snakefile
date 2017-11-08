from os.path import join

configfile: "/home/mcgaugheyd/git/EGA_EGAD00001002656_NGS_reanalyze/config.yaml"

def return_ID(wildcards):
    # returns the ID in the read group from the header
    import subprocess
    import glob
    cram_file = glob.glob('faux_cram/' + wildcards + '.*')[0]
    command = 'samtools view -H ' + cram_file + ' | grep ^@RG'
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
        lane_bam_files.append('temp/realigned/{}.ID{}.realigned.bam'.format(wildcards, rg_id))
    return(lane_bam_files)

def build_RG(wildcards):
    # builds a sam header for the lane bam
    import subprocess
 #   lane_file = 'temp/' + wildcards.sample + 'RG' + wildcards.rg_id + '.bam'
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

(SAMPLES, FILE_ENDINGS) = glob_wildcards(join('faux_cram/', '{sample}.ba{file_ending}'))
#SAMPLES, = glob_wildcards(join('faux_cram/', '{sample}.bam.cram'))
#SAMPLES = open('/home/mcgaugheyd/git/EGA_EGAD00001002656_NGS_reanalyze/all_cram_names.txt').read().splitlines()

rule all:
	input:
		expand('GVCFs/{sample}.g.vcf.gz', sample=SAMPLES),
		expand('GATK_metrics/{sample}.BQSRplots.pdf', sample=SAMPLES),
		'GATK_metrics/multiqc_report'
		#expand('DELETE.{sample}.DELETE', sample=SAMPLES)

rule globus_cram_transfer_from_Arges:
	input:
		'faux_cram/{sample}.bam.cram'
	output:
		temp('cram/{sample}.bam.cram')
	shell:
		"""
		globus_out=$(globus transfer --sync-level size \
			d0960b02-c5f7-11e5-9a36-22000b96db58:/Volumes/Arges/EGA_RD_WGS/{wildcards.sample}.bam.cram \
			e2620047-6d04-11e5-ba46-22000b92c6ec:/data/mcgaugheyd/{output})
		id=$(echo $globus_out | rev | cut -f1 -d ' ' | rev)
		globus task wait $id
		"""

rule globus_bam_transfer_from_Arges:
    input:
        'faux_cram/{sample}.bam'
    output:
        temp('cram/{sample}.bam')
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
		temp('temp/lane_bam/{sample}.ID{rg_id}.bam')
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
        temp('temp/lane_bam/{sample}.ID{rg_id}.bam')
    shell:
        """
        # ref cache goes to ~/ by default. TERRIBLE
        export REF_CACHE=/lscratch/$SLURM_JOB_ID/hts-refcache
		module load {config[samtools_version]}
        samtools view -b -r {wildcards.rg_id} {input} > {output}
        """
#	message:
#		'echo {read_group_IDs}'
#		'echo {output}'

rule align:
    input:
        'temp/lane_bam/{sample}.ID{rg_id}.bam'
    output:
        temp('temp/realigned/{sample}.ID{rg_id}.realigned.bam')
    threads: 8 
    params:
        samtools = '{config[samtools_version]}',
        bwa = '{config[bwa_version]}',
        ref = '{config[ref_genome]}'
    run:
        import subprocess
        rg = build_RG(str(input))
        call = 'module load ' + str(params.samtools) ; \
               'module load ' + str(params.bwa) ; \
               'samtools collate -uOn 128 ' + str(input) + ' /lscratch/$SLURM_JOB_ID/TMP_' + str(wildcards.sample) + '_ID' + str(wildcards.rg_id) + ' | \
                samtools fastq - | \
                bwa mem -M -t ' + str(threads) + ' -B 4 -O 6 -E 1 -M -p -R ' + str(rg) + \
                    ' ' + str(params.ref) + ' - | \
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

rule picard_clean_sam:
# "Soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads"
	input:
		'bam/{sample}.realigned.bam'
	output:
		temp('bam/{sample}.realigned.CleanSam.bam')
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
		'bam/{sample}.realigned.CleanSam.bam'
	output:
		temp('bam/{sample}.realigned.CleanSam.sorted.bam')
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
		'bam/{sample}.realigned.CleanSam.sorted.bam'
	output:
		bam = temp('bam/{sample}.realigned.CleanSam.sorted.markDup.bam'),
		metrics = 'GATK_metrics/{sample}.markDup.metrics'
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
		'bam/{sample}.realigned.CleanSam.sorted.markDup.bam'
	output:
		temp('bam/{sample}.realigned.CleanSam.sorted.markDup.bam.bai')
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
		bam = 'bam/{sample}.realigned.CleanSam.sorted.markDup.bam',
		bai = 'bam/{sample}.realigned.CleanSam.sorted.markDup.bam.bai'
	output:
		temp('bam/{sample}.forIndexRealigner.intervals')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g RealignerTargetCreator  \
			-R {config[ref_genome]}  \
			-I {input.bam} \
			--known {config[1000g_indels]} \
			--known {config[mills_gold_indels]} \
			-o {output}
		"""

rule gatk_indel_realigner:
# realigns indels to improve quality
	input:
		bam = 'bam/{sample}.realigned.CleanSam.sorted.markDup.bam',
		bai = 'bam/{sample}.realigned.CleanSam.sorted.markDup.bam.bai',
		targets = 'bam/{sample}.forIndexRealigner.intervals'
	output:
		temp('bam/{sample}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g IndelRealigner \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--knownAlleles {config[1000g_indels]} \
			--knownAlleles {config[mills_gold_indels]} \
			-targetIntervals {input.targets} \
			-o {output} 
		"""

rule gatk_base_recalibrator:
# recalculate base quality scores
	input:
		'bam/{sample}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam'
	output:
		'GATK_metrics/{sample}.recal_data.table1'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g BaseRecalibrator  \
			-R {config[ref_genome]} \
			-I {input} \
			--knownSites {config[1000g_indels]} \
			--knownSites {config[mills_gold_indels]} \
			--knownSites {config[dbsnp_var]} \
			-o {output}
		"""

rule gatk_print_reads:
# print out new band with recalibrated scoring
	input:
		bam = 'bam/{sample}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam',
		bqsr = 'GATK_metrics/{sample}.recal_data.table1'
	output:
		temp('bam/{sample}.realigned.CleanSam.sorted.markDup.gatk_realigner.recalibrated.bam')
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
	    bam = 'bam/{sample}.realigned.CleanSam.sorted.markDup.gatk_realigner.bam',
		bqsr = 'GATK_metrics/{sample}.recal_data.table1'
	output:
		'GATK_metrics/{sample}.recal_data.table2'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g BaseRecalibrator  \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--knownSites {config[1000g_indels]} \
			--knownSites {config[mills_gold_indels]} \
			--knownSites {config[dbsnp_var]} \
			-BQSR {input.bqsr} \
			-o {output}
			"""

rule gatk_analyze_covariates:
	input:
		one = 'GATK_metrics/{sample}.recal_data.table1',
		two = 'GATK_metrics/{sample}.recal_data.table2'
	output:
		'GATK_metrics/{sample}.BQSRplots.pdf'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g AnalyzeCovariates \
			-R {config[ref_genome]} \
			-before {input.one} \
			-after {input.two} \
			-plots {output}
		"""

rule gatk_haplotype_caller:
# call gvcf
	input:
		bam = 'bam/{sample}.realigned.CleanSam.sorted.markDup.gatk_realigner.recalibrated.bam',
		bqsr = 'GATK_metrics/{sample}.recal_data.table1' 
	output:
		'GVCFs/{sample}.g.vcf.gz'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g HaplotypeCaller  \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--emitRefConfidence GVCF \
			-BQSR {input.bqsr} \
			-o {output}
		"""

rule multiqc_gatk:
# run multiqc on recalibrator metrics
	input:
		expand('GATK_metrics/{sample}.recal_data.table1',sample=SAMPLES),
		expand('GATK_metrics/{sample}.recal_data.table2', sample=SAMPLES)
	output:
		'GATK_metrics/multiqc_report'
	shell:
		"""
		module load multiqc
		multiqc -f -o {output} GATK_metrics
		"""
