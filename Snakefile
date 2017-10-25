from os.path import join


def return_ID(wildcards):
    # returns the ID in the read group from the header
    import subprocess
    cram_file = 'cram/' + wildcards + '.bam.cram'
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

SAMPLES, = glob_wildcards(join('cram/', '{sample}.bam.cram'))

rule all:
    input:
        expand('bam/{sample}.realigned.bam', sample=SAMPLES)
rule split_cram_by_rg:
    input:
        'cram/{sample}.bam.cram'
    output:
        'temp/lane_bam/{sample}.ID{rg_id}.bam'
    shell:
        """
        # ref cache goes to ~/ by default. TERRIBLE
        export REF_CACHE=/lscratch/$SLURM_JOB_ID/hts-refcache
        
        samtools view -b -r {wildcards.rg_id} {input} > {output}
        """
    message:
        'echo {read_group_IDs}'
        'echo {output}'

rule align:
    input:
        'temp/lane_bam/{sample}.ID{rg_id}.bam'
    output:
        'temp/realigned/{sample}.ID{rg_id}.realigned.bam'
    threads: 8 
    params:
        ref = "/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta"
    run:
        import subprocess
        rg = build_RG(str(input))
        call = 'module load samtools/1.4; \
                module load bwa/0.7.12; \
                samtools collate -uOn 128 ' + str(input) + ' /lscratch/$SLURM_JOB_ID/TMP_' + str(wildcards.sample) + '_ID' + str(wildcards.rg_id) + ' | \
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
        'bam/{sample}.realigned.bam'
    threads: 4
    shell:
        """
        module load picard/2.9.2
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

