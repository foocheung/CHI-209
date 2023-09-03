from os import listdir
import glob
import os

configfile: "config.yaml"
#prin://cdn.sstatic.net/Img/teams/teams-illo-free-sidebar-promo.svg?v=47faa659a05et (config['samples'])
#SAMPLES = set(samples)
print (config['samples'])


rule all:
    input:
        expand("outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam", sample = config['samples']),
        expand("outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.bai", sample = config['samples']),
        expand("outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.dedup", sample = config['samples']),
        expand("outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.dedup.count", sample = config['samples'])
             
rule pass1:
        input:
            R1L1 = '/hpcdata/chi/PROJECTS/2022_CHI_PROPOSALS/RAMASWAMI_KSHV_CHI-209/RNASEQ/230821_VH00286_63_AAAVGTVHV/RAMA/FASTQ/{sample}_R1_001.fastq.gz', # may need adjustment if your fastq file name format is different
            R2L1 = '/hpcdata/chi/PROJECTS/2022_CHI_PROPOSALS/RAMASWAMI_KSHV_CHI-209/RNASEQ/230821_VH00286_63_AAAVGTVHV/RAMA/FASTQ/{sample}_R2_001.fastq.gz',
            refdir = directory('/hpcdata/sg/sg_data/CHI/PROJECTS/HSV/RNA_SEQ/RNASEQ/DATA/genome')
        params:
            outdir = 'outdir/{sample}_seq2_output'
        output:
            "outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam"
        threads: 1 # set the maximum number of available cores
        shell:
            'rm -rf {params.outdir} &&' # be careful with this. I don't know why, but Snakemake had problems without this cleaning.
            'mkdir {params.outdir} && ' # snakemake had problems finding output files with --outFileNamePrefix, so I used this approach instead
            'cd {params.outdir} && '
            '/sysapps/cluster/software/STAR/2.5.2a-goolf-1.7.20/bin/STAR --runThreadN {threads} '
            '--genomeDir {input.refdir} '
            '--readFilesIn {input.R1L1} {input.R2L1} '
            '--readFilesCommand zcat '
            '--outSAMtype BAM SortedByCoordinate && cd ..'


#rule sort:
#        input:
#            sortbam = 'outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam'
#        params:
#            outdir = 'outdir/{sample}_seq2_output'
#        output:
#            outdir = 'outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.sort'
#        threads:1  # set the maximum number of available cores
#        shell:
#            '/sysapps/cluster/software/SAMtools/1.9-goolf-1.7.20/bin/samtools  sort -@ 6 -n {input.sortbam} -o  {output.outdir}'

rule index:
        input:
            sortbam = 'outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam'
        params:
            outdir = 'outdir/{sample}_seq2_output'
        output:
            outdir = 'outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.bai'
        threads:1  # set the maximum number of available cores
        shell:
            '/sysapps/cluster/software/SAMtools/1.9-goolf-1.7.20/bin/samtools  index {input.sortbam} '


rule dedup:
        input:
            sortbam = 'outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam', 
            sortbami= 'outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.bai'
        params:
            outdir = 'outdir/{sample}_seq2_output'
        output:
            outdir = 'outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.dedup'
        threads:1  # set the maximum number of available cores
        shell:
            '/sysapps/cluster/software/Anaconda3/2022.10/envs/umi_tools-1.1.4/bin/umi_tools dedup --umi-separator=":" -I {input.sortbam} -S  {output.outdir}'



rule count:
        input:
            sortbam='outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.dedup'
        params:
            outdir ='outdir/{sample}_seq2_output'
        output:
            outdir ='outdir/{sample}_seq2_output/Aligned.sortedByCoord.out.bam.dedup.count'
        shell:
            '/sysapps/cluster/software/subread/2.0.0/bin/featureCounts --ignoreDup  -p  -T 10 -a /hpcdata/sg/sg_data/CHI/PROJECTS/RNA_SEQ/GALINA/GENOME/data2/Homo_sapiens.GRCh38.95.gtf -s 1 -o {output.outdir} {input.sortbam} ' 

 
