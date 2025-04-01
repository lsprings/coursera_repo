#!/bin/bash
#SBATCH -J mapping         #job name
#SBATCH -o mapping_o%j     #output file name (%j expands to jobID)
#SBATCH -e mapping_e%j   #error file
#SBATCH -n 24                   #total number of mpi tasks requested
#SBATCH -N 1                    #total number of nodes requested
#SBATCH -p skx               #queue (partition) -- normal, development, etc.
#SBATCH -t 14:30:00              #run time (hh:mm:ss)
#SBATCH -A TG-BIO240048
#SBATCH --mail-user=lsprings@utexas.edu #your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

date +"%c starting pipeline" >> ${1}_timelog.txt
module unload xalt
cd combined_datasets_3/combined_fastqs/
Ref='/work2/03177/lcs2347/stampede2/hg19/hg19.fa'
dbSNP='/work2/03177/lcs2347/stampede2/hg19/Homo_sapiens_assembly19.dbsnp138.myChNames.vcf.gz'

#cat /scratch/03177/lcs2347/batch3/trimmed_plus_qcs/${1}*1P.fastq.gz /scratch/03177/lcs2347/batch2/trimmed_plus_qcs/${1}*1P.fastq.gz > ${1}_combo_1P.fastq.gz
#cat /scratch/03177/lcs2347/batch3/trimmed_plus_qcs/${1}*2P.fastq.gz /scratch/03177/lcs2347/batch2/trimmed_plus_qcs/${1}*2P.fastq.gz > ${1}_combo_2P.fastq.gz
#### set to add read groups later ######
header=$(zcat ${1}_combo_1P* | head -n 1)
id=$(echo $header | head -n 1 | cut -f 1,4 -d":" | sed 's/@//')
pu=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//')

#apptainer exec /work2/03177/lcs2347/stampede3/containers/bwa_latest.sif bwa mem -t 24 -R $(echo "@RG\tID:${1}\tPL:illumina\tSM:${1}") $Ref /scratch/03177/lcs2347/combined_datasets_3/combined_fastqs/${1}_combo_1P.fastq.gz /scratch/03177/lcs2347/combined_datasets_3/combined_fastqs/${1}_combo_2P.fastq.gz > ${1}_combined_geno.sam
cd ../mapped_bams
apptainer exec /work2/03177/lcs2347/stampede3/containers/samtools_latest.sif samtools view -b -o ${1}_combined_geno.bam ${1}_combined_geno.sam
apptainer exec /work2/03177/lcs2347/stampede3/containers/samtools_latest.sif samtools sort -o ${1}_combined_geno_sorted.bam ${1}_combined_geno.bam
apptainer exec /work2/03177/lcs2347/stampede3/containers/samtools_latest.sif samtools index ${1}_combined_geno_sorted.bam
#mv *bam* ../mapped_bams
#mv *sam ../mapped_bams
date +"%c mapping complete" >> ../../${1}_timelog.txt

# MARK DUPLICATES and generate alignment stats 
cd ../mapped_bams
apptainer exec /work2/03177/lcs2347/stampede3/containers/gatk_latest.sif gatk MarkDuplicates -I ${1}_combined_geno_sorted.bam -O ${1}_combined_geno_sorted.markeddups.bam -M ${1}_combined_geno_sorted.markeddups.metrics
apptainer exec /work2/03177/lcs2347/stampede3/containers/samtools_latest.sif samtools flagstat ${1}_combined_geno_sorted.markeddups.bam > ${1}_combined_geno_sorted.markeddups.bam.flagstat
date +"%c markdups complete" >> ../../${1}_timelog.txt

#### ADD READGROUPS to replace read groups to add platform info ######
apptainer exec /work2/03177/lcs2347/stampede3/containers/gatk_latest.sif gatk AddOrReplaceReadGroups -I ${1}_combined_geno_sorted.markeddups.bam -O ${1}_combined_geno_sorted.markeddups.addRG.bam -RGID ${id} -RGLB combined -RGPL illumina -RGPU ${pu} -RGSM ${1}
date +"%c RG added" >> ../../${1}_timelog.txt

#### RECALIBRATING BAMS ######
apptainer exec /work2/03177/lcs2347/stampede3/containers/gatk_latest.sif gatk BaseRecalibrator -I ${1}_combined_geno_sorted.markeddups.addRG.bam -R ${Ref} --known-sites ${dbSNP} -O ${1}_combined_geno_sorted.markeddups.recal.table
apptainer exec /work2/03177/lcs2347/stampede3/containers/gatk_latest.sif gatk ApplyBQSR -R ${Ref} -I ${1}_combined_geno_sorted.markeddups.addRG.bam --bqsr-recal-file ${1}_combined_geno_sorted.markeddups.recal.table -O ${1}_combined_geno_sorted.markeddups.addRG.BQSR.bam
date +"%c BAM recalibrated" >> ../../${1}_timelog.txt

#### BASE CALLING ######
apptainer exec /work2/03177/lcs2347/stampede3/containers/gatk_latest.sif gatk --java-options "-Xmx36g" HaplotypeCaller -R ${Ref} -I ${1}_combined_geno_sorted.markeddups.addRG.BQSR.bam -O ${1}_combined_geno_sorted.markeddups.addRG.BQSR.g.vcf.gz -ERC GVCF
date +"%c haplotypes called" >> ../../${1}_timelog.txt

# with the combined fastqs, some of these took like 14 hrs
# Run this with 'sbatch JOB_NAME.sh SAMPLENAME', where SAMPLENAME refers to the alphanumeric code of each sample. Each individual will be run as its own job

