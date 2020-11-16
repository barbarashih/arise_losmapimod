#!/bin/bash
# Grid Engine options
#$ -N hisat2_alignment
#$ -cwd
#$ -pe sharedmem 1
#$ -l h_vmem=10G
#$ -l h_rt=40:00:00
# Initialise the modules framework
. /etc/profile.d/modules.sh
module add igmm/apps/samtools/1.6
SAMPLE_LIST="sample_annotation/ARISE_all.txt"
CURRENT_SAMPLE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
#TOOL_PICARD='private_modules/picard-tools.2.17.10/picard.jar'
TOOL_HISAT2='private_modules/hisat2-2.1.0/hisat2'
GENOME_FA='genome_builds/Homo_sapiens.GRCh38.dna.toplevel.fa'
KNOWN_SPLICE_SITES='genome_builds/splicesites.txt'
FASTQ_DIR='fastq'
TMP_DIR="tmp/${CURRENT_SAMPLE}"
OUT_DIR=$(pwd)
mkdir -p $TMP_DIR
mkdir -p "${OUT_DIR}/bam"
mkdir -p "${OUT_DIR}/novel_splice"
#mkdir -p "${OUT_DIR}/picard_summary"
mkdir -p "${OUT_DIR}/unmapped"
mkdir -p "${OUT_DIR}/flagstat"
mkdir -p "${OUT_DIR}/flagstat_filtered"
mkdir -p "${OUT_DIR}/hisat2_summary"

TMP_SAM="${TMP_DIR}/${CURRENT_SAMPLE}.sam"
TMP_BAM="${TMP_DIR}/${CURRENT_SAMPLE}.bam"
TMP_BAM_SORTED="${TMP_DIR}/${CURRENT_SAMPLE}_sorted.bam"
IN_FASTQ1="${FASTQ_DIR}/${CURRENT_SAMPLE}_R1_001.fastq.gz"
IN_FASTQ2="${FASTQ_DIR}/${CURRENT_SAMPLE}_R2_001.fastq.gz"
OUT_BAM="${OUT_DIR}/bam/${CURRENT_SAMPLE}.bam"
#OUT_SPLICE="${OUT_DIR}/novel_splice/${CURRENT_SAMPLE}.bed"
OUT_UNMAPPED="${OUT_DIR}/unmapped/${CURRENT_SAMPLE}.bam"
OUT_PICARD_SUMMARY="${OUT_DIR}/picard_summary/${CURRENT_SAMPLE}.txt"
OUT_FLAGSTAT="${OUT_DIR}/flagstat/${CURRENT_SAMPLE}.txt"
OUT_FLAGSTAT_filtered="${OUT_DIR}/flagstat_filtered/${CURRENT_SAMPLE}.txt"
OUT_HISAT2_SUMMARY="${OUT_DIR}/hisat2_summary/${CURRENT_SAMPLE}.txt"

${TOOL_HISAT2} -x $GENOME_FA -1 $IN_FASTQ1 -2 $IN_FASTQ2 -S $TMP_SAM --known-splicesite-infile $KNOWN_SPLICE_SITES --dta
samtools view -u $TMP_SAM | samtools sort -o $TMP_BAM
samtools flagstat $TMP_BAM > $OUT_FLAGSTAT
samtools view -b -f 4 $TMP_BAM > $OUT_UNMAPPED
samtools view -b -f 1 -F 12 $TMP_BAM > $OUT_BAM
#java -jar $TOOL_PICARD CleanSam INPUT=$TMP_BAM OUTPUT=$OUT_BAM
samtools index $OUT_BAM
samtools flagstat $OUT_BAM > $OUT_FLAGSTAT_filtered

#java -jar $TOOL_PICARD CleanSam INPUT=$TMP_BAM OUTPUT=$OUT_BAM
#java -jar $TOOL_PICARD CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=$GENOME_FA I=$TMP_SAM O=$OUT_PICARD_SUMMARY
rm -r $TMP_DIR
