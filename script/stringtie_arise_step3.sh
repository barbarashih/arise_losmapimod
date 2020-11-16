#!/bin/bash
# Grid Engine options
#$ -N stringtie_makegtf
#$ -cwd
#$ -pe sharedmem 2
#$ -l h_vmem=500M
#$ -l h_rt=40:00:0
# Initialise the modules framework
. /etc/profile.d/modules.sh
# This script extract gtf from each aligned .bam
# This script uses the merged gtf to re-extract read count from each bam

SAMPLE_LIST=sample_annotation/ARISE_all.txt
CURRENT_SAMPLE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
TOOL_STRINGTIE='private_modules/stringtie-1.3.4c.Linux_x86_64/stringtie'
GTF_REF='genome_builds/Homo_sapiens.GRCh38.92.gtf'
IN_DIR='bam'
OUT_DIR='gtf_round2'
GTF_MERGED='gtf_merged/merged_gtf.gtf'
mkdir -p $OUT_DIR
IN_BAM="${IN_DIR}/${CURRENT_SAMPLE}.bam"
OUT_GTF="${OUT_DIR}/${CURRENT_SAMPLE}.gtf"

${TOOL_STRINGTIE} $IN_BAM -p 2 -G $GTF_REF -e -o $OUT_GTF


