#!/bin/bash
# Grid Engine options
#$ -N stringtie_makegtf
#$ -cwd
#$ -pe sharedmem 2
#$ -l h_vmem=5G
#$ -l h_rt=40:00:0
# Initialise the modules framework
. /etc/profile.d/modules.sh
# This script extract gtf from each aligned .bam
# This script uses the merged gtf to re-extract read count from each bam

SAMPLE_LIST=sample_annotation/ARISE_all.txt
CURRENT_SAMPLE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
TOOL_STRINGTIE='private_modules/stringtie-1.3.4c.Linux_x86_64/stringtie'
GTF_REF='genome_builds/Homo_sapiens.GRCh38.92.gtf'
IN_BAM="bam/${CURRENT_SAMPLE}.bam"
OUT_DIR="ballgown/${CURRENT_SAMPLE}"
GTF_MERGED='gtf_merged/merged_gtf.gtf'
mkdir -p $OUT_DIR
OUT_BALLGOWN="$OUT_DIR/${CURRENT_SAMPLE}.gtf"

${TOOL_STRINGTIE} -e -B -p 8 -G $GTF_MERGED -o $OUT_BALLGOWN $IN_BAM

