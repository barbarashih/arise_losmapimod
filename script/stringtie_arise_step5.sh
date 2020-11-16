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
TOOL_STRINGTIE='private_modules/stringtie-1.3.4c.Linux_x86_64/stringtie'
GTF_REF='genome_builds/Homo_sapiens.GRCh38.92.gtf'
IN_DIR='bam'
OUT_DIR='gtf_round2'
GTF_MERGED='gtf_merged/merged_gtf.gtf'
mkdir -p $OUT_DIR
GTF_FP="/exports/cmvm/eddie/eb/groups/freeman_mabbott_arise/ARISE/file_lst/mergelist.txt"


${TOOL_STRINGTIE} --merge -p 2 -e -G $GTF_REF -o $GTF_MERGED $GTF_FP

