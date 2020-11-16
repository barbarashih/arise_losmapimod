#!/bin/bash
# Grid Engine options
#$ -N stringtie_merge_gtf
#$ -cwd
#$ -pe sharedmem 1
#$ -l h_vmem=1G
#$ -l h_rt=40:00:0
# Initialise the modules framework
. /etc/profile.d/modules.sh
# This script uses qsub to sort and index all bam files in the SRR_list
TOOL_STRINGTIE='private_modules/stringtie-1.3.4c.Linux_x86_64/stringtie'
GTF_REF='genome_builds/Homo_sapiens.GRCh38.92.gtf'
TMP_DIR="tmp/stringtie"
TMP_GTF_LST="file_lst/gtf_lst.txt"
IN_DIR='gtf_round1'
OUT_DIR='gtf_merged'
mkdir -p $TMP_DIR
mkdir -p $OUT_DIR $(dirname $TMP_GTF_LST)
OUT_GTF="${OUT_DIR}/merged_gtf.gtf"

find ${IN_DIR}/*.gtf > $TMP_GTF_LST

#### STEP 2: Merge GTF for individual sample (stringtie) ####
# Merge gtf files from each sample using the TEMP_GTF_LIST created in the loop in STEP 1.
${TOOL_STRINGTIE} --merge -G $GTF_REF -o $OUT_GTF -T 1 -i $TMP_GTF_LST
