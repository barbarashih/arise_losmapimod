#!/bin/bash
# Grid Engine options
#$ -N stringtie_preDE
#$ -cwd
#$ -pe sharedmem 1
#$ -l h_vmem=10G
#$ -l h_rt=40:00:00
# Initialise the modules framework

PREDE_INPUT='file_lst/preDe_input_lst.txt'
IN_GTF_DIR='gtf_round2'
STRINGTIE_PREDE='private_modules/stringtie-1.3.4c.Linux_x86_64/prepDE.py'
mkdir -p count_matrix

[ -e $PREDE_INPUT ] && rm $PREDE_INPUT
for file in $IN_GTF_DIR/*.gtf; do file_base=${file##*/}; echo -e ${file_base%.gtf} '\t' ${file} >> $PREDE_INPUT; done
python $STRINGTIE_PREDE -i $PREDE_INPUT -g count_matrix/gene_count_matrix.csv -t count_matrix/transcript_count_matrix.csv
python $STRINGTIE_PREDE -g count_matrix/gene_count_matrix.csv -t count_matrix/transcript_count_matrix_ballgown.csv
python $STRINGTIE_PREDE -i $PREDE_INPUT -g count_matrix/gene_count_matrix_cluster.csv -t count_matrix/transcript_count_matrix_cluster.csv -c
python $STRINGTIE_PREDE -g count_matrix/gene_count_matrix_cluster.csv -t count_matrix/transcript_count_matrix_cluster_ballgown.csv -c

