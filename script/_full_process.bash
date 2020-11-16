#### Requirement
# samtools/1.6
# picard-tools.2.17.10/picard.jar
# hisat2-2.1.0/hisat2
# stringtie-1.3.4c.Linux_x86_64

#### Download tools
mkdir -p genome_builds private_modules
mkdir private_modules/picard-tools.2.17.10
wget https://github.com/broadinstitute/picard/releases/download/2.17.10/picard.jar -P private_modules/picard-tools.2.17.10

# Download hisat from https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download and save files in private_modules
# Such that there is the 'private_modules/hisat2-2.1.0/hisat2'
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4c.Linux_x86_64.tar.gz -P private_modules
wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download -O private_modules/hisat2-2.1.0-Linux_x86_64.zip

cd private_modules
tar -xzf stringtie-1.3.4c.Linux_x86_64.tar.gz
unzip hisat2-2.1.0-Linux_x86_64.zip
cd ..
wget http://ccb.jhu.edu/software/stringtie/dl/prepDE.py -P private_modules/stringtie-1.3.4c.Linux_x86_64

#### Download reference genomes
wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz -P genome_builds
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -P genome_builds

cd genome_builds
gunzip Homo_sapiens.GRCh38.92.gtf.gz
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz

cd ..


#### Align with hisat
# build index
private_modules/hisat2-2.1.0/hisat2-build genome_builds/Homo_sapiens.GRCh38.dna.toplevel.fa genome_builds/Homo_sapiens.GRCh38.dna.toplevel.fa
# make splice site infile
hisat2_extract_splice_sites.py genome_builds/Homo_sapiens.GRCh38.92.gtf > genome_builds/splicesites.txt
# alignment
qsub -N hisat -t 1-126 script/hisat_alignment_arise.sh

#### Assemble
qsub -hold_jid hisat -N  stringtie1 -t 1-126 script/stringtie_arise_step1.sh
qsub -hold_jid stringtie1 -N stringtie2 script/stringtie_arise_step2.sh
qsub -hold_jid stringtie2 -N stringtie3 -t 1-126 script/stringtie_arise_step3.sh
qsub -hold_jid stringtie3 -N stringtie4 script/stringtie_arise_step4.sh
qsub -hold_jid stringtie4 -N stringtie4.1 -t 1-126 script/stringtie_makegtf_arise_step1.sh
qsub -hold_jid stringtie4.1 -N stringtie4.2 script/stringtie_makegtf_arise_step2.sh
qsub -hold_jid stringtie4.2 -N stringtie4.3 -t 1-126 script/stringtie_makegtf_arise_step3.sh
qsub -hold_jid stringtie4.3 -N stringtie5 script/stringtie_arise_step5.sh
qsub -hold_jid stringtie5 -N stringtie6 -t 1-126 script/stringtie_arise_step6.sh

##### Stats analysis
Rscript script/r_rnaseq_stats.r