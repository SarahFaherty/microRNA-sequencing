##############################################################################
# miRNA-seq library kit comparison using PBMC-isolated miRNA from			 # 
# Field-infected (n=4) and non-infected cattle (n=4)  						 # 
#     --- Linux bioinformatics workflow for pre-processing of data ---       #
##############################################################################
# Authors: Sarah Faherty O'Donnell
# Zenodo DOI badge:
# Version 
# Last updated on: 08/01/2018

###############################
# Generating md5 files in BYU supercomputer #
##############################
cd /fslgroup/fslg_dnasc/compute/170901_D00723_0220_ACBGBTANXX/Unaligned/Project/ 
for file in `find /fslgroup/fslg_dnasc/compute/170901_D00723_0220_ACBGBTANXX/Unaligned/ \
-name '*.fastq.gz'`; \
do md5sum $file >> md5byu.txt; \
done 

################################
# Download and files check sum #
################################
# Create and enter the data storage directory:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data
cd !$

# Download data from MSU server into Rodeo server directory made above
screen -D -R BYU_miRNA_download2
scp -r fslcollab164@scp.fsl.byu.edu:/fslhome/fslcollab164/fsl_groups/fslg_dnasc/compute/170901_D00723_0220_ACBGBTANXX/Unaligned/ .


# Change directory permissions:
chmod -R 755 $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data

# grep relevant lines from the md5 file
grep QIA md5byu.txt > md5byu_Project2.txt
wc -l md5byu_Project2.txt
cp md5byu_Project2.txt Project2

# First need to correct path from byu md5 files to contain only 
# file names
# Use awk to get columns from files and pipe to perl for editing 
# Searching and substituting using '///'
# piped to perl again to get rid of a forward slash infornt of file name
# piped column to a new file 
# used paste to combine columns and piped to a new file
# Repeat this for each project
cd $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project2
awk '{print $2}' md5byu_Project2.txt | perl -p -e 's/fslgroup.*(QIA.*\.gz)$/$1/' | perl -p -e s'/\///' > byumd5_filenames.txt 
paste <(awk '{print $1}' md5byu_Project2.txt) <(awk '{print $1}' byumd5_filenames.txt) > md5byu_Project2_corrected.txt 

# Check md5sum files  
md5sum -c md5byu_Project2_corrected.txt >> md5_UCD_Project2.txt

# Check that all files passed the check:
grep -c 'OK' md5_UCD_Project2.txt

# Change directory permissions to read and execute only: ???
chmod -R 555 $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data


###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/miRNAseqValidation/quality_check/pre-filtering
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/pre-filtering \
--noextract --nogroup -t 2 \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project1/NEXTflex-v10_S8_L001_R1_001.fastq.gz

# Create a bash script to perform FastQC quality check on all fastq.gz files in the various project directories using *:
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project*/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/pre-filtering $file" \
>> fastqc.sh; done

# Split and run all scripts on Stampede:
split -d -l 24 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/pre-filtering/tmp; \
done

for file in \
`find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; \
done

for file in \
`find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -rf $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/pre-filtering/tmp

#############################################
# Trimming of adapter sequence within reads #
#############################################

# Required software is cutadapt (version 1.10).
# Consult manual for details:
# https://cutadapt.readthedocs.io/en/stable/guide.html
# download and install TrimGalore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz

# Create and enter working directory:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/miRNAseqValidation/quality_check/trimming
cd !$

############################
# Path to location of project 1 (NEXTFlex) raw data is \
# ~/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project1
###########################

# Create bash script to trim the Illumina RNA 3’ Adapter (RA3) of each
# FASTQ file while keeping the sequencing lane information:
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project1 \
-name *fastq.gz`; \
do outfile=`basename $file | perl -p -e 's/.fastq.gz//'`; \
echo "cutadapt -a TGGAATTCTCGGGTGCCAAGG -O 10 --match-read-wildcards \
--discard-untrimmed -m 17 -o ./${outfile}_trim.fastq.gz $file" \
>> cutadapt.sh; \
done

# Split and run scripts on Stampede:
split -d -l 24 cutadapt.sh cutadapt.sh.
for script in `ls cutadapt.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done
# I stopped here on the 20.11.2017
# 

# Check if all files were processed:
grep -c 'Finished' cutadapt.sh.00.nohup
grep -c 'Finished' cutadapt.sh.01.nohup

# Generate a master file containing cutadapt trimming stats results:
for file in `ls $HOME/scratch/miRNAseqValidation/fastq_sequence/*.nohup`; \
do grep -oP "E\d+_\w+_L00(5|6)\_R1_001.fastq.gz" $file \
>> $HOME/scratch/miRNAseqValidation/fastq_sequence/filename.txt; \
grep "Total reads processed\:" $file | \
perl -p -e 's/\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)/$1/' \
>> $HOME/scratch/miRNAseqValidation/fastq_sequence/processed.txt; \
grep "Reads with adapters\:" $file | \
perl -p -e 's/\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)\s.*/$1/' \
>> $HOME/scratch/miRNAseqValidation/fastq_sequence/reads_with_adapters.txt; \
grep "Reads that were too short\:" $file | \
perl -p -e 's/\w*\s\w*\s\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)\s.*/$1/' \
>> $HOME/scratch/miRNAseqValidation/fastq_sequence/short.txt; \
grep "Reads written (passing filters)\:" $file | \
perl -p -e 's/\w*\s\w*\s.\w*\s\w*.\:\s*(\d*.\d*.\d*)\s.*/$1/' \
>> $HOME/scratch/miRNAseqValidation/fastq_sequence/trimmed.txt; \
paste filename.txt processed.txt reads_with_adapters.txt short.txt trimmed.txt \
> $HOME/scratch/miRNAseqValidation/fastq_sequence/trimmed_stats.txt; \
done

echo -e "Sample\tTotal reads processed\tReads with \
adapters\tReads that were too short\tReads written (passing filters)\t" \
| cat - $HOME/scratch/miRNAseqValidation/fastq_sequence/trimmed_stats.txt \
> $HOME/scratch/miRNAseqValidation/fastq_sequence/trimming_stats.txt

rm -r filename.txt processed.txt reads_with_adapters.txt \
short.txt trimmed.txt trimmed_stats.txt

###############################################
# FastQC quality check of trimmed FASTQ files #
###############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter working directory:
mkdir $HOME/scratch/miRNAseqValidation/quality_check/post_filtering
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/miRNAseqValidation/quality_check/post_filtering \
--noextract --nogroup -t 2 \
$HOME/scratch/miRNAseqValidation/fastq_sequence/E10_CGTACG_L005_R1_001_trim.fastq.gz

# Create bash script to perform FastQC quality check on all trim.fastq.gz files:
for file in `find $HOME/scratch/miRNAseqValidation/fastq_sequence \
-name *trim.fastq.gz`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/scratch/miRNAseqValidation/quality_check/post_filtering $file" \
>> fastqc.sh; done

# Split and run all scripts on Stampede:
split -d -l 24 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/miRNAseqValidation/quality_check/post_filtering/tmp; \
done

for file in \
`find $HOME/scratch/miRNAseqValidation/quality_check/post_filtering/tmp \
-name summary.txt`; do more $file >> reports_post-filtering.txt; \
done

for file in \
`find $HOME/scratch/miRNAseqValidation/quality_check/post_filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r tmp

#######################################
# Reference genome UMD3.1 preparation #
#######################################

# Create and enter the reference genome directory:
mkdir -p /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/source_file
cd !$

# Download the reference genome UMD3.1 from NCBI into Stampede:
nohup wget -r -nd \
"ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/Assembled_chromosomes/seq/bt_ref_Bos_taurus_UMD_3.1_*.fa.gz" &

# Uncompress the reference genome UMD3.1 from NCBI:
gunzip -c bt_ref_Bos_taurus_UMD_3.1_*.fa.gz > Btau_UMD3.1_multi.fa
gunzip bt_ref_Bos_taurus_UMD_3.1_*.fa.gz

# Modify headers in the reference genome UMD3.1 from NCBI:
for file in `ls *.fa`; \
do perl -p -i -e \
's/^>(.*)(Bos taurus breed Hereford chromosome )(.{1,2})(\,.*)$/>chr$3 $1$2$3$4/' \
$file; \
done

for file in `ls *.fa`; \
do perl -p -i -e 's/^>(.*)(Bos taurus mitochondrion)(\,.*)$/>chrMT $1$2$3/' \
$file; \
done

# Create and enter the reference genome annotation directory:
mkdir -p /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file
cd !$

# Download the miRNA annotation file from miRBase (based on reference 
# genome UMD3.1 from NCBI):
wget ftp://mirbase.org/pub/mirbase/21/genomes/bta.gff3

# Convert the GFF3 annotation file from miRBase to GTF format:
perl /home/nnalpas/SVN/gff2gtf.pl -i \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file/bta.gff3 \
-o /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file/Btau_miRNA2016.gtf
grep -P "\tpre-miRNA\t" Btau_miRNA2016.gtf >> Btau_pre-miRNA2016.gtf
grep -P "\tmiRNA\t" Btau_miRNA2016.gtf >> Btau_mature-miRNA2016.gtf

##########################################################
# Preparation of Bos taurus miRNA sequences from miRBase #
##########################################################

# Go to working directory:
cd /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta
# Download the various FASTA files for mature, high confidence mature,
# precursor(hairpin), and high confidence precursor miRNA sequences
# from miRBase (version 21):
wget ftp://mirbase.org/pub/mirbase/21/hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/21/high_conf_hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/21/mature.fa.gz
wget ftp://mirbase.org/pub/mirbase/21/high_conf_mature.fa.gz

# Uncompress the miRNA FASTA files from miRBase:
for file in \
`ls /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta*.gz`; \
do gzip -d $file; \
done

# Combine information from miRNA annotation file and mature + high confidence
# mature FASTA sequences obtained from miRBase:
perl /home/nnalpas/Scripts/miRNA_info_grepping.pl -fasta \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/mature.fa \
-gff \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file/bta.gff3 \
-output mature_miRNA_Btaurus.txt

perl /home/nnalpas/Scripts/miRNA_info_grepping.pl -fasta \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
-gff \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file/bta.gff3 \
-output high_conf_mature_miRNA_Btaurus.txt

# Create the mature miRNA FASTA file for Bos taurus sequences only:
perl /home/nnalpas/SVN/Fasta_keep_value.pl -fasta \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/mature.fa \
-keep Bos -output bta_mature-miRNA.fa

# Create the mature miRNA FASTA file for other all species:
perl /home/nnalpas/SVN/Fasta_ignore_value.pl -fasta \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/mature.fa \
-ignore Bos -output other_mature-miRNA.fa

# Create the precursor (hairpin) miRNAs FASTA file for Bos taurus sequences only:
perl /home/nnalpas/SVN/Fasta_keep_value.pl -fasta \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/hairpin.fa \
-keep Bos -output bta_hairpin-miRNA.fa

######################
# Following analyses #
######################

# Continue pipeline to generate counts per miRNA via two different methods,
# consult the appropriate pipelines:

# Method 1: Novoalign-featureCounts softwares,
# see pipeline "BioValidation-miRNA_Novoalign-featureCounts.sh"

# Method 2: miRdeep2 software,
# see pipeline "BioValidation-miRNA_miRdeep2.sh"
