##############################################################################
# 				Project2: QIAGEN small RNA library prep kit 				               # 
# miRNA-seq library kit comparison using PBMC-isolated miRNA from		      	 # 
# Field-infected (n=4) and non-infected cattle (n=4)                         #
#     --- Linux bioinformatics workflow for pre-processing of data ---       #
##############################################################################

# Authors: Sarah Faherty O'Donnell, Carolina N. Correia.
# Last updated on: 20/01/2021

#############################################
# Generating md5 files in BYU supercomputer #
#############################################
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

# Required software is FastQC v0.11.9, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p /home/workspace/sfahertyodonnell/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/pre-filtering/Project2
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o /home/workspace/sfahertyodonnell/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/pre-filtering/Project2 \
--noextract --nogroup -t 10 \
/home/workspace/sfahertyodonnell/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project2/QIAseq3-1_S9_L003_R1_001.fastq.gz

# Create a bash script to perform FastQC quality check on all fastq.gz files in the various project directories using *:
for file in `find /home/workspace/sfahertyodonnell/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project2/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/sfahertyodonnell/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/pre-filtering/Project2/ $file" \
>> fastqc.sh; done

# Check the number of files in fastqc.sh is correct
# Given we have 8 Project2 libraries across 2 lanes, wc should be 16 (+ 1 undetermined files = 17)
wc -l fastqc.sh

# Run script on Rodeo
chmod 755 fastqc.sh
nohup ./fastqc.sh > fastqc.sh.nohup & 

# Deleted all the HTML files:
rm -r *.html

# Check if all the files were processed:
for file in `ls fastqc.sh.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Transfer compressed folders to personal laptop via SCP
# and process output files with fastqcr R package:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/sfahertyodonnell/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/pre-filtering/Project2/tmp .

# Remove temporary folder and its files:
rm -r tmp








#############################################
# Trimming of adapter sequence within reads #
#############################################

# Required software is cutadapt (version 1.10).
# Consult manual for details:
# https://cutadapt.readthedocs.io/en/stable/guide.html
# Check what version of cutadapt we have on Rodeo 
# Verion as of the 18/01/2018 is Cutadapt v1.14
cutadapt --v 

# Check if this is the latest version online
# If not, do we require the newer version for this analysis?
# Check the changes between the two versions and make decision

# Create and enter working directory:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2
cd !$

#############################################
# Use cutadapt to trim one file in Project2 #
#############################################
# Path to location of project 2 (QIAseq) raw data is \
# ~/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project2/
#############################################
# Consult QIAseq manual for details on trimming:
# file:///C:/Users/Sarah/Downloads/HB-2157-002_1103418%20HB%20QIAseq%20miRNA%20Lib%20Kit%201116%20WW.pdf
# 3' adapter sequence for QIAseq kit is (5’-3’) AACTGTAGGCACCATCAAT
# Found at this link:
# https://www.qiagen.com/kr/resources/faq?id=f12b85b4-df4f-43b5-9e82-a4fd0ddbdcc0&lang=en&Print=1
# -O 10 is the minimum overlap length of specified adapter seq and adapter seq \
# in reads
# --discard-untrimmed
# -m 17 will discard reads less than 17bp
# This resulted in the following summary in test_trim: 
#Total reads processed:              11,386,093
#Reads with adapters:                 9,870,860 (86.7%)
#Reads that were too short:             891,212 (7.8%)
#Reads written (passing filters):     8,979,648 (78.9%)


cutadapt -a AACTGTAGGCACCATCAAT -O 10 --discard-untrimmed -m 17 \
 -o test_trim.fastq.gz \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project2/QIAseq3-1_S9_L003_R1_001.fastq.gz

# If happy with the trimming parameters, remove the test_trim file
rm test_trim.fastq.gz

# Create bash script to trim the Illumina RNA 3’ Adapter (RA3) of each
# FASTQ file while keeping the sequencing lane information.
# Make sure to include all parameters of the code above
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/raw_data/Unaligned/Project2 \
-name *fastq.gz`; \
do outfile=`basename $file | perl -p -e 's/.fastq.gz//'`; \
echo "cutadapt -a AACTGTAGGCACCATCAAT -O 10 \
--discard-untrimmed -m 17 -o ./${outfile}_trim.fastq.gz $file" \
>> cutadapt.sh; \
done

# Run scripts on Rodeo:
# cutadapt.sh.nohup is where the output summaries from cutadapt get sent
chmod 755 cutadapt.sh
nohup ./cutadapt.sh > cutadapt.sh.nohup &


# Check status of trimming:
tail cutadapt.sh.nohup

# Check if all files were processed:
grep -c 'Finished' cutadapt.sh.nohup

# Generate a master file containing cutadapt trimming stats results:
# I only have one .nohup file in this project so I don't need a for loop
# grep -oP uses perl reg expression for grep to get the name of the files
# this command will get the name of each file starting with QIA
# and finishing with trim.fastq.gz
# . in perl means anything
# + in perl means for one or more characters
# because its not a for loop we don't need to append to the .txt file

grep -oP "QIA.+trim.fastq.gz" cutadapt.sh.nohup \
> $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/filename.txt \
grep "Total reads processed\:" cutadapt.sh.nohup | \
perl -p -e 's/\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)/$1/' \
> $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/processed.txt \
grep "Reads with adapters\:" cutadapt.sh.nohup | \
perl -p -e 's/\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)\s.*/$1/' \
> $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/reads_with_adapters.txt \
grep "Reads that were too short\:" cutadapt.sh.nohup | \
perl -p -e 's/\w*\s\w*\s\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)\s.*/$1/' \
> $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/short.txt \
grep "Reads written (passing filters)\:" cutadapt.sh.nohup | \
perl -p -e 's/\w*\s\w*\s.\w*\s\w*.\:\s*(\d*.\d*.\d*)\s.*/$1/' \
> $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/trimmed.txt \
paste filename.txt processed.txt reads_with_adapters.txt short.txt trimmed.txt \
> $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/trimmed_stats.txt
echo -e "Sample\tTotal reads processed\tReads with \
adapters\tReads that were too short\tReads written (passing filters)\t" \
> $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/headers.txt \
| cat headers.txt trimmed_stats.txt > trimming_stats.txt

# Check the correct number of files are in the new txt files
# should be 16 for each apart from trimming_stats.txt which
# should have 17 because of the header line
wc -l filename.txt
wc -l processed.txt
wc -l reads_with_adapters.txt
wc -l short.txt
wc -l trimmed.txt
wc -l trimmed_stats.txt

rm -r filename.txt processed.txt reads_with_adapters.txt \
short.txt trimmed.txt trimmed_stats.txt headers.txt


###############################################
# FastQC quality check of trimmed FASTQ files #
###############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter working directory:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/post_filtering/Project2
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/post_filtering/Project2 \
--noextract --nogroup -t 2 \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/QIAseq3-1_S9_L003_R1_001_trim.fastq.gz

# Create bash script to perform FastQC quality check on all trim.fastq.gz files:
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project2/ \
-name *trim.fastq.gz`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/post_filtering/Project2 $file" \
>> fastqc.sh; done

# Run all scripts on Rodeo:
chmod 755 fastqc.sh
nohup ./fastqc.sh > fastqc.sh.nohup &

# Check if correct number of files are being processed
wc -l fastqc.sh

# Check if all the files were processed:
for file in `ls fastqc.sh.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done
wc -l failed_fastqc.txt

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/post_filtering/Project2/tmp; \
done

for file in \
`find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/post_filtering/Project2/tmp \
-name summary.txt`; do more $file >> reports_post-filtering.txt; \
done

for file in \
`find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/post_filtering/Project2/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r tmp

#######################################
# Reference genome UMD3.1 preparation #
#######################################

#  ******GENOME ALREADY PREPARED WHEN ANALYSING PROJECT1 DATA*****			#
# Don't need to do this again, but this was how it was done...for reference #

# As we use miRbase as an annotation source for miRNA analyses
# We must use the same assembly version that they did
#
##gff-version 3
##date 2014-6-22
#
# Chromosomal coordinates of Bos taurus microRNAs
# microRNAs:               miRBase v21
# genome-build-id:         UMD3.1
# genome-build-accession:  NCBI_Assembly:GCA_000003055.3 
# ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/ARCHIVE/ANNOTATION_RELEASE.103/Assembled_chromosomes/seq/
# This arhived assembly didn't have the mitochondrial chromosomes

# Create and enter the reference genome directory:
mkdir -p /workspace/genomes/bostaurus/UMD3.1_NCBI/source_file
cd !$

# Download the reference genome UMD3.1 (archived annotation_release.103) from NCBI into Rodeo:
nohup wget -r -nd \
"ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/ARCHIVE/ANNOTATION_RELEASE.103/Assembled_chromosomes/seq/bt_ref_Bos_taurus_UMD_3.1_*.fa.gz" &

# Uncompress the reference genome UMD3.1 from NCBI:
# gunzip -c writes on standard output, keeping original files unchanged
# All uncompressed chromosomes will be collated in Btau_UMD3.1_multi.fa
# Individual chromosomes are then also uncompressed
gunzip -c bt_ref_Bos_taurus_UMD_3.1_*.fa.gz > Btau_UMD3.1_multi.fa
gunzip bt_ref_Bos_taurus_UMD_3.1_*.fa.gz

# Modify headers in the reference genome UMD3.1 from NCBI:
# This moves chr to the start of the header
# perl: -p assumes loop like -n but prints line also, like sed
# perl: -i edits <> files in place
# perl: -e one line of program
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
mkdir -p /home/workspace/genomes/bostaurus/UMD3.1_NCBI/annotation_file
cd !$

# Download the miRNA annotation file from miRBase (based on reference 
# genome UMD3.1 from NCBI):
wget ftp://mirbase.org/pub/mirbase/21/genomes/bta.gff3

# Convert the GFF3 annotation file from miRBase to GTF format
# using the perl script gff2gtf.pl found in scripts directory on Rodeo:
cd /home/workspace/genomes/bostaurus/UMD3.1_NCBI/annotation_file
perl $HOME/BTB_SFI_Project/WP2/Scripts/miRNA_seq_scripts/Perl_scripts/gff2gtf.pl -i \
bta.gff3 \
-o Btau_miRNA2018.gtf .
grep -P "\tpre-miRNA\t" Btau_miRNA2018.gtf >> Btau_pre-miRNA2018.gtf
grep -P "\tmiRNA\t" Btau_miRNA2018.gtf >> Btau_mature-miRNA2018.gtf


##########################################################
# Preparation of Bos taurus miRNA sequences from miRBase #
##########################################################

# Go to working directory:

mkdir /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/
cd !$
# Download the various FASTA files for mature, high confidence mature,
# precursor(hairpin), and high confidence precursor miRNA sequences
# from miRBase (version 21):
wget ftp://mirbase.org/pub/mirbase/21/hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/21/high_conf_hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/21/mature.fa.gz
wget ftp://mirbase.org/pub/mirbase/21/high_conf_mature.fa.gz

# In this directory (MiRBase_fasta)
# Uncompress the miRNA FASTA files from miRBase:
for file in \
`ls *.gz`; \
do gzip -d $file; \
done

# Combine information from miRNA annotation file and mature + high confidence
# mature FASTA sequences obtained from miRBase:
perl $HOME/BTB_SFI_Project/WP2/Scripts/miRNA_seq_scripts/Perl_scripts/miRNA_info_grepping.pl -fasta \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/mature.fa \
-gff \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/annotation_file/bta.gff3 \
-output mature_miRNA_Btaurus.txt

perl $HOME/BTB_SFI_Project/WP2/Scripts/miRNA_seq_scripts/Perl_scripts/miRNA_info_grepping.pl -fasta \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
-gff \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/annotation_file/bta.gff3 \
-output high_conf_mature_miRNA_Btaurus.txt

# Create the mature miRNA FASTA file for Bos taurus sequences only:
perl $HOME/BTB_SFI_Project/WP2/Scripts/miRNA_seq_scripts/Perl_scripts/Fasta_keep_value.pl -fasta \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/mature.fa \
-keep Bos -output bta_mature-miRNA.fa

# Create the mature miRNA FASTA file for other all species:
perl $HOME/BTB_SFI_Project/WP2/Scripts/miRNA_seq_scripts/Perl_scripts/Fasta_ignore_value.pl -fasta \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/mature.fa \
-ignore Bos -output other_mature-miRNA.fa

# Create the precursor (hairpin) miRNAs FASTA file for Bos taurus sequences only:
perl $HOME/BTB_SFI_Project/WP2/Scripts/miRNA_seq_scripts/Perl_scripts/Fasta_keep_value.pl -fasta \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/hairpin.fa \
-keep Bos -output bta_hairpin-miRNA.fa

######################
# Following analyses #
######################

# Continue pipeline to generate counts per miRNA via two different methods,
# consult the appropriate pipelines:

# Method 1: miRdeep2 software,
# see pipeline "BioValidation-miRNA_miRdeep2.sh"