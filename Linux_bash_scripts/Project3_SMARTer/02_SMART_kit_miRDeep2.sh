##############################################################################
#						Project3: SMARTer miRNA_LibraryKit_Study 			 #
#	 miRNA-seq library kit comparison using PBMC-isolated miRNA from		 # 
# 			Field-infected (n=4) and non-infected cattle (n=4)  			 #
#    --- Linux bioinformatics workflow for known and novel miRNAs  ---       #
#                     -- Method : miRDeep2 software --                       #
##############################################################################
# Authors: Sarah Faherty O'Donnell 
# Zenodo DOI badge:
# Version 
# Last updated on: 01/02/2018

########################################################################
# Merge and uncompress miRNA-seq FASTQ files to be used with miRDeep2  #
########################################################################

# Create and enter temporary working directory:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project3/tmp
cd !$

# Create bash script to uncompress and merge trim.fast.gz from lanes
# 005 and 006 for each library:
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project3 \
-name *_L005_R1_001_trim.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_L005_)/_L006_/'`; \
sample=`basename $file | perl -p -e 's/(SMART\w+\d-\d_\w\d+_)\w+\d+_\w+\d+_\d*_trim.fastq.gz/$1/'`; \
echo "zcat $file $file2 > ${sample}merged.trim.fastq" \
>> uncompress_merge.sh; \
done

# Run script:
chmod 755 uncompress_merge.sh
nohup ./uncompress_merge.sh &

#######################################
# Index reference genome using Bowtie #
#######################################

# Rodeo version: bowtie-1.1.0
# Bowtie uses the Burrows–Wheeler transform 
# consult manual for details:
# http://bowtie-bio.sourceforge.net/manual.shtml

# Create and enter the Index reference genome directory:
mkdir -p /home/workspace/genomes/bostaurus/UMD3.1_NCBI/bowtie-1.1.0_index
cd !$

# Index the reference genome UMD3.1 using Bowtie:
nohup bowtie-build \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/source_file/Btau_UMD3.1_multi.fa \
Btau_UMD3.1_multi_index &

# Test installation of Bowtie index
# using sequence from ncbi chr1 bos taurus fasta
bowtie -c Btau_UMD3.1_multi_index  AGTACTTTGACTCCTCTAACTAGGCAAGGTTTTGACTGAAA

#############################################################
# Preprocessing of miRNA-seq data using miRDeep2: mapper.pl #
#############################################################

# Required software is miRDeep2 v.2.0.0.8, consult manual/tutorial for details:
# https://www.mdc-berlin.de/8551903/en/
# Check installation of software
mapper.pl
quantifier.pl
miRDeep2.pl

# Create and enter the miRDeep2 directory for mapping work:
mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper
cd !$

# Create symbolic links to FASTQ files:
for file in \
`ls $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/quality_check/trimming/Project3/tmp/*_merged.trim.fastq`; \
do ln -s \
$file $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/`basename $file`; \
done

# Run mapper.pl in one FASTQ file to see if it's working well:
mapper.pl SMARTer4-1_S17_merged.trim.fastq -e -h -m -p \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/bowtie-1.1.0_index/Btau_UMD3.1_multi_index \
-s test_collapsed.fa -t test.arf -v

# If this works, delete the test files:
rm test.arf test_collapsed.fa 
rm bowtie.log
rm -r mapper_logs

# Create bash script to map miRNA reads to the reference genome:
for file in `ls *_merged.trim.fastq`; \
do outfile=`basename $file | perl -p -e 's/_merged\.trim\.fastq//'`; \
echo "mapper.pl $file -e -h -m -p \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/bowtie-1.1.0_index/Btau_UMD3.1_multi_index \
-s ${outfile}_collapsed.fa -t ${outfile}.arf -v" \
>> mapper.sh; \
done

# Run all scripts on Rodeo:
chmod 755 mapper.sh
nohup ./mapper.sh > mapper.sh.nohup &

# Generate a master file containing mapping statistics
for file in \
`ls $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/mapper*.nohup`; \
do grep -oP "log_\d+" $file >> ./log_id.txt; \
grep "total:" $file >> ./totals.txt; \
paste log_id.txt totals.txt > ./stats.txt; \
done

# Using awk, keep only columns of interest from stats.txt
awk '{print $1, $3, $4, $5, $6, $7}' stats.txt > stats2.txt

# adding header to stats2.txt and save mapper_stats.txt
echo -e "Log_Id\tInput reads\tTotal Mapped Reads\tTotal Unmapped Reads\tPercentage Mapped Reads\t Percentage Unmapped Reads" \
| cat - ./stats2.txt > ./mapper_stats.txt

# delete unnecessary files 
rm -r log_id.txt totals.txt stats.txt stats2.txt

# Transfer mapper_stats.txt to laptop using Winscp

################################################################
# Quantification of known miRNAs using miRDeep2: quantifier.pl #
################################################################

# Create and enter the working directory:
mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier
cd !$

# Run quantifier.pl in one FASTA file to see if it's working well:
quantifier.pl -p \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/SMARTer4-1_S17_collapsed.fa -t bta

# Create a shell script to quantify the mapped mature miRNAs:
# [You will need the mature and precursor (hairpin) miRNA FASTA files for
# Bos taurus sequences only. Please refer to the 
# BioValidation-miRNA-seq_QC_filter.sh script]
for file in \
`ls $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier/$outfile; \
cd $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier/$outfile; \
quantifier.pl -p \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $file -t bta" >> quantifier.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 1 quantifier.sh quantifier.sh.
for script in `ls quantifier.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier/high_confidence
cd !$
# Create a shell script to quantify the mapped high confidence mature and
# precursor (hairpin) miRNAs:
for file in \
`ls $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier/high_confidence/$outfile; \
cd $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier/high_confidence/$outfile; \
quantifier.pl -p \
/home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
-m /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_hairpin.fa \
-r $file -t bta" >> quantifier_high_conf.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 1 quantifier_high_conf.sh quantifier_high_conf.sh.
for script in `ls quantifier_high_conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all read counts files from regular and high confidence miRNAs for
# transfering into laptop using WinSCP:
mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/regular
cd !$
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier/SMART* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/regular/${outfile}_expressed.csv; \
done

mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/high_conf
cd !$
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/quantifier/high_confidence/SMART* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/.*\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/high_conf/${outfile}_expressed_high_conf.csv; \
done

########################################################################
# Identification of known and novel miRNAs using miRDeep2: miRDeep2.pl #
######################################################################## 

# Create and enter working directory:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep
cd !$

# Copy and modify the required fasta files since miRdeep2 software
# requires no space in headers and no characters other than acgtunACGTUN in
# the sequences:
cp /home/workspace/genomes/bostaurus/UMD3.1_NCBI/source_file/Btau_UMD3.1_multi.fa \
./Btau_UMD3.1_multi.fa
cp /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
./bta_mature-miRNA.fa
cp /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
./high_conf_mature.fa
cp /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/other_mature-miRNA.fa \
./other_mature-miRNA.fa
cp /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
./bta_hairpin-miRNA.fa
cp /home/workspace/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_hairpin.fa \
./high_conf_hairpin.fa

perl -p -i -e 's/^(>.*?)\s.*$/$1/' \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/Btau_UMD3.1_multi.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/bta_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*$)/$1_$2/' \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_conf_mature.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/other_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/bta_hairpin-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) (.*?) (.*?) (.*?) .*$/$1_$2/' \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_conf_hairpin.fa
sed -e '/^[^>]/s/[^ATGCatgc]/N/g' high_conf_hairpin.fa > high_conf_hairpin_clean.fa

# Run mirdeep.pl in one FASTA file to see if it's working well:
miRDeep2.pl $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/SMARTer4-1_S17_collapsed.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/SMARTer4-1_S17.arf  \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/bta_mature-miRNA.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/other_mature-miRNA.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/bta_hairpin-miRNA.fa -t Cow

# Create bash script for identification and quantification of known and
# novel miRNAs:
for file in \
`ls $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/$outfile; \
cd $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/$outfile; \
miRDeep2.pl $file \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/${outfile}.arf \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/bta_mature-miRNA.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/other_mature-miRNA.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/bta_hairpin-miRNA.fa -t Cow" \
>> miRdeep2.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 1 miRdeep2.sh miRdeep2.sh.
for script in `ls miRdeep2.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_confidence
cd !$

# Create bash script for identification and quantification of
# high confidence miRNAs:
for file in \
`ls $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo \
"mkdir $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_confidence/$outfile; \
cd $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_confidence/$outfile; \
miRDeep2.pl $file \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mapper/${outfile}.arf \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_conf_mature.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/other_mature-miRNA.fa \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_conf_hairpin_clean.fa \
-t Cow" \
>> miRdeep2_high-conf.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 1 miRdeep2_high-conf.sh miRdeep2_high-conf.sh.
for script in `ls miRdeep2_high-conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all read count files from regular and high confidence miRNAs for
# transfering into laptop using WinSCP:
mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/mirdeep.pl/regular
cd !$
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/SMART* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*mirdeep\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/mirdeep.pl/regular/${outfile}_exp_mirdeep.csv; \
done

mkdir -p $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/mirdeep.pl/high_conf
cd !$
for file in `find $HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/mirdeep2/Project3/mirdeep/high_confidence/SMART* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*mirdeep\/.*\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/BTB_SFI_Project/WP2/miRNA_LibraryKit_Study/Counts/mirdeep2/Project3/mirdeep.pl/high_conf/${outfile}_exp_mirdeep_hc.csv; \
done

########################################
# R analysis of gene counts with edgeR #
########################################

# Subsequent sense genes analyses were performed using the R statistical
# and the edgeR package. Please refer to file:
# BioValidation-miRNAseq_edgeR_pipeline.R

