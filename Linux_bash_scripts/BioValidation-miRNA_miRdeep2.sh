##############################################################################
# miRNA-seq Biological Validation of Blood Serum from Field-infected Animals #
#    --- Linux bioinformatics workflow for known and novel miRNAs  ---       #
#                     -- Method 2: miRDeep2 software --                      #
##############################################################################
# Authors: Nalpas, N.C.; Correia, C.N. (2014) 
# Zenodo DOI badge: http://dx.doi.org/10.5281/zenodo.16164
# Version 1.2.0
# Last updated on: 14/07/2016

########################################################################
# Merge and uncompress miRNA-seq FASTQ files to be used with miRDeep2  #
########################################################################

# Create and enter temporary workiing directory:
mkdir $HOME/scratch/miRNAseqValidation/fastq_sequence/tmp
cd !$

# Create bash script to uncompress and merge trim.fast.gz from lanes
# 005 and 006 for each library:
for file in `find $HOME/scratch/miRNAseqValidation/fastq_sequence \
-name *_L005_R1_001_trim.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_L005_)/_L006_/'`; \
sample=`basename $file | perl -p -e 's/(E\d+_)\w+_L\d+_R\d_\d*_trim.fastq.gz/$1/'`; \
echo "zcat $file $file2 > ${sample}trim.fastq" \
>> uncompress_merge.sh; \
done

# Run script:
chmod 755 uncompress_merge.sh
nohup ./uncompress_merge.sh &

#######################################
# Index reference genome using Bowtie #
#######################################

# Required software is Bowtie 1.1.0, consult manual for details:
# http://bowtie-bio.sourceforge.net/manual.shtml

# Create and enter the Index reference genome directory:
mkdir -p /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0
cd !$

# Index the reference genome UMD3.1 using Bowtie:
nohup bowtie-build \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/source_file/Btau_UMD3.1_multi.fa \
Btau_UMD3.1_multi &

#############################################################
# Preprocessing of miRNA-seq data using miRDeep2: mapper.pl #
#############################################################

# Required software is miRDeep2 v.2.0.0.8, consult manual/tutorial for details:
# https://www.mdc-berlin.de/8551903/en/

# Create and enter the miRDeep2 directory for mapping work:
mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/mapper
cd !$

# Create symbolic links to FASTQ files:
for file in \
`ls $HOME/scratch/miRNAseqValidation/fastq_sequence/tmp/*_trim.fastq`; \
do ln -s \
$file $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/`basename $file`; \
done

# Run mapper.pl in one FASTQ file to see if it's working well:
mapper.pl E10_trim.fastq -e -h -m -o 3 -l 17 -r 50 -q -v -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0/Btau_UMD3.1_multi \
-s test_collapsed.fa -t test.arf

# Create bash script to map miRNA reads to the reference genome:
for file in `ls *_trim.fastq`; \
do outfile=`basename $file | perl -p -e 's/_trim\.fastq//'`; \
echo "mapper.pl $file -e -h -m -o 3 -l 17 -r 50 -q -v -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0/Btau_UMD3.1_multi \
-s ${outfile}_collapsed.fa -t ${outfile}.arf" \
>> mapper.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 12 mapper.sh mapper.sh.
for script in `ls mapper.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

################################################################
# Quantification of known miRNAs using miRDeep2: quantifier.pl #
################################################################

# Create and enter the working directory:
mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier
cd !$

# Run quantifier.pl in one FASTA file to see if it's working well:
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/E10_collapsed.fa -t bta

# Create a shell script to quantify the mapped mature miRNAs:
# [You will need the mature and precursor (hairpin) miRNA FASTA files for
# Bos taurus sequences only. Please refer to the 
# BioValidation-miRNA-seq_QC_filter.sh script]
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/$outfile; \
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $file -t bta" >> quantifier.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 quantifier.sh quantifier.sh.
for script in `ls quantifier.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence
cd !$
# Create a shell script to quantify the mapped high confidence mature and
# precursor (hairpin) miRNAs:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/$outfile; \
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_hairpin.fa \
-r $file -t bta" >> quantifier_high_conf.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 quantifier_high_conf.sh quantifier_high_conf.sh.
for script in `ls quantifier_high_conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all read counts files from regular and high confidence miRNAs for
# transfering into laptop using WinSCP:
mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep2/regular
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep2/regular/${outfile}_expressed.csv; \
done

mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep2/high_conf
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/.*\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep2/high_conf/${outfile}_expressed_high_conf.csv; \
done

########################################################################
# Identification of known and novel miRNAs using miRDeep2: miRDeep2.pl #
######################################################################## 

# Create and enter working directory:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep
cd !$

# Copy and modify the required fasta files since miRdeep2 software
# requires no space in headers and no characters other than acgtunACGTUN in
# the sequences:
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/source_file/Btau_UMD3.1_multi.fa \
./Btau_UMD3.1_multi.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
./bta_mature-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
./high_conf_mature.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/other_mature-miRNA.fa \
./other_mature-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
./bta_hairpin-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_hairpin.fa \
./high_conf_hairpin.fa

perl -p -i -e 's/^(>.*?)\s.*$/$1/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*$)/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_mature.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_hairpin-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) (.*?) (.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_hairpin.fa
sed -e '/^[^>]/s/[^ATGCatgc]/N/g' high_conf_hairpin.fa > high_conf_hairpin_clean.fa

# Run mirdeep.pl in one FASTA file to see if it's working well:
miRDeep2.pl $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/E10_collapsed.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mapper/E10.arf \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_hairpin-miRNA.fa -t Cow

# Create bash script for identification and quantification of known and
# novel miRNAs:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/$outfile; \
miRDeep2.pl $file \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mapper/${outfile}.arf \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_hairpin-miRNA.fa -t Cow" \
>> miRdeep2.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 miRdeep2.sh miRdeep2.sh.
for script in `ls miRdeep2.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence
cd !$

# Create bash script for identification and quantification of
# high confidence miRNAs:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo \
"mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence/$outfile; \
miRDeep2.pl $file \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mapper/${outfile}.arf \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_mature.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_hairpin_clean.fa \
-t Cow" \
>> miRdeep2_high-conf.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 miRdeep2_high-conf.sh miRdeep2_high-conf.sh.
for script in `ls miRdeep2_high-conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all read counts files from regular and high confidence miRNAs for
# transfering into laptop using WinSCP:
mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep2/mirdeep.pl/regular
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*mirdeep\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep2/mirdeep.pl/regular/${outfile}_exp_mirdeep.csv; \
done

mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep2/high_conf
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*mirdeep\/.*\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep2/mirdeep.pl/high_conf/${outfile}_exp_mirdeep_hc.csv; \
done

########################################
# R analysis of gene counts with edgeR #
########################################

# Subsequent sense genes analyses were performed using the R statistical
# and the edgeR package. Please refer to file:
# BioValidation-miRNAseq_edgeR_pipeline.R

