##########################################################################
#                         Project3: SMARTer_kit						        		 #
#   miRNA-seq library kit comparison using PBMC-isolated miRNA from			 # 
#               Field-infected (n=4) and non-infected cattle (n=4)  	   #
#                                                                        #
#              --- R workflow for the miRDeep2 approach ---              #
#                                 Part 1                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: 

# Authors of current version (2.0.0): Sarah Faherty O'Donnell
# DOI badge of current version:
# Last updated on 27/01/2018

############################################
# Installing and loading required packages
############################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("biobroom")

#install.packages("here")
#install.packages("tidyverse")
#install.packages("Cairo")
#install.packages("extrafont")
#install.packages("devtools")

library(here)
library(limma)
library(edgeR)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(biobroom)
library(Cairo)
library(extrafont)

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()
# Output from here():
# [1] "C:/Users/Sarah/Dropbox/@WP2/miRNA_seqLibraykit/miSeq_Bioinformatics/R_scripts/Project3_SMARTer"


# Define variables for subdirectories
# single / doesn't work in windows
# used // to have a sigle / after "quant_mature_counts"
countsDir <- here("quant_mature_counts//")
imgDir <- here("Figures//")
tablesDir <- here("Tables//")

# Define the method (and project) used for miRNA identification
method <- "SMARTer_miRDeep2"

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

########################################
# 03 Import miRDeep2 quantifier counts #
########################################

# Create a vector of file names
files <- list.files(path        = countsDir,
                    pattern     = "^SMART",
                    all.files   = TRUE,
                    full.names  = FALSE,
                    recursive   = FALSE,
                    ignore.case = FALSE)

files
names(files) <- files
length(files)

# Create a dataframe with raw counts for all samples
files %>%
  map_df(~ read_tsv(paste0(countsDir, .x), col_names = TRUE),
         .id = "filename") %>%
  dplyr::rename(gene_name = `#miRNA`, precursor_name = precursor) %>%
  dplyr::select(-c(total, seq, `seq(norm)`)) %>%
  spread(filename, read_count) -> rawCounts

# Check data frame
rawCounts
# To view a tabulated (excel -like) form of the rawCounts:
View(rawCounts)

# Clean column names and add a prefix to avoid sample names starting
# with numbers
colnames(rawCounts) %<>%
  str_replace("_expressed.csv", "") %>%
  str_replace("-\\d", "")

# Check data frame
rawCounts
View(rawCounts)

###############################################################
# 04 Gene annotation using information obtained from GTF file #
###############################################################

# Read in the annotation information
miRNA_info <- read.table(file   = "miRNA_Btaurus.txt",
                         header = TRUE,
                         sep    = "\t",
                         quote  = "")

head(miRNA_info)
dim(miRNA_info)

# Determine which miRNAs have identical mature sequence but
# originate from different precursors and add it to annotation
miRNA_info %<>%
  dplyr::group_by(sequence) %>%
  dplyr::summarise(identical_sequence = paste(gene_id, collapse = ",")) %>%
  dplyr::inner_join(miRNA_info, by = "sequence") %>%
  dplyr::select(starts_with("gene"), chromosome, contains("position"),
                strand, sequence, starts_with("precursor"),
                identical_sequence)

miRNA_info

# Merge gene annotation with the counts and output data
annotCounts <- dplyr::inner_join(x = miRNA_info,
                                 y = rawCounts,
                                 by = c("gene_name", "precursor_name"))

View(annotCounts)

annotCounts %>%
  write_csv(file.path(paste0(tablesDir, method, "_Raw-counts.csv",
                             sep = "")),
            col_names = TRUE)

# Add genes as rownames
annotCounts %<>%
  data.frame() %>%
  column_to_rownames(var = "gene_id")

head(annotCounts)
dim(annotCounts)

##############################################
# 05 Create groups and variables for DGElist #
##############################################

# Select counts only
annotCounts %>%
  dplyr::select(starts_with("SMART")) -> counts

head(counts)
dim(counts)

# Define control and treatment groups
colnames(counts) %>%
  as.tibble() %>%
  dplyr::select(sample = value) %>%
  dplyr::mutate(group = sample) -> group

group$group %<>%
  str_replace("SMARTer4_", "") %>%
  str_replace("S(21|22|23|24)", "Control") %>%
  str_replace("S(17|18|19|20)", "Infected")


group %<>%
  data.frame() %>%
  column_to_rownames(var = "sample")

head(group)
dim(group)
# Checking that group rownames match the counts columns
identical(rownames(group), colnames(counts))

# Select annotation only
annotCounts %>%
  dplyr::select(-starts_with("SMART")) -> gene_annot

head(gene_annot)
dim(gene_annot)

#####################
# 06 Create DGElist #
#####################

# Create DGElist
dgelist <- DGEList(counts       = counts,
                   group        = group$group,
                   genes        = gene_annot,
                   lib.size     = NULL,
                   norm.factors = NULL,
                   remove.zeros = FALSE)

names(dgelist)
dim(dgelist)
head(dgelist$counts)
head(dgelist$samples)
head(dgelist$genes)

################################################
# 07 Additional sample information for DGElist #
################################################

# Convert treatment group to factor
dgelist$samples$group %<>%
  factor(levels = c("Control", "Infected"))

# Check order of factor levels
levels(dgelist$samples$group)

# Check data frame
head(dgelist$samples)

# Output sample information
dgelist$samples %>%
  rownames_to_column(var="Sample_IDs") %>%
  write_csv(file.path(paste0(tablesDir, method, "_samples-info.csv", sep = "")),
            col_names = TRUE)
################################################
# 08 Density plot: raw gene counts per library #
################################################

# Tidy DGElist and plot data
dgelist %>%
  tidy() %>%
  ggplot() +
  geom_density(aes(x     = log10(count + 1),
                   group = sample), size = 0.1) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(method) +
  ylab("Density of raw gene counts per sample") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_raw


density_raw

# Export image
ggsave(paste0(method, "_Raw-density.pdf", sep = ""),
       plot      = density_raw,
       device    = cairo_pdf,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 10,
       path      = imgDir)

###########################################
# 09 Remove zero and lowly expressed tags #
###########################################

# Filter non-expressed tags (all genes that have zero counts in all samples)
dgelist_no_zeros <- dgelist[rowSums(dgelist$counts) > 0, ]
dim(dgelist_no_zeros$counts)
head(dgelist_no_zeros$counts)
colnames(dgelist_no_zeros$counts)

# Filter lowly expressed tags, retaining only tags with
# more than 50 counts per million in 4 or more libraries
# (4 libraries correspond to 4 animals in either group
# (infected or control groups)
dgelist_filt <- dgelist_no_zeros[rowSums(cpm(dgelist_no_zeros) > 50) >= 4, ]
dim(dgelist_filt$counts)
head(dgelist_filt$counts)

# Ouptut filtered counts
dgelist_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "miRBaseID") %>%
  write_csv(file.path(paste0(tablesDir, method, "_Filt_counts.csv",
                             sep = "")),
            col_names = TRUE)

##############################
# 10 Recompute library sizes #
##############################

dgelist_filt$samples$lib.size <- colSums(dgelist_filt$counts)
head(dgelist_filt$samples)
head(dgelist$samples)

###########################################################################
# 11 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after
# calculating normalisation factors
dgelist_norm <- calcNormFactors(dgelist_filt, method = "TMM")
head(dgelist_norm$samples)

#######################
# 12 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 13 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 02_SMARTer_miRDeep2.R