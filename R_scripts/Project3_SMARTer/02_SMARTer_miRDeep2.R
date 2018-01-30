##########################################################################
#                         Project3: SMARTer_kit						        		 #
#   miRNA-seq library kit comparison using PBMC-isolated miRNA from			 # 
#               Field-infected (n=4) and non-infected cattle (n=4)  	   #
#                                                                        #
#              --- R workflow for the miRDeep2 approach ---              #
#                                 Part 2                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: 

# Authors of current version (2.0.0): Sarah Faherty O'Donnell
# DOI badge of current version:
# Last updated on 27/01/2018

############################################
# 14 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(statmod)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(biobroom)
library(Cairo)
library(extrafont)
library(forcats)
library(ggridges)
library(PerformanceAnalytics)
library(ggrepel)
library(cowplot)

# Uncomment functions below to install packages in case you don't have them

#install.packages("ggridges")
#install.packages("PerformanceAnalytics")
#install.packages("ggrepel")
#install.packages("cowplot")
#install.packages("statmod")
##################################
# 15 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_SMARTer_miRDeep2.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

##########################################################
# 16 Tidy DGElist for exploratory data analysis plotting #
##########################################################

tidy_dgelist <- tidy(dgelist_filt, addSamples = TRUE)

tidy_dgelist
View(tidy_dgelist)

# Convert samples to factors
tidy_dgelist$sample %<>%
  factor(levels = c("SMARTer4_S17", "SMARTer4_S18", "SMARTer4_S19", "SMARTer4_S20",
                    "SMARTer4_S21", "SMARTer4_S22", "SMARTer4_S23", "SMARTer4_S24"))
# CHeck samples factors
levels(tidy_dgelist$sample)

# Check group factors
levels(tidy_dgelist$group)

# Check data frame
tidy_dgelist

########################################################
# 17 Plot: density of filtered gene counts per library #
########################################################

ggplot(tidy_dgelist, aes(x = log10(count + 1),
                         y = sample)) +
  scale_y_discrete(limits = rev(levels(tidy_dgelist$sample))) +
  geom_density_ridges(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("condition",
                    values = c("#b2b2b2", rep("#e06377", 6))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("Density of filtered gene counts per sample (", method, ")")) +
  ylab("sample") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_filt


density_filt

# Export high quality image
ggsave(paste(method, "_density-filt.pdf", sep = ""),
       plot      = density_filt,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

##########################################################
# 18 Plot: Gene expression correlation within condition groups #
##########################################################

# Check gene expression correlation within infected cohort (libraries 17-20) and 
# within control cohort (libraries 21-24) to 
# identify potential outliers (CPM values not required with Spearman)


png(file = file.path(paste0(imgDir, method, "_Infected", "_cor.png", sep = "")),
    width = 1500, height = 1500, units = "px")
chart.Correlation(R = log(x = (dgelist_norm$counts[, grep(
  pattern = "S(17|18|19|20)$", x = colnames(dgelist_norm$counts),
  perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
  main = " Infected group expression correlation")
dev.off()

png(file = file.path(paste0(imgDir, method, "_Control", "_cor.png", sep = "")),
    width = 1500, height = 1500, units = "px")
chart.Correlation(R = log(x = (dgelist_norm$counts[, grep(
  pattern = "S(21|22|23|24)$", x = colnames(dgelist_norm$counts),
  perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
  main = " Control group expression correlation")
dev.off()

##################################################
# 19 Define function for getting MDS coordinates #
##################################################

getBCVcoord <- function(dgelst, time_pattrn) {
  
  mds <- plotMDS.DGEList(x = dgelst[ , grepl(paste(time_pattrn,
                                                   collapse = "|"),
                                             x = colnames(dgelst))],
                         plot = FALSE,
                         method = "bcv")
  
  mds_coord <- mds$cmdscale.out # Get coords to plot with ggplot2
  
  mds_coord %<>% # Tidy coords
    tidy() %>%
    dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
    dplyr::mutate(group = sample)
  
  mds_coord$group %<>% # Clean group info for plotting
    str_replace("SMARTer4_", "") %>%
    str_replace("S(21|22|23|24)", "Control") %>%
    str_replace("S(17|18|19|20)", "Infected")
  
  return(mds_coord)
  
}

##########################
# 20 Get MDS coordinates #
##########################

# Get tidy MDS coordinates for all samples 
all_coord <- getBCVcoord(dgelist_norm, c("_S17", "_S18", "_S19",
                                         "_S20", "_S21", "_S22", "_S23", "_S24"))


#####################################
# 21 Plot: MDS plot #
#####################################

MDS_all <- ggplot(all_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Condition",
                      values = c("#b2b2b2" , "#e06377")) +
  scale_shape_manual("Condition",
                     values = c(19, 17)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("SMARTer kit") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

MDS_all

# Export high quality image
ggsave(filename  = paste(method, "_MDS_all.pdf", sep = ""),
       plot      = MDS_all,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 10,
       units     = "in")

#################################################
# 22 Create a design matrix                    #
#################################################

# Create a design matrix with condition as blocking factor
block_condition <- model.matrix(~ group,
                                data = dgelist_norm$samples)

dim(block_condition)
dim(dgelist_norm$samples)
block_condition

# Rename design matrix columns for simplicity
colnames(block_condition) %<>%
  str_replace("group", "")


head(block_condition)

# Output the design matrix info
as.data.frame(block_condition) %>%
  rownames_to_column(var = "Samples") %>%
  write_csv(path = file.path(paste0(tablesDir, method, "_design-matrix.csv")),
            col_names = TRUE) 



#########################################
# 23 Estimate the dispersion parameters #
#########################################

# Common and trended dispersions are estimated with the
# Cox-Reid method and tagwise dispersions with the
# empirical Bayes method
dgelist_disp <- estimateDisp.DGEList(y       = dgelist_norm,
                                     design  = block_condition,
                                     robust  = TRUE,
                                     verbose = TRUE)

names(dgelist_disp)
dim(dgelist_disp)
# Check the calculated dispersion
dgelist_disp$common.dispersion

# Check the calculated dispersion's square root,
# which corresponds to the biological coefficient of variation (BCV)
# Compare the different kits BCVs
sqrt(dgelist_disp$common.dispersion)
# BCV of each gene (useful for possible reference genes)
sqrt(dgelist_disp$tagwise.dispersion)
min(sqrt(dgelist_disp$tagwise.dispersion))
max(sqrt(dgelist_disp$tagwise.dispersion))

# Create a (dataframe) matrix of the tagwise dispersion associated with gene information
Tagwisedisp <- as.data.frame(cbind(dgelist_disp$genes,
                                   dgelist_disp$tagwise.dispersion))
head(Tagwisedisp)
dim(Tagwisedisp)

# Output tagwise dispersion values with gene info
Tagwisedisp %>%
  rownames_to_column(var = "miRBase_ID") %>%
  write_csv(path = file.path(paste0(tablesDir, method, "_Tagwise_dispersion.csv")),
            col_names = TRUE)

################################
# 24 Plot: BCV and dispersions #
################################

# Create a dataframe with the dispersion values
names(dgelist_disp)

Disp <- as.data.frame(cbind(dgelist_disp$genes,
                            dgelist_disp$tagwise.dispersion,
                            dgelist_disp$common.dispersion,
                            dgelist_disp$trended.dispersion,
                            dgelist_disp$AveLogCPM))

head(Disp)
dim(Disp)

colnames(Disp) %<>%
  str_replace("dgelist_disp\\$", "")

Disp %<>%
  dplyr::mutate(type_point = "Tagwise dispersion") %>%
  dplyr::mutate(type_hline = "Common dispersion") %>%
  dplyr::mutate(type_smooth = "Trended dispersion")

# Plot all dispersions
ggplot(Disp) +
  geom_point(aes(x = AveLogCPM,
                 y = sqrt(tagwise.dispersion),
                 fill = type_point),
             alpha = 0.5) +
  geom_hline(aes(yintercept = sqrt(common.dispersion),
                 colour = type_hline)) +
  geom_smooth(aes(x = AveLogCPM,
                  y = sqrt(trended.dispersion),
                  colour = type_smooth),
              linetype = 2) +
  scale_fill_manual("", values = c("black")) +
  scale_colour_manual("", values = c("red", "blue")) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("Estimated dispersions for SMARTer Kit") +
  xlab(expression(paste("Average ", log[2],"CPM"))) +
  ylab("Biological Coefficient of Variation") -> dgelist_BCV

dgelist_BCV

# Output high resolution plot
ggsave(paste(method, "_BCV.pdf", sep = ""),
       plot      = dgelist_BCV,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 10,
       units     = "in")

#######################
# 25 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 26 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 3 of this analysis #
######################################

# File: 03_SMARTer_miRDeep2.R