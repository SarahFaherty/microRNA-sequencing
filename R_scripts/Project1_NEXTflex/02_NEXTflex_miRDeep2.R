##########################################################################
#                         Project1: NEXTflex_kit						        		 #
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
# Last updated on 23/01/2018

############################################
# 14 Load and/or install required packages #
############################################

library(here)
library(edgeR)
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

##################################
# 15 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_NEXTflex_miRDeep2.RData")

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
  factor(levels = c("NEXTflex_S1", "NEXTflex_S2", "NEXTflex_S3", "NEXTflex_S4",
                    "NEXTflex_S5", "NEXTflex_S6", "NEXTflex_S7", "NEXTflex_S8"))
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

# Check gene expression correlation within infected cohort (libraries 1-4) and 
# within control cohort (libraries 5-8) to 
# identify potential outliers (CPM values not required with Spearman)


png(file = file.path(paste0(imgDir, method, "_Infected", "_cor.png", sep = "")),
      width = 1500, height = 1500, units = "px")
chart.Correlation(R = log(x = (dgelist_norm$counts[, grep(
    pattern = "S(1|2|3|4)$", x = colnames(dgelist_norm$counts),
    perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
    main = " Infected group expression correlation")
dev.off()

png(file = file.path(paste0(imgDir, method, "_Control", "_cor.png", sep = "")),
    width = 1500, height = 1500, units = "px")
chart.Correlation(R = log(x = (dgelist_norm$counts[, grep(
  pattern = "S(5|6|7|8)$", x = colnames(dgelist_norm$counts),
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
    str_replace("NEXTflex_", "") %>%
    str_replace("S(5|6|7|8)", "Control") %>%
    str_replace("S(1|2|3|4)", "Infected")
  
  return(mds_coord)
  
}

##########################
# 20 Get MDS coordinates #
##########################

# Get tidy MDS coordinates for all samples 
all_coord <- getBCVcoord(dgelist_norm, c("_S1", "_S2", "_S3",
                                         "_S4", "_S5", "_S6", "_S7", "_S8"))


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
  ggtitle("NEXTflex kit") +
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


#######################
# 26 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 27 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 3 of this analysis #
######################################

# File: 03_NEXTflex_miRDeep2.R