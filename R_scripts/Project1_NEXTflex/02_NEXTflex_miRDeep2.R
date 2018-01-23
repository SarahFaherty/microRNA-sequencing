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

# Clean animal IDs
tidy_dgelist$animal %<>%
  stringr::str_replace("A", "") %>%
  fct_inorder()

# Change time point info for labels
tidy_dgelist$time.point %<>%
  str_replace("pre1", "-1 wk") %>%
  str_replace("pre2", "-2 wk") %>%
  str_replace("W1$", "+1 wk") %>%
  str_replace("W2", "+2 wk") %>%
  str_replace("W6", "+6 wk") %>%
  str_replace("W10", "+10 wk") %>%
  str_replace("W12", "+12 wk") %>%
  factor(levels = c("-2 wk", "-1 wk", "+1 wk", "+2 wk",
                    "+6 wk", "+10 wk", "+12 wk"))

# Combine animal and time point info for
# plotting labels
tidy_dgelist %<>%
  dplyr::mutate(labels = paste0(time.point, "_", animal))

tidy_dgelist$labels %<>%
  factor(levels = c("-2 wk_6511", "-2 wk_6514", "-2 wk_6520", "-2 wk_6522", "-2 wk_6526",
                    "-2 wk_6635", "-2 wk_6636", "-2 wk_6637", "-2 wk_6644", "-2 wk_6698",
                    "-1 wk_6511", "-1 wk_6514", "-1 wk_6520", "-1 wk_6522", "-1 wk_6526",
                    "-1 wk_6635", "-1 wk_6636", "-1 wk_6637", "-1 wk_6644", "-1 wk_6698",
                    "+1 wk_6511", "+1 wk_6514", "+1 wk_6520", "+1 wk_6522", "+1 wk_6526",
                    "+1 wk_6635", "+1 wk_6636", "+1 wk_6637", "+1 wk_6644", "+1 wk_6698",
                    "+2 wk_6511", "+2 wk_6514", "+2 wk_6520", "+2 wk_6522", "+2 wk_6526",
                    "+2 wk_6635", "+2 wk_6636", "+2 wk_6637", "+2 wk_6644", "+2 wk_6698",
                    "+6 wk_6511", "+6 wk_6514", "+6 wk_6520", "+6 wk_6522", "+6 wk_6526",
                    "+6 wk_6635", "+6 wk_6636", "+6 wk_6637", "+6 wk_6644", "+6 wk_6698",
                    "+10 wk_6511", "+10 wk_6514", "+10 wk_6520", "+10 wk_6522", "+10 wk_6526",
                    "+10 wk_6635", "+10 wk_6636", "+10 wk_6637", "+10 wk_6644", "+10 wk_6698",
                    "+12 wk_6511", "+12 wk_6514", "+12 wk_6520", "+12 wk_6522", "+12 wk_6526",
                    "+12 wk_6635", "+12 wk_6636", "+12 wk_6637", "+12 wk_6644", "+12 wk_6698"))

# Check factors
levels(tidy_dgelist$labels)

# Check data frame
tidy_dgelist

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