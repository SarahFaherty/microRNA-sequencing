##########################################################################
#                         Project1: NEXTflex_kit						        		 #
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
# Last updated on 22/01/2018

############################################
# Installing and loading required packages
############################################
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
install.packages("here")
install.packages("tidyverse")
install.packages("Cairo")
install.packages("extrafont")

# This paskage wasn't installed due to version of R
# install.packages("biobroom")
# package 'biobroom' is not available (for R version 3.4.3)

library(here)
library(limma)
library(edgeR)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
#library(biobroom)
library(Cairo)
library(extrafont)

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()
# Output from here():
# [1] "C:/Users/Sarah/Documents"

# Define variables for subdirectories
countsDir <- here("quant_mature_counts/")
imgDir <- here("Figures/")
tablesDir <- here("Tables/")

# Define the method used for miRNA identification
method <- "miRDeep2"

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()
