##########################################################################
#                         Project1: NEXTflex_kit						        		 #
#   miRNA-seq library kit comparison using PBMC-isolated miRNA from			 # 
#               Field-infected (n=4) and non-infected cattle (n=4)  	   #
#                                                                        #
#              --- R workflow for the miRDeep2 approach ---              #
#                                 Part 3                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge:

# Authors of current version (2.0.0): Sarah Faherty O'Donnell
# DOI badge of current version:
# Last updated on 25/01/2018

############################################
# 27 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(Cairo)
library(extrafont)
library(VennDiagram)
library(forcats)
library(treemap)

# Uncomment functions below to install packages in case you don't have them

#install.packages("VennDiagram")
#install.packages("treemap")

##################################
# 28 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_NEXTflex_miRDeep2.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

####################
# 29 Fit GLM model #
####################

# Fit a negative binomial generalized log-linear model
# using the design matrix and calculated dispersions
dgelist_fit <- glmFit(y = dgelist_disp,
                      design = block_condition)

names(dgelist_fit)
colnames(dgelist_fit$design)

#####################################################
# 30 Determine differential expression by fitting a #
# negative binomial GLM with Likelihood Ratio Tests #
#####################################################

# Test for differential expression between control (intercept) and infected,
# using the coefficients from degelist_fit$design

Infected <- glmLRT(dgelist_fit, coef = "Infected")
testDE_Infected <- topTags(object        = Infected,
                     n             = "inf",
                     adjust.method = "BH",
                     sort.by       = "PValue")

dim(testDE_Infected$table)
head(testDE_Infected$table)




#####################################
# 31 Output all genes tested for DE #
#####################################

# Infected
testDE_Infected$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_AllGenes.csv")),
            col_names = TRUE)

##############################################
# 32 Filter genes considered DE (FDR < 0.05) #
##############################################

# Infected
testDE_Infected$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> Infected_FDR_05

Infected_FDR_05

### Output lists of genes considered DE (FDR < 0.05)
DElists <- list(Infected_FDR_05)
DEfile <- c(paste0(c("Infected_FDR_05"), "_genes.csv"))
DEpaths <- file.path(paste0(tablesDir, method, "_", DEfile))

pwalk(list(DElists, DEpaths),
      write_csv,
      col_names = TRUE)

###############################################
# 33 Plot: Treemaps of DE genes (FDR < 0.05) #
###############################################

# Get numbers of up and down regulated genes

list_DE <- list(Infected_FDR_05)
names(list_DE) <- c("Infected")

Up_Down <- map_df(list_DE,
                  ~ dplyr::count(.x,
                                 up = sum(logFC > 0),
                                 down = sum(logFC < 0),
                                 zero = sum(logFC == 0)),
                  .id = "Condition")

# Check data frame
Up_Down

# Time point as factor
Up_Down$Condition %<>%
  factor() %>%
  fct_inorder()

# Plotting labels
Up_Down %<>% dplyr::mutate(labelsUp = paste(up, "Increased", sep = ' '))
Up_Down %<>% dplyr::mutate(labelsDown = paste(down,"Decreased", sep = ' '))

# Check data frame
Up_Down

# Plot chart increased expression
# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, method, "_tree_up.pdf")),
          width    = 8,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
treemap(Up_Down,
        index             = "labelsUp", "labelsDown",
        vSize             = "up", "down",
        type              = "index",
        palette           = "PRGn",
        title             = "NEXTflex kit: Differentially expressed miRNA",
        fontsize.title    = 14,
        fontfamily.title  = "Calibri",
        fontfamily.labels = "Calibri",
        fontsize.labels   = 16)
dev.off()


#######################
# 34 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 35 Save R session info #
##########################

devtools::session_info()

#######
# END #
#######