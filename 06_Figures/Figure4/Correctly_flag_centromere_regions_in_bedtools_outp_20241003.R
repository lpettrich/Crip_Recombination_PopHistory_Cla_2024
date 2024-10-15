###########################################################
# #
# Statistics recombination spots and Cla-element location #
#   Script by Laura Pettrich  #
#   REVISED   #
# Juli 2024   #
###########################################################
# this scripts flags all inserts and recombinaton events on arm1, arm2 and in the centromere
# analysis will be limited to insertions only present in same area

## Clean environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
#library(chromoMap)
#library(scales)
library(ggplot2)

# Set directory
# setwd("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")
setwd("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/bedtools-closest/run12/")


getwd()

#######################################################################################
# READ IN  DATA
#######################################################################################
#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
# Here you define vectors for every data table you have
# each population
pop <- c("MF_Chr1", "MF_Chr2", "MF_Chr3", "MF_Chr4", 
         "MG_Chr1", "MG_Chr2", "MG_Chr3", "MG_Chr4", 
         "NMF_Chr1", "NMF_Chr2", "NMF_Chr3", "NMF_Chr4", 
         "SI_Chr1", "SI_Chr2", "SI_Chr3", "SI_Chr4", 
         "SS_Chr1", "SS_Chr2", "SS_Chr3", "SS_Chr4")

# each chromosome
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

# each individual
ind_numbers <- 1:4
pop_prefix <- c("MF", "MG", "NMF", "SI", "SS")
ind_prefix <- paste(rep(pop_prefix, each = 4), 1:4, sep = "")
# Find indices where "MG" appears and change numbers
mg_indices <- grep("^MG", ind_prefix)
ind_prefix[mg_indices] <- paste("MG", 2:5, sep = "")
print(ind_prefix)

chr_suffixes <- paste("Chr", 1:4, sep = "")

ind <- paste(rep(ind_prefix, each = 4), chr_suffixes, sep = "_")
print(ind)



# Read in centromere ranges
cenrange <- read.csv("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/CentromereRange.csv", sep = ",", header = TRUE)
head(cenrange)

#-----------------------------------
# SPLIT CHROMOSOME ARMS
#-----------------------------------
#---------------------------------------------------------------------------------------------------------

# MF

# Read in bedtools closest results (10kb)
dfm <- read.csv("bedtools_closest_ouput_100kb.csv", sep = ",", header = TRUE)
head(dfm)

# Select intersting columns
dfm <- dfm %>%
  select("chromosome", "start", "end", "rho", 
         "chromosome_ins", "start_ins", "end_ins", "no_insert", "dist_to_Cla", 
         "ind_chr", "pop")
# Create the chr_arms column, marking Chr4 with NA and applying the filter to other chromosomes
dfm <- dfm %>%
  filter(pop == "MF") %>%
  rowwise() %>%
  mutate(
    chr_arms = case_when(
      chromosome == "Chr4" ~ NA_character_,  # Mark everything from Chr4 with NA
      
      # Tag as arm1 if start, end, start_ins, end_ins are all less than censtart
      (start < cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end < cenrange$censtart[chromosome == cenrange$Chromosome] &
         start_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome]) ~ "arm1",
      
      # Tag as arm2 if start, end, start_ins, end_ins are all greater than cenend
      (start > cenrange$cenend[chromosome == cenrange$Chromosome] & 
         end > cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome] &
         end_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "arm2",
      
      # Tag centromer if start is greater or equal and cenend is smaller or equal
      (start >= cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end <= cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins >= cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins <= cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "centromere",
      
      TRUE ~ NA_character_  # This will leave NA if neither condition is met
    )
  ) %>%
  ungroup()

# View the updated dataframe
head(dfm)

write.table(dfm, "bedtools_closest_ouput_100kb_arms_MF_20241003.csv", sep = ",", col.names = TRUE)
dfm <- read.csv("bedtools_closest_ouput_100kb_arms_MF_20241003.csv", sep = ",", header = TRUE)

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

# MG

# Read in bedtools closest results (10kb)
dfm <- read.csv("bedtools_closest_ouput_100kb.csv", sep = ",", header = TRUE)
head(dfm)

# Select intersting columns
dfm <- dfm %>%
  select("chromosome", "start", "end", "rho", 
         "chromosome_ins", "start_ins", "end_ins", "no_insert", "dist_to_Cla", 
         "ind_chr", "pop")
# Create the chr_arms column, marking Chr4 with NA and applying the filter to other chromosomes
dfm <- dfm %>%
  filter(pop == "MG") %>%
  rowwise() %>%
  mutate(
    chr_arms = case_when(
      chromosome == "Chr4" ~ NA_character_,  # Mark everything from Chr4 with NA
      
      # Tag as arm1 if start, end, start_ins, end_ins are all less than censtart
      (start < cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end < cenrange$censtart[chromosome == cenrange$Chromosome] &
         start_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome]) ~ "arm1",
      
      # Tag as arm2 if start, end, start_ins, end_ins are all greater than cenend
      (start > cenrange$cenend[chromosome == cenrange$Chromosome] & 
         end > cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome] &
         end_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "arm2",
      
      # Tag centromer if start is greater or equal and cenend is smaller or equal
      (start >= cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end <= cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins >= cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins <= cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "centromere",
      
      TRUE ~ NA_character_  # This will leave NA if neither condition is met
    )
  ) %>%
  ungroup()

# View the updated dataframe
head(dfm)

write.table(dfm, "bedtools_closest_ouput_100kb_arms_MG_20241003.csv", sep = ",", col.names = TRUE)
dfm <- read.csv("bedtools_closest_ouput_100kb_arms_MG_20241003.csv", sep = ",", header = TRUE)

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

# NMF

# Read in bedtools closest results (10kb)
dfm <- read.csv("bedtools_closest_ouput_100kb.csv", sep = ",", header = TRUE)
head(dfm)

# Select intersting columns
dfm <- dfm %>%
  select("chromosome", "start", "end", "rho", 
         "chromosome_ins", "start_ins", "end_ins", "no_insert", "dist_to_Cla", 
         "ind_chr", "pop")
# Create the chr_arms column, marking Chr4 with NA and applying the filter to other chromosomes
dfm <- dfm %>%
  filter(pop == "NMF") %>%
  rowwise() %>%
  mutate(
    chr_arms = case_when(
      chromosome == "Chr4" ~ NA_character_,  # Mark everything from Chr4 with NA
      
      # Tag as arm1 if start, end, start_ins, end_ins are all less than censtart
      (start < cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end < cenrange$censtart[chromosome == cenrange$Chromosome] &
         start_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome]) ~ "arm1",
      
      # Tag as arm2 if start, end, start_ins, end_ins are all greater than cenend
      (start > cenrange$cenend[chromosome == cenrange$Chromosome] & 
         end > cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome] &
         end_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "arm2",
      
      # Tag centromer if start is greater or equal and cenend is smaller or equal
      (start >= cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end <= cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins >= cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins <= cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "centromere",
      
      TRUE ~ NA_character_  # This will leave NA if neither condition is met
    )
  ) %>%
  ungroup()

# View the updated dataframe
head(dfm)

write.table(dfm, "bedtools_closest_ouput_100kb_arms_NMF_20241003.csv", sep = ",", col.names = TRUE)
dfm <- read.csv("bedtools_closest_ouput_100kb_arms_NMF_20241003.csv", sep = ",", header = TRUE)

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

# SI

# Read in bedtools closest results (10kb)
dfm <- read.csv("bedtools_closest_ouput_100kb.csv", sep = ",", header = TRUE)
head(dfm)

# Select intersting columns
dfm <- dfm %>%
  select("chromosome", "start", "end", "rho", 
         "chromosome_ins", "start_ins", "end_ins", "no_insert", "dist_to_Cla", 
         "ind_chr", "pop")
# Create the chr_arms column, marking Chr4 with NA and applying the filter to other chromosomes
dfm <- dfm %>%
  filter(pop == "SI") %>%
  rowwise() %>%
  mutate(
    chr_arms = case_when(
      chromosome == "Chr4" ~ NA_character_,  # Mark everything from Chr4 with NA
      
      # Tag as arm1 if start, end, start_ins, end_ins are all less than censtart
      (start < cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end < cenrange$censtart[chromosome == cenrange$Chromosome] &
         start_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome]) ~ "arm1",
      
      # Tag as arm2 if start, end, start_ins, end_ins are all greater than cenend
      (start > cenrange$cenend[chromosome == cenrange$Chromosome] & 
         end > cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome] &
         end_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "arm2",
      
      # Tag centromer if start is greater or equal and cenend is smaller or equal
      (start >= cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end <= cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins >= cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins <= cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "centromere",
      
      TRUE ~ NA_character_  # This will leave NA if neither condition is met
    )
  ) %>%
  ungroup()

# View the updated dataframe
head(dfm)

write.table(dfm, "bedtools_closest_ouput_100kb_arms_SI_20241003.csv", sep = ",", col.names = TRUE)
dfm <- read.csv("bedtools_closest_ouput_100kb_arms_SI_20241003.csv", sep = ",", header = TRUE)

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

# SS

# Read in bedtools closest results (10kb)
dfm <- read.csv("bedtools_closest_ouput_100kb.csv", sep = ",", header = TRUE)
head(dfm)

# Select intersting columns
dfm <- dfm %>%
  select("chromosome", "start", "end", "rho", 
         "chromosome_ins", "start_ins", "end_ins", "no_insert", "dist_to_Cla", 
         "ind_chr", "pop")
# Create the chr_arms column, marking Chr4 with NA and applying the filter to other chromosomes
dfm <- dfm %>%
  filter(pop == "SS") %>%
  rowwise() %>%
  mutate(
    chr_arms = case_when(
      chromosome == "Chr4" ~ NA_character_,  # Mark everything from Chr4 with NA
      
      # Tag as arm1 if start, end, start_ins, end_ins are all less than censtart
      (start < cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end < cenrange$censtart[chromosome == cenrange$Chromosome] &
         start_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins < cenrange$censtart[chromosome_ins == cenrange$Chromosome]) ~ "arm1",
      
      # Tag as arm2 if start, end, start_ins, end_ins are all greater than cenend
      (start > cenrange$cenend[chromosome == cenrange$Chromosome] & 
         end > cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome] &
         end_ins > cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "arm2",
      
      # Tag centromer if start is greater or equal and cenend is smaller or equal
      (start >= cenrange$censtart[chromosome == cenrange$Chromosome] & 
         end <= cenrange$cenend[chromosome == cenrange$Chromosome] &
         start_ins >= cenrange$censtart[chromosome_ins == cenrange$Chromosome] &
         end_ins <= cenrange$cenend[chromosome_ins == cenrange$Chromosome]) ~ "centromere",
      
      TRUE ~ NA_character_  # This will leave NA if neither condition is met
    )
  ) %>%
  ungroup()

# View the updated dataframe
head(dfm)

write.table(dfm, "bedtools_closest_ouput_100kb_arms_SS_20241003.csv", sep = ",", col.names = TRUE)
dfm <- read.csv("bedtools_closest_ouput_100kb_arms_SS_20241003.csv", sep = ",", header = TRUE)

#---------------------------------------------------------------------------------------------------------