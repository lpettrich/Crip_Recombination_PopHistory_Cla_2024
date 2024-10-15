# Which insertions are on arm1, arm2 and the centromeres?
# Clean environment
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

# Read in bedtools closest results (10kb)
dfm <- read.csv("bedtools_closest_ouput_100kb.csv", sep = ",", header = TRUE)
head(dfm)

# Select intersting columns
dfm <- dfm %>%
  select("chromosome", "start", "end", "rho", 
         "chromosome_ins", "start_ins", "end_ins", "no_insert", "dist_to_Cla", 
         "ind_chr", "pop")


# Read in centromere ranges
cenrange <- read.csv("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/CentromereRange.csv", sep = ",", header = TRUE)
head(cenrange)

#-----------------------------------
# SPLIT CHROMOSOME ARMS
#-----------------------------------
# in the supplement directory
dfm <- read.table("~/recombination-map/MELT/MELT-run2/Cla1/AllPop/2-GroupAnalysis/Cla1.master_AllPop.bedgraph")

head(dfm)

Cla_colnames = c("chromosome", "start", "end", "no_insert")
dfm <- dfm %>%
  setNames(Cla_colnames)

dfm <- dfm %>% 
  mutate(insert_size = end - start)

head(dfm)


# Create the chr_arms column, marking Chr4 with complete and applying the filter to other chromosomes
# Ensure the cenrange data is merged before mutate to avoid indexing issues
dfm <- dfm %>%
  left_join(cenrange, by = c("chromosome" = "Chromosome")) %>%
  rowwise() %>%
  mutate(
    chr_arms = case_when(
      chromosome == "Chr4" ~ "complete",  # Mark everything from Chr4 with NA
      
      # Tag as arm1 if start and end are all less than censtart
      start < censtart & end < censtart ~ "arm1",
      
      # Tag as arm2 if start and end are all greater than cenend
      start > cenend & end > cenend ~ "arm2",
      
      # Tag centromere
      start >= censtart & end <= cenend ~ "centromere",
      
      TRUE ~ NA_character_  # Leave NA if neither condition is met
    )
  ) %>%
  ungroup()


#-----------------------------------
# DATA ANALYSIS
#-----------------------------------

head(dfm)


# Read in inserts
ins_sub <- read.csv("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/bedtools-closest/run12/MELT_sharedsubsetins.csv", sep = ",", header = TRUE)

ins_all <- read.csv("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/bedtools-closest/run12/MELT_sharedallins.csv", sep = ",", header = TRUE)

ins_uni <- read.csv("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/bedtools-closest/run12/MELT_uniqueins.csv", sep = ",", header = TRUE)

inserts <- rbind(ins_sub, ins_all, ins_uni)

head(inserts)


# Filter inserts of dfm
dfm <- left_join(dfm, inserts, by = c("no_insert"))

head(dfm)


# Add unique and shared flag
# Create a presence mapping based on inserts
presence_mapping <- dfm %>%
  mutate(Presence = case_when(
    !grepl(",", pops) ~ "Unique",  # Only one population is present (no commas)
    grepl(",", pops)  ~ "Shared",  # Multiple populations are present (comma-separated)
    TRUE ~ NA_character_  # Ensure NA for any other cases
  ))


# After filtering for unique and shared, the number should stay the same
length(unique(dfm$no_insert))
length(unique(inserts$no_insert))
length(unique(presence_mapping$no_insert))





#-----------------------------------------------------------------
# count how many on chromosomes and centromeres
# Coded with help of ChatGPT
# Group by 'no_insert' and count how many distinct populations each insert appears in
presence_mapping <- presence_mapping %>%
  mutate(chr_arms_m = paste(chromosome, chr_arms, sep = "_"))

head(presence_mapping)
write.table(presence_mapping, "Complete_list_on_unique_shared_on_which_arm_per_pop.csv", sep = ",", col.names = TRUE)


arms_insert_counts <- presence_mapping %>%
  group_by(chr_arms_m, Presence) %>%
  summarise(cla_count = n_distinct(unique(no_insert))) 

print(arms_insert_counts, n = 24)
sum(arms_insert_counts$cla_count)

write.table(arms_insert_counts, "MELT_countins_shared_unique_per_arm.csv", sep = ",", col.names = TRUE)


# count instertions in centromere range
sum_cent <- arms_insert_counts %>%
  filter(chr_arms_m %in% c("Chr1_centromere", "Chr2_centromere", "Chr3_centromere"))

print(sum_cent)

# count insertions on arms
sum_cent %>% 
  group_by(Presence) %>%
  reframe(sum(cla_count))


sum_arms <- arms_insert_counts %>%
  filter(chr_arms_m %in% c("Chr1_arm1", "Chr2_arm1", "Chr3_arm1",
                           "Chr1_arm2", "Chr2_arm2", "Chr3_arm2"))

print(sum_arms)

sum_arms %>% 
  group_by(Presence) %>%
  reframe(sum(cla_count))

# Chr4
arms_insert_counts %>% 
  filter(chr_arms_m == "Chr4_complete") %>% 
  group_by(Presence) %>%
  reframe(sum(cla_count))

# Total
arms_insert_counts %>% 
  group_by(Presence) %>%
  reframe(sum(cla_count))

arms_insert_counts %>%
  filter(chr_arms_m %in% c("Scaffold1_NA", "Scaffold3_NA", "Scaffold7_NA")) %>%
  group_by(Presence) %>%
  reframe(sum(cla_count))





