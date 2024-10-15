###########################################################
#                                                         #
# Statistics recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           REVISED                       #
#                         Juli 2024                       #
###########################################################

## Clean environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
#library(chromoMap)
#library(scales)
#library(ggplot2)

# Set directory
#setwd("/home/alle/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/")
setwd("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/")


getwd()

#-----------------------------------
# EXCLUDE CENTROMERE REGIONS
#-----------------------------------
# List of filenames
filenames <- c(
  'RepeatOBserver/run01/Summary_output/output_data/Criparius_H0-AT_Chr1_Shannon_div.txt',
  'RepeatOBserver/run01/Summary_output/output_data/Criparius_H0-AT_Chr2_Shannon_div.txt',
  'RepeatOBserver/run01/Summary_output/output_data/Criparius_H0-AT_Chr3_Shannon_div.txt',
  'RepeatOBserver/run01/Summary_output/output_data/Criparius_H0-AT_Chr4_Shannon_div.txt'
)

# Function to read file and add a new column
read_and_add_column <- function(filepath) {
  # Extract the filename from the filepath
  filename <- basename(filepath)
  
  # Extract chromosome info from filename
  chromosome <- strsplit(filename, "_")[[1]][3]
  
  # Read the file into a dataframe
  df <- read.table(filepath, header=FALSE, sep=" ")
  
  # Add a new column with the chromosome info
  df <- df %>%
    mutate(Chromosome = chromosome) %>%
    rename(Genome_position = V1,
           Shannon_div = V2)
  
  return(df)
}

# Read all files and combine into a single dataframe
Shannon_div_total <- do.call(rbind, lapply(filenames, read_and_add_column))

# Plot it
Shannon_div_total %>%
  ggplot(aes(x=Genome_position, y=Shannon_div)) +
  geom_point() +
  facet_wrap(~Chromosome, scales = "free")+
  theme_classic()


# Roll mean for Shannon div data
Shannon_div_total$Shannon_roll_mean <- zoo::rollapply(Shannon_div_total$Shannon_div, width = 100, FUN=mean, fill = NA, partial=(100/2))


# Plot it
Shannon_div_total %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean)) +
  geom_point() +
  facet_wrap(~Chromosome, scales = "free")+
  theme_classic()

# Define Peak regions!
# Sample data 
summary(Shannon_div_total)
summary(Shannon_div_total$Shannon_roll_mean)

h_mean <- Shannon_div_total %>% 
  #group_by(Chromosome)  %>% 
  summarise(N = n(),
            mean_h = mean(Shannon_roll_mean),
            SE.low = mean_h - (sd(Shannon_roll_mean))/sqrt(N),
            SE.high = mean_h + (sd(Shannon_roll_mean))/sqrt(N))
head(h_mean)

# Define a threshold
threshold <- 4.8


# Create a new column to mark peak points
Shannon_div_total <- Shannon_div_total%>%
  mutate(is_peak = Shannon_roll_mean < threshold)

# Plot it
Shannon_div_total %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5) +
  facet_wrap(~Chromosome, scales = "free")+
  theme_classic()

# Centromere prediction of RepOBserver based on Shannon Div
cenpred <- read.table("RepeatOBserver/run01/Summary_output/Shannon_div/Criparius_H0-AT_Centromere_summary_Shannon_500.txt")
cenpred <- cenpred %>%
  rename(
    Analysis = V1,
    Centromere_Prediciton = V2,
    Chr_Length = V3,
    Species = V4,
    Chromosome = V5)
  
cenpred

# Centromere prediction of RepOBserver based on Histogram
cenhist <- read.table("RepeatOBserver/run01/Summary_output/histograms/Criparius_H0-AT_Centromere_histograms_summary.txt")
cenhist <- cenhist %>%
  rename(
    Species = V1,
    Chromosome = V2,
    Centromere_Prediciton_hist = V3,
    Chr_Length = V4)

cenhist

# Join datafrrames
Shannon_div_total <- Shannon_div_total %>%
  left_join(cenpred, by = "Chromosome")

Shannon_div_total <- Shannon_div_total %>%
  left_join(cenhist, by = "Chromosome")

head(Shannon_div_total)

# Mark it in the plot
Shannon_div_total %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5) +
  geom_vline(aes(xintercept = Centromere_Prediciton), color = "red", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = Centromere_Prediciton-2e6), color = "grey", linetype = "dashed", size = 0.5) + 
  geom_vline(aes(xintercept = Centromere_Prediciton+2e6), color = "grey", linetype = "dashed", size = 0.5) +
  facet_wrap(~Chromosome, scales = "free") +
  theme_classic() 

# Extract centromer regions
peak <- Shannon_div_total %>%
  filter(is_peak == TRUE) 

# Extract first value of the first peak that is <4.8 and extract last value of the second peak that is <4.8
# Chr1
peak_chr1 <- peak %>%
  filter(Chromosome == "Chr1") 

peak_chr1 <- peak_chr1 %>%
  filter(Genome_position > 24007501) %>%
  filter(Genome_position < 35892501) # a bit narrower, excluding a small extra peak: 34727501

peak_chr1 %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5)
  

# Chr2 
peak_chr2 <- peak %>%
  filter(Chromosome == "Chr2") %>%
  filter(Genome_position > 27232501) %>%
  filter(Genome_position < 36202501) 

peak_chr2 %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5)

# Chr3
peak_chr3 <- peak %>%
  filter(Chromosome == "Chr3") %>%
  filter(Genome_position > 20902501) %>%
  filter(Genome_position < 27762501) 

peak_chr3 %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5)

# Chr4 --> signal not strong enough to be sure that it is centromeres
peak_chr4 <- peak %>%
  filter(Chromosome == "Chr4") %>%
  filter(Genome_position > 7722501) %>%
  filter(Genome_position < 9572501) 

peak_chr4 %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5)

# each chromosome
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

# Add to big plot
cenrange <- data_frame(
  censtart = c(24007501, 27232501, 20902501, 7722501),
  cenend = c(35892501, 36202501, 27762501, 9572501),
  Chromosome = chr
)

Shannon_div_total <- Shannon_div_total %>%
  left_join(cenrange, by = "Chromosome")
head(Shannon_div_total)

# Mark it in the plot
Shannon_div_total %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5) +
  geom_vline(aes(xintercept = Centromere_Prediciton), color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = Centromere_Prediciton_hist), color = "grey", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = censtart), color = "red", linetype = "dashed", size = 0.5) + 
  geom_vline(aes(xintercept = cenend), color = "red", linetype = "dashed", size = 0.5) +
  facet_wrap(~Chromosome, scales = "free") +
  theme_classic()

# Save a table of the range
cenrange <- cenrange  %>%
  mutate(length = cenend - censtart)

print(cenrange)
# Centromere Range:
# censtart   cenend Chromosome   length
# 24007501 35892501 Chr1       11885000
# 27232501 36202501 Chr2        8970000
# 20902501 27762501 Chr3        6860000
#  7722501  9572501 Chr4        1850000

# save as file
write.table(cenrange, file = "CentromereRange.csv", sep = ",", col.names = TRUE)



# Exclude centromere range from analysis!
# We will use the complete Chr4 for the anlyisis because the cenromere position is uncertain
# Filter out the rows that fall within the specified centromere ranges
Shannon_div_filtered <- Shannon_div_total %>%
  rowwise() %>%
  filter(
    Chromosome == "Chr4" | 
      if_else(
        Chromosome %in% cenrange$Chromosome,
        !(Genome_position >= cenrange$censtart[Chromosome == cenrange$Chromosome] &
            Genome_position <= cenrange$cenend[Chromosome == cenrange$Chromosome]),
        TRUE
      )
  ) %>%
  ungroup()

# Check the filtered dataframe
head(Shannon_div_filtered)

# Mark it in the plot
Shannon_div_filtered %>%
  ggplot(aes(x=Genome_position, y=Shannon_roll_mean, color = is_peak)) +
  geom_point(size=0.5) +
  geom_vline(aes(xintercept = Centromere_Prediciton), color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = Centromere_Prediciton_hist), color = "grey", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = censtart), color = "red", linetype = "dashed", size = 0.5) + 
  geom_vline(aes(xintercept = cenend), color = "red", linetype = "dashed", size = 0.5) +
  facet_wrap(~Chromosome, scales = "free") +
  theme_classic()


# save as file
write.table(Shannon_div_filtered, file = "Shannon_div_roll_mean_filtered_regions_used_for_decay.csv", sep = ",", col.names = TRUE)
write.table(Shannon_div_total, file = "Shannon_div_total_roll_mean.csv", sep = ",", col.names = TRUE)
