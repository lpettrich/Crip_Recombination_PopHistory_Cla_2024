###########################################################
#                                                         #
# Statistics recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           January 2024                  #
###########################################################

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
df1 <- read.csv("bedtools_closest_ouput_100kb_arms_MF.csv", sep = ",", header = TRUE)
df2 <- read.csv("bedtools_closest_ouput_100kb_arms_MG.csv", sep = ",", header = TRUE)
df3 <- read.csv("bedtools_closest_ouput_100kb_arms_NMF.csv", sep = ",", header = TRUE)
df4 <- read.csv("bedtools_closest_ouput_100kb_arms_SI.csv", sep = ",", header = TRUE)
df5 <- read.csv("bedtools_closest_ouput_100kb_arms_SS.csv", sep = ",", header = TRUE)

dfm <- rbind(df1, df2, df3, df4, df5)

head(dfm)

dfm <- dfm %>% 
  mutate(insert_size = end_ins - start_ins)

# Filter to only have <500 bp insert size
dfm <- dfm %>%
  filter(insert_size >= 500)

#-----------------------------------
# DATA ANALYSIS
#-----------------------------------
# Select intersting columns
dfm <- dfm %>%
  select("chromosome", "start", "end", "rho", 
         "no_insert", "dist_to_Cla", "ind_chr", "pop", "chr_arms")


head(dfm)

# NOW WE CALCULATE THE MEAN FOR EACH WINDOW OF EACH POPULATION
dfm_mean <- dfm %>%
  group_by(chromosome, start, end, pop, ind_chr) %>%
  summarise(
    N = n(),
    mean_rho = (rho),
    mean_dist_to_Cla = (dist_to_Cla),
    sd.rho = sd(rho),
    sd.dist.Cla = sd(dist_to_Cla))

head(dfm_mean)

## dfm_mean is NOT the data frame used for plotting

#######################################################################################
# DECAY PLOTS
#######################################################################################

# First look
dfm %>% 
  filter(chromosome == "Chr1") %>%
  filter(chr_arms == "arm1") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = pop)) + geom_point(size=.5) 

# First look
dfm %>% 
  filter(chromosome == "Chr1") %>%
  filter(chr_arms == "arm2") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = pop)) + geom_point(size=.5) 

# First look
dfm %>% 
  filter(chromosome == "Chr2") %>%
  filter(chr_arms == "arm1") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = pop)) + geom_point(size=.5) 

# First look
dfm %>% 
  filter(chromosome == "Chr2") %>%
  filter(chr_arms == "arm2") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = pop)) + geom_point(size=.5)

# First look
dfm %>% 
  filter(chromosome == "Chr3") %>%
  filter(chr_arms == "arm1") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = pop)) + geom_point(size=.5) 

# First look
dfm %>% 
  filter(chromosome == "Chr3") %>%
  filter(chr_arms == "arm2") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = pop)) + geom_point(size=.5)

#####################
# Decay Plots   #
#####################
#-------------------------------------------------------------------
# According to https://github.com/jarobin/CaliforniaCondorGenome2021/blob/main/gc_cpg_recombinationrate.sh
#-------------------------------------------------------------------

#####################################################################################
#####################################################################################
#-------------------------------------------------------------------
# Remove incomplete windows
#-------------------------------------------------------------------
window_size <- 100e3
dfm <- dfm %>%
  mutate(wnd_size = end-start)

unique(dfm$wnd_size) # window size should be only 100,000

dfm <- dfm %>% 
  filter(wnd_size == window_size)

unique(dfm$wnd_size) # window size should be only 100,000
# it worked!

#-------------------------------------------------------------------
# Merge info on chr and arms
#-------------------------------------------------------------------
head(dfm)

dfm <- dfm %>%
  mutate(chromosome_arm = paste(chromosome, chr_arms, sep = "_")) %>%
  mutate(chromosome_arm_pop = paste(chromosome, chr_arms, pop, sep = "_"))

#-------------------------------------------------------------------
# Check data
#-------------------------------------------------------------------
head(dfm)
tail(dfm)
unique(dfm$chromosome_arm)


max(dfm$dist_to_Cla)
min(dfm$dist_to_Cla)

summary(dfm$dist_to_Cla)
boxplot(dfm$dist_to_Cla)


head(dfm)
#-------------------------------------------------------------------------------------------------------------
################################################################
#                   PEARSON or SPEARMAN?                       #
################################################################
# Coding for statistics was improved with ChatGPT
pearson_result <- cor.test(dfm$rho, dfm$dist_to_Cla, method = "pearson")
print(pearson_result)


#spearman_result <- cor.test(dfm$rho, dfm$dist_to_Cla, method = "spearman")
#print(spearman_result)

################################################################
#                   PEARSON CORRELTAION PLOT                   #
################################################################
# Now we want to test every population separately (biological replicates)
# Every chromosome will be tested separately


# Split dfm_mean by Pop and Chr
# Create an empty list to store the split dataframes
split_dataframes <- list()

pop <- unique(dfm$chromosome_arm_pop)

# Use a for loop to split the dataframe by 'Group' and store in the list
for (group_value in pop) {
  subset_df <- subset(dfm, chromosome_arm_pop == group_value)
  
  # Check if subset_df is empty and assign NULL if true
  if (nrow(subset_df) == 0) {
    split_dataframes[[group_value]] <- NULL
  } else {
    split_dataframes[[group_value]] <- subset_df
  }
}

# Now split_dataframes may contain NULL for empty subsets


# Store the split dataframes as separate dataframes with distinct variable names
for (group_value in pop) {
  # Define variable names for the new dataframes
  new_df_name <- paste("df_", group_value, sep = "")
  # Assign the dataframe from the list to the variable name
  assign(new_df_name, split_dataframes[[group_value]])
}

# Create an empty list to store correlation results
correlation_results <- list()

# PEARSON CORRELATION
# 1st way of saving
# Perform Pearson correlation tests for each dataframe
for (df_name in names(split_dataframes)) {
  # Get the current dataframe
  current_df <- split_dataframes[[df_name]]
  # Calculate Pearson correlation
  correlation_result <- cor.test(current_df$rho, current_df$dist_to_Cla, method = "pearson")
  # Store the result in the named list with dataframe names as keys
  correlation_results[[df_name]] <- correlation_result
}

# Print the correlation results for each dataframe
for (df_name in names(correlation_results)) {
  cat("Correlation Test for DataFrame:", df_name, "\n")
  print(correlation_results[[df_name]])
}


# 2nd way of saving
# Create an empty dataframe to store results
results_df <- data.frame(Dataframe = character(), Correlation = numeric(), pvalue = numeric())

# Perform Pearson correlation tests for each dataframe
for (df_name in names(split_dataframes)) {
  # Get the current dataframe
  current_df <- split_dataframes[[df_name]]
  # Calculate Pearson correlation
  correlation_result <- cor.test(current_df$rho, current_df$dist_to_Cla, method = "pearson")
  # Create a new row with dataframe name and correlation result
  new_row <- data.frame(Dataframe = df_name, Correlation = correlation_result$estimate, 
                        pvalue = correlation_result$p.value)
  # Add the row to the results dataframe
  results_df <- rbind(results_df, new_row)
}

# Print the results dataframe
print(results_df)

str(results_df)

# Add asterisks to indicate significance levels
results_df$Significance <- ifelse(results_df$pvalue < 0.001, "***",
                                  ifelse(results_df$pvalue < 0.01, "**",
                                         ifelse(results_df$pvalue < 0.05, "*",
                                                "")))



# Plot correlation
# Define a theme
mytheme <- theme(axis.text.x = element_text(size = 20, color = "black"),
                 axis.text.y = element_text(size = 20, color = "black"),
                 axis.title.y = element_text(size = 20,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 20,face = "bold", color = "black"),
                 title = element_text(size = 20, color = "black"),
                 text = element_text(size = 20, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("white"),
                 panel.grid.minor = element_line("white"),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.text = element_text(size = 20, color = "black"),
                 legend.title = element_text(size = 20, color = "black"), 
                 legend.background = element_rect(fill="white"),
                 legend.position="right") 

level_order <- c("Chr1_arm1_MF", "Chr1_arm1_MG", "Chr1_arm1_NMF", "Chr1_arm1_SI", "Chr1_arm1_SS",
                 "Chr1_arm2_MF", "Chr1_arm2_MG", "Chr1_arm2_NMF", "Chr1_arm2_SI", "Chr1_arm2_SS",
                 "Chr1_NA_MF", "Chr1_NA_MG", "Chr1_NA_NMF", "Chr1_NA_SI", "Chr1_NA_SS",
                 "Chr2_arm1_MF", "Chr2_arm1_MG", "Chr2_arm1_NMF", "Chr2_arm1_SI", "Chr2_arm1_SS",
                 "Chr2_arm2_MF", "Chr2_arm2_MG", "Chr2_arm2_NMF", "Chr2_arm2_SI", "Chr2_arm2_SS",
                 "Chr2_NA_MF", "Chr2_NA_MG", "Chr2_NA_NMF", "Chr2_NA_SI", "Chr2_NA_SS",
                 "Chr3_arm1_MF", "Chr3_arm1_MG", "Chr3_arm1_NMF", "Chr3_arm1_SI", "Chr3_arm1_SS",
                 "Chr3_arm2_MF", "Chr3_arm2_MG", "Chr3_arm2_NMF", "Chr3_arm2_SI", "Chr3_arm2_SS",
                 "Chr3_NA_MF", "Chr3_NA_MG", "Chr3_NA_NMF", "Chr3_NA_SI", "Chr3_NA_SS",
                 "Chr4_NA_MF", "Chr4_NA_MG", "Chr4_NA_SI", "Chr4_NA_SS")


ggplot(results_df, aes(x = Dataframe, y = Correlation)) +
  geom_point() +
  scale_y_continuous(limits = c(-1,1)) +
  scale_x_discrete(limits = level_order) +
  geom_text(aes(label = paste(round(Correlation, 2),  "\n", results_df$Significance)), vjust = -0.5) +  # Add correlation labels
  labs(title = "Correlation Plot", x = "Dataframe", y = "Correlation") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# Plot as heatmap
results_df

results_df <- results_df %>%
  separate(Dataframe, into = c("chr", "arm", "pop"), sep = "_", remove = FALSE) %>%
  mutate(chr_arm = paste(chr, arm, sep = "_"))  # Combine Chr and arm into a new column

# Save results of Pearson correlation
write.table(results_df, 
            "Pearson_Correlation_mean_rho_vs_distance-to-Cla_run12_more500bp.txt", sep = ";")

# Create a correlation matrix
cor_matrix <- matrix(results_df$Correlation, nrow = 1)

# Create a heatmap-like correlation plot
p <- ggplot() +
  geom_tile(data = results_df, aes(x = factor(pop), y = 1, fill = Correlation)) +  # Refactor Dataframe for each facet
  geom_text(data = results_df, aes(x = factor(pop), y = 1, label = paste(round(Correlation, 2), "\n", results_df$Significance)), 
            vjust = 1, hjust = 0.5, size = 3) +
  scale_fill_gradient2(low = "blue3", mid = "ivory", high = "red3", limits = c(-1, 1)) +  # Define color scale
  labs(x = "", y = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(~ chr_arm, scales = "free_x")  # Facet based on chromosome, allowing x-axis to vary per facet


p
#ggsave("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/png-files/heatmap.png", 
#       plot = last_plot(), width = 1000, height = 50, units = "px")

library(cowplot)
save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_more500bp.png", p, nrow = 1, ncol = 2, base_asp = 3.8, dpi = 500, bg = "white")

save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_more500bp.svg", p, nrow = 1, ncol = 2, base_asp = 3.8, dpi = 500, bg = "white")


# Drop values were arm is NA
results_df2 <- results_df %>%
  filter(arm != "NA")

p <- ggplot() +
  geom_tile(data = results_df2, aes(x = factor(pop), y = 1, fill = Correlation)) +  # Refactor Dataframe for each facet
  geom_text(data = results_df2, aes(x = factor(pop), y = 1, label = paste(round(Correlation, 2), "\n", results_df2$Significance)), 
            vjust = 1, hjust = 0.5, size = 7) +
  scale_fill_gradient2(low = "blue3", mid = "ivory", high = "red3", limits = c(-1, 1)) +  # Define color scale
  labs(x = "", y = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(~ chr_arm, scales = "fixed", ncol = 2)  # Facet based on chromosome, allowing x-axis to vary per facet


p

save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_more500bp_noNAarm.png", p, nrow = 1.8, ncol = 2, base_asp = 3.8, dpi = 500, bg = "white")

save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_more500bp_noNAarm.svg", p, nrow = 3, ncol = 3, base_asp = 3.8, dpi = 500, bg = "white")

