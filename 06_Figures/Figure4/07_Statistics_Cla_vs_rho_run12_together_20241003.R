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
library(patchwork)
library(cowplot)

# Set directory
# setwd("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")
setwd("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/bedtools-closest/run12/")


getwd()

# Read data
df1 <- read.table("Pearson_Correlation_mean_rho_vs_distance-to-Cla_run12_less500bp_20241003.txt", sep = ";")
df2 <- read.table("Pearson_Correlation_mean_rho_vs_distance-to-Cla_run12_more500bp_20241003.txt", sep = ";")

df1$filter <- "<500bp"
df2$filter <- ">=500bp"

results_df <- rbind(df1, df2)


# Define a theme
mytheme <- theme(axis.text.x = element_text(size = 16, color = "black"),
                 axis.text.y = element_text(size = 16, color = "black"),
                 axis.title.y = element_text(size = 16,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 16,face = "bold", color = "black"),
                 title = element_text(size = 16, color = "black"),
                 text = element_text(size = 16, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("white"),
                 panel.grid.minor = element_line("white"),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.text = element_text(size = 16, color = "black"),
                 legend.title = element_text(size = 16, color = "black"), 
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

# Create a combined column for facet labels
results_df$facet_label <- paste(results_df$chr, results_df$arm, sep = "\n")  # Combine Chr and arm into a single label


# Create a heatmap-like correlation plot
p <- ggplot() +
  geom_tile(data = results_df, aes(x = factor(pop), y = 1, fill = Correlation)) +  # Refactor Dataframe for each facet
  geom_text(data = results_df, aes(x = factor(pop), y = 1, label = paste(round(Correlation, 2), "\n", results_df$Significance)), 
            vjust = 1, hjust = 0.5, size = 5)  +
  scale_fill_gradient2(low = "blue3", mid = "ivory", high = "red3", limits = c(-1, 1)) +  # Define color scale
  labs(x = "", y = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(facet_label ~ filter, scales = "fixed", switch = "y")   +  # Facet based on chromosome and arm
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1))   # Angle y-axis facet labels 


p


save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_20241003.png", p, nrow = 3, ncol = 2, base_asp = 1.2, dpi = 500, bg = "white")

save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_20241003.svg", p, nrow = 3, ncol = 2, base_asp = 1.2, dpi = 500, bg = "white")




# Drop values were arm is NA
results_df$arm[results_df$chr == "Chr4"] <- "comp"

results_df$arm[results_df$arm == "centromere"] <- "cent"

results_df$facet_label <- paste(results_df$chr, results_df$arm, sep = "\n")  # Combine Chr and arm into a single label

results_df2 <- results_df %>%
  filter(arm != "NA")



# Create a heatmap-like correlation plot
p <- ggplot() +
  geom_tile(data = results_df2, aes(x = factor(pop), y = 1, fill = Correlation)) +  # Refactor Dataframe for each facet
  geom_text(data = results_df2, aes(x = factor(pop), y = 1, label = paste(round(Correlation, 2), "\n", results_df2$Significance)), 
            vjust = 1, hjust = 0.5, size = 5) +
  scale_fill_gradient2(low = "blue3", mid = "ivory", high = "red3", limits = c(-1, 1)) +  # Define color scale
  labs(x = "", y = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(facet_label ~ filter, scales = "fixed", switch = "y")   +  # Facet based on chromosome and arm
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1))  # Angle y-axis facet labels


p

save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_noNAarm_20241003.png", p, nrow = 3, ncol = 2, base_asp = 1.2, dpi = 500, bg = "white")

save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/heatmap_run12_noNAarm_20241003.svg", p, nrow = 3, ncol = 2, base_asp = 1.2, dpi = 500, bg = "white")

