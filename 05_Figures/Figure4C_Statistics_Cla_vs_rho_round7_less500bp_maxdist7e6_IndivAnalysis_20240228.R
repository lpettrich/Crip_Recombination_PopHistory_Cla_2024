###########################################################
#                                                         #
# Statistics recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           January 2024                  #
###########################################################

## Clean environment
rm(list = ls())

#3 Load libraries
library(tidyverse)
library(dplyr)
#library(chromoMap)
#library(scales)
library(ggplot2)

# Set directory
#setwd("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")
setwd("~/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")
getwd()


#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
# Here you define vectors for every data table you have
# This info will be added a column later
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

# Here you tell R the path where it should look for dataframes
#common_path = "C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv/02_less500bp//"
common_path = "~/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv/02_less500bp/"

# Read in bedtools closest results
# These dataframes will be used for analysis
# You add each population seperately
## MF
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "\\bMF_*", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_1

data_lst_1 <- lapply(files_to_read_1, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})


# Read in bedtools closest results
## MG
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "MG_*", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_2

data_lst_2 <- lapply(files_to_read_2, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## NMF
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "NMF_*", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_3

data_lst_3 <- lapply(files_to_read_3, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## SI
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "SI_*", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_4

data_lst_4 <- lapply(files_to_read_4, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## SS
files_to_read_5 = list.files(
  path = common_path,        # directory to search within
  pattern = "SS_*", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_5

data_lst_5 <- lapply(files_to_read_5, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

#-------------------------------------------------------------------
# Rename columns
#-------------------------------------------------------------------
# Rename data list names based on files_to_read
names(data_lst_1) <- ind[1:16]
names(data_lst_2) <- ind[17:32]
names(data_lst_3) <- ind[33:48]
names(data_lst_4) <- ind[49:64]
names(data_lst_5) <- ind[65:80]

# Check if renaming was done correctly
print(names(data_lst_1))
print(names(data_lst_2))
print(names(data_lst_3))
print(names(data_lst_4))
print(names(data_lst_5))

# Drop dataframes with NULL enttries
for (i in 1:5) {
  obj_name <- paste0("data_lst_", i)
  assign(obj_name, Filter(function(df) !is.null(df), get(obj_name)))
}


# Create a list of all the data frames
data_lst <- list(data_lst_1, data_lst_2, data_lst_3, data_lst_4, data_lst_5)

# Add a "ind" column to each data frame with the corresponding 'ind' value
for (j in 1:5) {
  for (i in 1:length(data_lst[[j]])) {
    data_lst[[j]][[i]]$ind <- names(data_lst[[j]])[i]
  }
}

# Create separate data lists (data_lst_1, data_lst_2, etc.) from data_lst
for (j in 1:5) {
  assign(paste0("data_lst_", j), data_lst[[j]])
}

# Remove all dataframes with no Cla-Element in adjacent regions (value -1 in column v6 - v8)
# Alternatively, if you want to remove the entire dataframe from the list
data_lst_1 <- data_lst_1[sapply(data_lst_1, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_2 <- data_lst_2[sapply(data_lst_2, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_3 <- data_lst_3[sapply(data_lst_3, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_4 <- data_lst_4[sapply(data_lst_4, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_5 <- data_lst_5[sapply(data_lst_5, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]



# Mutate dataframe 
df1 <- bind_rows(data_lst_1)
df2 <- bind_rows(data_lst_2)
df3 <- bind_rows(data_lst_3)
df4 <- bind_rows(data_lst_4)
df5 <- bind_rows(data_lst_5)

# Add column with Pop name
df1$pop <- "MF"
df2$pop <- "MG"
df3$pop <- "NMF"
df4$pop <- "SI"
df5$pop <- "SS"

# Bind data frames together
dfm <- rbind(df1, df2, df3, df4, df5)
head(dfm)
tail(dfm)
# columns means:
# "chromosome", "start", "end", "rho", 
# "chromosome", "start", "end", "insertions", "strand"
# "distance", "ind", "pop"

dfm <- dfm %>%
  select(V1:4, V8, V11, ind, pop)

# Rename columns in data frame
# https://stackoverflow.com/questions/28648513/changing-column-names-in-a-list-of-data-frames-in-r
new_colnames <- c("chromosome", "start", "end", "rho", 
                  "no_insert", "dist_to_Cla", "Ind_Chr", "Pop")

dfm <- dfm %>%
  rename_with(~ new_colnames)
head(dfm)

# NOW WE NEED TO CALCULATE THE MEAN FOR EACH WINDOW OF EACH POPULATION
dfm_mean <- dfm %>%
  group_by(chromosome, start, end, Pop) %>%
  summarise(
    N = n(),
    mean_rho = (rho),
    mean_dist_to_Cla = (dist_to_Cla),
    sd.rho = sd(rho),
    sd.dist.Cla = sd(dist_to_Cla),
    Pop_Chr = paste0(Pop, "_", chromosome))

head(dfm_mean)
str(dfm_mean)
summary(dfm_mean)

## dfm_mean is data frame used for plotting

# Recombination compared to distance to Cla
# First plot to get overview
unique(dfm_mean$Pop_Chr)

test <- dfm_mean %>% 
  filter(Pop_Chr == "MF_Chr1")

str(test)
plot(test$mean_rho ~ test$mean_dist_to_Cla)

dfm_mean %>% 
  #filter(Pop_Chr == "SS_Chr1") %>%
  filter(chromosome == "Chr1") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop_Chr)) + geom_point(size=.5) +
  scale_y_continuous(limits = c(0, max(dfm_mean$mean_rho)))


dfm_mean %>% 
  #filter(Pop_Chr == "MF_Chr1") %>%
  filter(chromosome == "Chr2") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop_Chr)) + geom_point(size=.5)

dfm_mean %>% 
  #filter(Pop_Chr == "MF_Chr1") %>%
  filter(chromosome == "Chr3") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop_Chr)) + geom_point(size=.5)

dfm_mean %>% 
  #filter(Pop_Chr == "MF_Chr1") %>%
  filter(chromosome == "Chr4") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop_Chr)) + geom_point(size=.5)

# Not really decisive
# Make it better!
# According to https://github.com/jarobin/CaliforniaCondorGenome2021/blob/main/gc_cpg_recombinationrate.sh

# Remove incomplete windows
dfm_mean <- dfm_mean %>%
  mutate(wnd_size = end-start)

unique(dfm_mean$wnd_size) # window size should be only 10,000

window_size = 10e3

dfm_mean <- dfm_mean %>% 
  filter(wnd_size == window_size)

# Remove values that are too distant
# read-in chromosome and scaffold data as comparison
#chr.data <- read.table("/home/alle/recombination-map/MELT/crip4.0_genome_chr_length.csv", sep = ";", header = TRUE)
#chr.data

# Summary of Cla-distance
summary(dfm_mean$mean_dist_to_Cla) 
boxplot(dfm_mean$mean_dist_to_Cla)
quantile(dfm_mean$mean_dist_to_Cla, prob=c(.25, .5,.75), type=1)
# 3rd quartile is 4,890,780
# set max dist to reasonable limit

max_dist=7e6
dfm_mean <- dfm_mean %>% 
  filter(mean_dist_to_Cla <= max_dist)
summary(dfm_mean$mean_dist_to_Cla) 
boxplot(dfm_mean$mean_dist_to_Cla)



#-------------------------------------------------------------------------------------------------------------
################################################################
#                   PEARSON or SPEARMAN?                       #
################################################################
# Coding for statistics was improved with ChatGPT
# Assuming your data frames are named 'gene_data' and 'recombination_data'
pearson_result <- cor.test(dfm_mean$mean_rho, dfm_mean$mean_dist_to_Cla, method = "pearson")
print(pearson_result)

# Assuming your data frames are named 'gene_data' and 'recombination_data'
spearman_result <- cor.test(dfm_mean$mean_rho, dfm_mean$mean_dist_to_Cla, method = "spearman")
print(spearman_result)

################################################################
#                   PEARSON CORRELTAION PLOT                   #
################################################################
# Now we want to test every population separately (biological replicates)
# Every chromosome will be tested separately

# Calculate mean per window and population
dfm_mean

# Split dfm_mean by Pop and Chr
# Create an empty list to store the split dataframes
split_dataframes <- list()

# Use a for loop to split the dataframe by 'Group' and store in the list
for (group_value in pop) {
  subset_df <- subset(dfm_mean, Pop_Chr == group_value)
  
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
  correlation_result <- cor.test(current_df$mean_rho, current_df$mean_dist_to_Cla, method = "pearson")
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
  correlation_result <- cor.test(current_df$mean_rho, current_df$mean_dist_to_Cla, method = "pearson")
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

# Save results of Pearson correlation
write.table(results_df, 
            "Pearson_Correlation_mean_rho_vs_distance-to-Cla.txt", sep = ";")

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

level_order <- c(
  "MG_Chr1", "NMF_Chr1", "MF_Chr1", "SI_Chr1", "SS_Chr1",
  "MG_Chr2", "NMF_Chr2", "MF_Chr2", "SI_Chr2", "SS_Chr2",
  "MG_Chr3", "NMF_Chr3", "MF_Chr3", "SI_Chr3", "SS_Chr3",
  "MG_Chr4", "NMF_Chr4", "MF_Chr4", "SI_Chr4", "SS_Chr4"
)


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
# Create a correlation matrix
cor_matrix <- matrix(results_df$Correlation, nrow = 1)

# Create a heatmap-like correlation plot
p <- ggplot() +
  geom_tile(data = results_df, aes(x = Dataframe, y = 1, fill = Correlation)) +
  geom_text(data = results_df, aes(x = Dataframe, y = 1, label = paste(round(Correlation, 2), "\n", results_df$Significance)), 
            vjust = 1, hjust = 0.5, size = 6) +
  #geom_text(aes(label = "NA"), 
  #          x = corrNA$Dataframe, y = 1, 
  #          vjust = 1, hjust = 0.5) +
  scale_fill_gradient2(low = "blue3", mid = "ivory", high = "red3", limits = c(-1, 1)) +  # Define color scale
  scale_x_discrete(limits = level_order, expand = c(0.075,0.075)) +
  labs(x = "", y = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  # Rotate x-axis labels for readability

p
#ggsave("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/png-files/heatmap.png", 
#       plot = last_plot(), width = 1000, height = 50, units = "px")

library(cowplot)
save_plot("../../../Plots/png-files/heatmap7e6_round7_less500bp.png", p, nrow = 1, ncol = 2, base_asp = 3.8, dpi = 500, bg = "white")

save_plot("../../../Plots/svg-files/heatmap7e6_round7_less500bp.svg", p, nrow = 1, ncol = 2, base_asp = 3.8, dpi = 500, bg = "white")

