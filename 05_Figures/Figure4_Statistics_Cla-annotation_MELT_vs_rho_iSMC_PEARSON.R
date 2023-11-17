###########################################################
#                                                         #
# Statistics recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           August 2023                   #
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
setwd("/home/alle/recombination-map")
getwd()

#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
pop <- c("MF_Chr1", "MF_Chr2", "MF_Chr3", "MF_Chr4", 
         "MG_Chr1", "MG_Chr2", "MG_Chr3", "MG_Chr4", 
         "NMF_Chr1", "NMF_Chr2", "NMF_Chr3", "NMF_Chr4", 
         "SI_Chr1", "SI_Chr2", "SI_Chr3", "SI_Chr4", 
         "SS_Chr1", "SS_Chr2", "SS_Chr3", "SS_Chr4")
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

common_path = "/home/alle/recombination-map/comparison-Cla-rho/01-bedtools_closest"

# Read in bedtools closest results
## MF
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "\\bMF_*", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_1

data_lst_1 = lapply(files_to_read_1, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## MG
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "MG_*", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_2

data_lst_2 = lapply(files_to_read_2, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## NMF
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "NMF_*", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_3

data_lst_3 = lapply(files_to_read_3, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## SI
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "SI_*", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_4

data_lst_4 = lapply(files_to_read_4, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## SS
files_to_read_5 = list.files(
  path = common_path,        # directory to search within
  pattern = "SS_*", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_5

data_lst_5 = lapply(files_to_read_5, read.table, header = FALSE)  # read all the matching files

# V10 is autocorrected to TRUE instead of keeping the T in data_lst_2 - 5
# Define data types and read in the data again
# Assuming data_lst_1 is a list of data frames
# data_lst_1 has the correct data types
data_types <- sapply(data_lst_1, function(df) sapply(df, class))
data_types

data_lst_2 = lapply(files_to_read_2, read.table, header = FALSE, colClasses = data_types[1:19])
data_lst_3 = lapply(files_to_read_3, read.table, header = FALSE, colClasses = data_types[1:19])
data_lst_4 = lapply(files_to_read_4, read.table, header = FALSE, colClasses = data_types[1:19])
data_lst_5 = lapply(files_to_read_5, read.table, header = FALSE, colClasses = data_types[1:19])


#-------------------------------------------------------------------
# Rename columns
#-------------------------------------------------------------------
# Rename data list names based on files_to_read
names(data_lst_1) <- pop[1:4]
names(data_lst_2) <- pop[5:8]
names(data_lst_3) <- pop[9:12]
names(data_lst_4) <- pop[13:16]
names(data_lst_5) <- pop[17:20]

# Check if renaming was done correctly
print(names(data_lst_1))
print(names(data_lst_2))
print(names(data_lst_3))
print(names(data_lst_4))
print(names(data_lst_5))

# Create a list of data frames
data_lst <- list(data_lst_1, data_lst_2, data_lst_3, data_lst_4, data_lst_5)

# Add a "pop" column and label data frames
# Add a "pop" column to each data frame with the corresponding 'pop' value
for (j in 1:5) {
  for (i in 1:length(data_lst[[j]])) {
    data_lst[[j]][[i]]$pop <- names(data_lst[[j]])[i]
  }
}

# Create separate data lists (data_lst_1, data_lst_2, etc.) from data_lst
for (j in 1:5) {
  assign(paste0("data_lst_", j), data_lst[[j]])
}


# Only column 1 - 4 and 19 and 20 relevant for us ("chromosome", "start", "end", "rho", "dist_to_Cla")
# Mutate dataframe and only keep these columns
df1 <- bind_rows(data_lst_1)
df2 <- bind_rows(data_lst_2)
df3 <- bind_rows(data_lst_3)
df4 <- bind_rows(data_lst_4)
df5 <- bind_rows(data_lst_5)

# SOLVED BY DEFINING DATA TYPES WHILE READING IN DATA
#---------------------------------------------------------------
# Resolve error that is caused by different data categories
#data_lst_2 <- lapply(data_lst_2, function(df) {
#  df$V10 <- as.character(df$V10)
#  return(df)
#})


# Retry bind_rows
#df2 <- bind_rows(data_lst_2)
#-----------------------------------------------------------------

# Bind data frames together
dfm <- rbind(df1, df2, df3, df4, df5)

dfm <- dfm %>%
  select(V1:4, V19, pop)

# Rename columns in data frame
# https://stackoverflow.com/questions/28648513/changing-column-names-in-a-list-of-data-frames-in-r
new_colnames <- c("chromosome", "start", "end", "rho", "dist_to_Cla", "Pop_Chr")

dfm <- dfm %>%
  rename_with(~ new_colnames)

## dfm is data frame used for plotting

# Recombination compared to distance to Cla
# First plot to get overview
unique(dfm$Pop_Chr)

test <- dfm %>% 
  filter(Pop_Chr == "MF_Chr1")

str(test)
plot(test$rho ~ test$dist_to_Cla)

dfm %>% 
  #filter(Pop_Chr == "SS_Chr1") %>%
  filter(chromosome == "Chr1") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = Pop_Chr)) + geom_point(size=.5) +
  scale_y_continuous(limits = c(0, max(dfm$rho)))


dfm %>% 
  #filter(Pop_Chr == "MF_Chr1") %>%
  filter(chromosome == "Chr2") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = Pop_Chr)) + geom_point(size=.5)

dfm %>% 
  #filter(Pop_Chr == "MF_Chr1") %>%
  filter(chromosome == "Chr3") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = Pop_Chr)) + geom_point(size=.5)

dfm %>% 
  #filter(Pop_Chr == "MF_Chr1") %>%
  filter(chromosome == "Chr4") %>%
  ggplot(aes(x = dist_to_Cla, y = rho, colour = Pop_Chr)) + geom_point(size=.5)

# Not really decisive
# Make it better!
# According to https://github.com/jarobin/CaliforniaCondorGenome2021/blob/main/gc_cpg_recombinationrate.sh

# Remove incomplete windows
dfm <- dfm %>%
  mutate(wnd_size = end-start)

unique(dfm$wnd_size) # window size should be only 10,000

window_size = 10e3

dfm <- dfm %>% 
  filter(wnd_size == window_size)

# Remove values that are too distant
# read-in chromosome and scaffold data as cmparison
chr.data <- read.table("/home/alle/recombination-map/MELT/crip4.0_genome_chr_length.csv", sep = ";", header = TRUE)
chr.data

# Summary of Cla-distance
summary(dfm$dist_to_Cla) 
boxplot(dfm$dist_to_Cla)
quantile(dfm$dist_to_Cla, prob=c(.25, .5,.75), type=1)
# 3rd quartile is 4,890,780
# set max dist to reasonable limit

max_dist=20e6
dfm <- dfm %>% 
  filter(dist_to_Cla <= max_dist)
summary(dfm$dist_to_Cla) 
boxplot(dfm$dist_to_Cla)



#-------------------------------------------------------------------------------------------------------------
################################################################
#                   PEARSON or SPEARMAN?                       #
################################################################
# Coding for statistics was improved with ChatGPT
# Assuming your data frames are named 'gene_data' and 'recombination_data'
pearson_result <- cor.test(dfm$rho, dfm$dist_to_Cla, method = "pearson")
print(pearson_result)

# Assuming your data frames are named 'gene_data' and 'recombination_data'
spearman_result <- cor.test(dfm$rho, dfm$dist_to_Cla, method = "spearman")
print(spearman_result)

################################################################
#                   PEARSON CORRELTAION PLOT                   #
################################################################
# Now we want to test every population separately (biological replicates)
# Every chromosome will be tested separately

# Split dfm by Pop and Chr
# Create an empty list to store the split dataframes
split_dataframes <- list()
# Use a for loop to split the dataframe by 'Group' and store in the list
for (group_value in pop) {
  subset_df <- subset(dfm, Pop_Chr == group_value)
  split_dataframes[[group_value]] <- subset_df
}

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

# Save results of Pearson correlation
write.table(results_df, 
            "Pearson_Correlation_rho_vs_distance-to-Cla.txt", sep = ";")

# Plot correlation
# Define a theme
mytheme <- theme(axis.text.x = element_text(size = 14, color = "black"),
                 axis.text.y = element_text(size = 16, color = "black"),
                 axis.title.y = element_text(size = 18,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 18,face = "bold", color = "black"),
                 title = element_text(size = 17, color = "black"),
                 text = element_text(size=17, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("white"),
                 panel.grid.minor = element_line("white"),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 14, color = "black"), 
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
# Create a correlation matrix
cor_matrix <- matrix(results_df$Correlation, nrow = 1)

# Create a heatmap-like correlation plot
ggplot() +
  geom_tile(data = results_df, aes(x = Dataframe, y = 1, fill = Correlation)) +
  geom_text(data = results_df, aes(x = Dataframe, y = 1, label = paste(round(Correlation, 2), "\n", results_df$Significance)), 
            vjust = 1, hjust = 0.5) +
  scale_fill_gradient2(low = "blue3", mid = "ivory", high = "red3", limits = c(-1, 1)) +  # Define color scale
  scale_x_discrete(limits = level_order, expand = c(0.075,0.075)) +
  labs(x = "", y = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  # Rotate x-axis labels for readability

# Scatterplots of data
# Create an empty list to store ggplot plots
plot_list <- list()

for (df_name in names(split_dataframes)) {
  current_df <- split_dataframes[[df_name]]
  group_value <- unique(current_df$Pop_Chr)
  
  # Create a ggplot plot
  current_plot <- current_df %>%
    ggplot(aes(x = dist_to_Cla, y = rho, colour = Pop_Chr)) +
    geom_point(size = 0.5) +
    labs(title = paste("Scatterplot for", group_value)) +
    scale_y_continuous(limits = c(0, max(dfm$rho))) +
    mytheme
  
  # Add the plot to the list with a unique name
  plot_name <- paste("plot_", group_value, sep = "")
  plot_list[[plot_name]] <- current_plot
}

for (plot_name in names(plot_list)) {
  print(plot_list[[plot_name]])
}



################################################################
#                   Linear Regression                          #
################################################################
#---------------------------------------------------------------
# Linear regression joining populations together
#---------------------------------------------------------------
# Modify dataframe
head(dfm)
# Split dfm by chromosome
# Create an empty list to store the split dataframes
split_dataframes2 <- list()
# Use a for loop to split the dataframe by 'Group' and store in the list
for (group_value in chr) {
  subset_df2 <- subset(dfm, chromosome == group_value)
  split_dataframes2[[group_value]] <- subset_df2
}


# Calculate linear regression
# https://www.statology.org/bootstrapping-in-r/
# https://www.utstat.toronto.edu/~brunner/oldclass/appliedf12/lectures/2101f12BootstrapR.pdf
# https://rpubs.com/emmalaughlin/1002583
# https://towardsdatascience.com/bootstrap-regression-in-r-98bfe4ff5007
# https://stats.stackexchange.com/questions/316483/manually-bootstrapping-linear-regression-in-r

lr_results <- list()

for (group_value in chr) {
  subset_lr <- split_dataframes2[[group_value]]
  lr_results[[group_value]] <- lm(dist_to_Cla~rho, subset_lr)
  cat("Linear Regression for DataFrame:", group_value, "\n")
  print(summary(lr_results[[group_value]]))
  }

# Plots model results
for (group_value in chr) {
  plot(lr_results[[group_value]]) #plotting the model
  }

# Plots with ggplot: scatterplot with fitted regression line
ggplot(split_dataframes2[["Chr1"]], aes(x = dist_to_Cla, y = rho)) +
  geom_point() +
  stat_smooth(method = "lm", se=TRUE)

ggplot(split_dataframes2[["Chr2"]], aes(x = dist_to_Cla, y = rho)) +
  geom_point() +
  stat_smooth(method = "lm", se=TRUE)

ggplot(split_dataframes2[["Chr3"]], aes(x = dist_to_Cla, y = rho)) +
  geom_point() +
  stat_smooth(method = "lm", se=TRUE)

ggplot(split_dataframes2[["Chr4"]], aes(x = dist_to_Cla, y = rho)) +
  geom_point() +
  stat_smooth(method = "lm", se=TRUE)

# Bootstrapping
# Containers for the coefficients
sample_coef_intercept <- list()
sample_coef_x1 <- list()
model_bootstrap <- list()

for (group_value in chr) {
  for (i in 1:100) {
    subset_lr <- split_dataframes2[[group_value]]
    
    #Creating a resampled dataset from the sample data
    sample_d = subset_lr[sample(1:nrow(subset_lr), 
                                nrow(subset_lr), replace = TRUE), ]
    
    #Running the regression on these data
    model_bootstrap[[group_value]] <- lm(dist_to_Cla~rho, data = sample_d)
    
    #Saving the coefficients
    sample_coef_intercept[[group_value]] <-
      c(sample_coef_intercept, model_bootstrap$coefficients[1])
    
    sample_coef_x1[[group_value]] <-
      c(sample_coef_x1, model_bootstrap$coefficients[2])
  }
}

plot(model_bootstrap[["Chr1"]])
plot(model_bootstrap[["Chr2"]])
plot(model_bootstrap[["Chr3"]])
plot(model_bootstrap[["Chr4"]])

# Bootrapping worked but plots look a bit weird!
# Check Condor-paper, maybe that's why they defined binned mean values


#---------------------------------------------------
# Regression done according to Condor-paper
#---------------------------------------------------
max_dist <- 20e6
window_size <- 10e3

# Generate bins according to max-dist and window size
bin_size <- 10e3
bin_starts <- seq(0, max_dist-bin_size, by=10e3)
bin_ends <- bin_starts+bin_size

# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- list()

for (group_value in chr) {
  
  subset_df2 <- subset(dfm, chromosome == group_value)
  split_dataframes2[[group_value]] <- subset_df2
  
  for (i in 1:length(bin_starts)){
  # Get windows within distance bin
  tt=subset_df2[which(subset_df2$dist_to_Cla >= bin_starts[i] & 
                        subset_df2$dist_to_Cla < bin_ends[i]),]
  # Get mean rho (weighted)
  if (dim(tt)[1] > 0){
    n_sites <- tt$end-tt$start
    total_rho <- tt$rho*n_sites
    mean_rho_per_dist_bin[[group_value]][i] <- sum(total_rho)/sum(n_sites)
  } 
  else {
    mean_rho_per_dist_bin[[group_value]][i] <- NA
    }
  }
}


# Boootstrap (resample rows from df with replacement)
n_boot <- 10
boot_mat <- list()
mean_rho_per_dist_bin_boot <- list()

for (group_value in chr) {
  
  subset_df2 <- split_dataframes2[[group_value]] 
  boot_mat[[group_value]] <- matrix(nrow=n_boot, ncol=length(bin_starts))
  
  for (b in 1:n_boot){
  print(b)
  newrows <- sample(1:dim(subset_df2)[1], size=dim(subset_df2)[1], replace=T)
  df_boot <- subset_df2[newrows,]
  
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df_boot[which(subset_df2$dist_to_Cla >= bin_starts[i] & 
                       subset_df2$dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites <- tt$end-tt$start
      total_rho <- tt$rho*n_sites
      mean_rho_per_dist_bin_boot[[group_value]][i] <- sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin_boot[[group_value]][i] <- NA
    }
  }
  boot_mat[[group_value]][b,] <- mean_rho_per_dist_bin_boot[[group_value]]
  }	
}

scale=1e3
myylab=expression(paste("Mean ", rho, "/bp x 10"^-3, sep=""))

# Set plot y-axis limits by getting span of mean rho values and adding 5% padding
yrange=max(rbind(mean_rho_per_dist_bin[["Chr1"]], boot_mat[["Chr1"]]), 
           na.rm = T) - min(rbind(mean_rho_per_dist_bin[["Chr1"]], boot_mat[["Chr1"]]), 
                            na.rm = T)
ymax=max(rbind(mean_rho_per_dist_bin[["Chr1"]], boot_mat[["Chr1"]]), 
         na.rm = T) + .05*yrange
ymin=min(rbind(mean_rho_per_dist_bin[["Chr1"]], boot_mat[["Chr1"]]), 
         na.rm = T) - .05*yrange

# Initialize plot
par(mar=c(5,4,2,0.5))
plot(1:length(mean_rho_per_dist_bin[["Chr1"]]), 
     scale*mean_rho_per_dist_bin[["Chr1"]], 
     type="n", 
     ylim=scale*c(ymin, ymax), xlab="", ylab="")
title(xlab="Distance to nearest\nCla element (kb)", line=3.5)
title(ylab=myylab, line=2.5)

# Add ribbon showing range of bootstrap min and max rho values per bin
ribbon_col="#c8c8c8"
ribbon_xs=c(1:dim(boot_mat)[2], dim(boot_mat)[2]:1)
ribbon_ys=c(apply(boot_mat, 2, function(x) max(x, na.rm = TRUE)), rev(apply(boot_mat, 2, function(x) min(x, na.rm = TRUE))))
polygon(ribbon_xs, scale*ribbon_ys, col=ribbon_col, border=NA)

# Add the line for the actual mean rho per distance bin
lines(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, lwd=2)

###############################################################################
#Test with only Chr1
###############################################################################
max_dist <- 1e6
window_size <- 10e3

# Generate bins according to max-dist and window size
bin_size <- 10e3
bin_starts <- seq(0, max_dist-bin_size, by=10e3)
bin_ends <- bin_starts+bin_size

# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr1")) {
  
  subset_df2 <- subset(dfm, chromosome == group_value)
  subset_df2 <- subset_df2 %>%
    filter(Pop_Chr == "MF_Chr1")
  
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=subset_df2[which(subset_df2$dist_to_Cla >= bin_starts[i] & 
                          subset_df2$dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites <- tt$end-tt$start
      total_rho <- tt$rho*n_sites
      mean_rho_per_dist_bin[i] <- sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin[i] <- NA
    }
  }
}


# Boootstrap (resample rows from df with replacement)
n_boot <- 100
boot_mat <- matrix(nrow=n_boot, ncol=length(bin_starts))
mean_rho_per_dist_bin_boot <- NULL


for (b in 1:n_boot){
    print(b)
    newrows <- sample(1:dim(subset_df2)[1], size=dim(subset_df2)[1], replace=T)
    df_boot <- subset_df2[newrows,]
    
    for (i in 1:length(bin_starts)){
      # Get windows within distance bin
      tt=df_boot[which(subset_df2$dist_to_Cla >= bin_starts[i] & 
                         subset_df2$dist_to_Cla < bin_ends[i]),]
      # Get mean rho (weighted)
      if (dim(tt)[1] > 0){
        n_sites <- tt$end-tt$start
        total_rho <- tt$rho*n_sites
        mean_rho_per_dist_bin_boot[i] <- sum(total_rho)/sum(n_sites)
      } 
      else {
        mean_rho_per_dist_bin_boot[i] <- NA
      }
    }
    boot_mat[b,] <- mean_rho_per_dist_bin_boot
}	



scale=1e3
myylab=expression(paste("Mean ", rho, "/bp x 10"^-4, sep=""))

# Set plot y-axis limits by getting span of mean rho values and adding 5% padding
yrange=max(rbind(mean_rho_per_dist_bin, boot_mat), 
           na.rm = T) - min(rbind(mean_rho_per_dist_bin, boot_mat), 
                            na.rm = T)
ymax=max(rbind(mean_rho_per_dist_bin, boot_mat), 
         na.rm = T) + .05*yrange
ymin=min(rbind(mean_rho_per_dist_bin, boot_mat), 
         na.rm = T) - .05*yrange

# Initialize plot
par(mar=c(5,4,2,0.5))
plot(1:length(mean_rho_per_dist_bin), 
     mean_rho_per_dist_bin, 
     type="n", 
     ylim=c(ymin, ymax), xlab="", ylab="")
title(xlab="Distance to nearest\nCla element (kb)", line=3.5)
title(ylab=myylab, line=2.5)

# Add ribbon showing range of bootstrap min and max rho values per bin
ribbon_col <- "#c8c8c8"
ribbon_xs <- c(1:dim(boot_mat)[2], dim(boot_mat)[2]:1)
ribbon_ys <- c(apply(boot_mat, 2, function(x) max(x, na.rm = TRUE)), 
            rev(apply(boot_mat, 2, function(x) min(x, na.rm = TRUE))))
polygon(ribbon_xs, ribbon_ys, col=ribbon_col, border=NA)

# Add the line for the actual mean rho per distance bin
lines(1:length(mean_rho_per_dist_bin), mean_rho_per_dist_bin, lwd=2)


# Move Everything Of Regression To Separate Script To Test With Mean


