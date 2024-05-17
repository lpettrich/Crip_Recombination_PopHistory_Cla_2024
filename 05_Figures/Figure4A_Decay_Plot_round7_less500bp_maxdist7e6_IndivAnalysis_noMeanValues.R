###########################################################
#                                                         #
# Statistics recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           August 2023                   #
###########################################################

## Clean environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)

# Set directory
# setwd("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")
setwd("~/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")


getwd()

#######################################################################################
# READ IN AND MODIFY DATA
#######################################################################################
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
# common_path = "C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv/02_less500bp"
common_path = "~/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv/02_less500bp"

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
  group_by(chromosome, start, end, Pop, Ind_Chr) %>%
  summarise(
    N = n(),
    mean_rho = (rho),
    mean_dist_to_Cla = (dist_to_Cla),
    sd.rho = sd(rho),
    sd.dist.Cla = sd(dist_to_Cla))
head(dfm_mean)

## dfm_mean is data frame used for plotting

#######################################################################################
# DATA ANALYSIS
#######################################################################################

# First look
dfm_mean %>% 
  #filter(Pop_Chr == "SS_Chr1") %>%
  filter(chromosome == "Chr1") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop)) + geom_point(size=.5) 


dfm_mean %>% 
  #filter(Pop_Chr == "SS_Chr1") %>%
  filter(chromosome == "Chr2") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop)) + geom_point(size=.5) 

dfm_mean %>% 
  #filter(Pop_Chr == "SS_Chr1") %>%
  filter(chromosome == "Chr3") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop)) + geom_point(size=.5) 

dfm_mean %>% 
  #filter(Pop_Chr == "SS_Chr1") %>%
  filter(chromosome == "Chr4") %>%
  ggplot(aes(x = mean_dist_to_Cla, y = mean_rho, colour = Pop)) + geom_point(size=.5) 

ggplot(dfm_mean, aes(x = mean_dist_to_Cla, y = mean_rho)) +
  geom_point() +
  facet_wrap(~chromosome, nrow = 2)

# Not really decisive
# Make it better!

#####################
#     Decay Plots   #
#####################
#-------------------------------------------------------------------
# According to https://github.com/jarobin/CaliforniaCondorGenome2021/blob/main/gc_cpg_recombinationrate.sh
#-------------------------------------------------------------------

#####################################################################################
#####################################################################################
#-------------------------------------------------------------------
# Remove incomplete windows
#-------------------------------------------------------------------
window_size <- 10e3
dfm_mean <- dfm_mean %>%
  mutate(wnd_size = end-start)

unique(dfm_mean$wnd_size) # window size should be only 10,000

dfm_mean <- dfm_mean %>% 
  filter(wnd_size == window_size)

unique(dfm_mean$wnd_size) # window size should be only 10,000
# it worked!

#-------------------------------------------------------------------
# Remove values that are too distant
#-------------------------------------------------------------------
max(dfm_mean$mean_dist_to_Cla)
min(dfm_mean$mean_dist_to_Cla)

summary(dfm_mean$mean_dist_to_Cla)
quantile(dfm_mean$mean_dist_to_Cla)
quantile(dfm_mean$mean_dist_to_Cla, prob=c(.25, .5,.75), type=1)

# How many entries would be removed with threshold?
dfm_mean %>%
  group_by(chromosome) %>%
  summarize(n_entries = n(),
            n_gt20e6 = sum(mean_dist_to_Cla >= 7e6),
            p_gt20e6  = n_gt20e6 / n_entries)

dfm_mean %>%
  group_by(chromosome) %>%
  summarize(mean_max_dist = max(mean_dist_to_Cla))
#####################################################################################
#####################################################################################

#-------------------------------------------------------------------
# Remove values that are too distant
#-------------------------------------------------------------------
# read-in chromosome and scaffold data as comparison
#chr.data <- read.table("/home/alle/recombination-map/MELT/crip4.0_genome_chr_length.csv", sep = ";", header = TRUE)
#chr.data



# we set it to 7,000,000
max_dist=7e6
dfm_mean <- dfm_mean %>% 
  filter(mean_dist_to_Cla <= max_dist)
summary(dfm_mean$mean_dist_to_Cla) 
boxplot(dfm_mean$mean_dist_to_Cla)


#---------------------------------------------------
# Regression done according to Condor-paper
#---------------------------------------------------

###############################################################################
#Test with several populations
###############################################################################
#max_dist <- 7e6
#window_size <- 10e3

# Generate bins according to max-dist and window size
# Generate bins 
bin_size=10e3
bin_starts=seq(0, max_dist-bin_size, by=10e3)
bin_ends=bin_starts+bin_size

# CHROMOSOME 1
# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr1")) {
  
  df <- subset(dfm_mean, chromosome == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df[which(df$mean_dist_to_Cla >= bin_starts[i] & df$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin[i]=NA
    }
  }
}


# Boootstrap (resample rows from df with replacement)
n_boot=100
boot_mat=matrix(nrow=n_boot, ncol=length(bin_starts))

for (b in 1:n_boot){
  print(b)
  newrows=sample(1:dim(df)[1], size=dim(df)[1], replace=T)
  df_boot=df[newrows,]
  mean_rho_per_dist_bin_boot=NULL
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df_boot[which(df_boot$mean_dist_to_Cla >= bin_starts[i] & df_boot$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin_boot[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin_boot[i]=NA
    }
  }
  boot_mat[b,]=mean_rho_per_dist_bin_boot
}	



scale=1e3
myylab=expression(paste("Mean ", rho, "/bp x 10"^-3, sep=""))

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
plot(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, type="n", ylim=scale*c(ymin, ymax), xlab="", ylab="")
title(xlab="Distance to nearest Cla-element\n(in 10 kb windows)", line=3.5)
title(ylab=myylab, line=2.5)

# Add ribbon showing range of bootstrap min and max rho values per bin
ribbon_col="#c8c8c8"
ribbon_xs=c(1:dim(boot_mat)[2], dim(boot_mat)[2]:1)
ribbon_ys=c(apply(boot_mat, 2, function(x) max(x, na.rm = TRUE)), rev(apply(boot_mat, 2, function(x) min(x, na.rm = TRUE))))
polygon(ribbon_xs, scale*ribbon_ys, col=ribbon_col, border=NA)

# Add the line for the actual mean rho per distance bin
lines(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, lwd=2)

# AS GGPLOT:
# Define theme
mytheme <- theme(axis.text.x = element_text(size = 14, color = "black"),
                 axis.text.y = element_text(size = 14, color = "black"),
                 axis.title.y = element_text(size = 14,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 14,face = "bold", color = "black"),
                 title = element_text(size = 12, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("lightgrey"),
                 panel.grid.minor = element_line("white"),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.background = element_rect(fill="white"),
                 legend.position="bottom")

# Combine mean_rho_per_dist_bin and boot_mat into a data frame
plot_data <- data.frame(
  x = 1:length(mean_rho_per_dist_bin),
  y = scale * mean_rho_per_dist_bin,
  ymin = scale * apply(boot_mat, 2, function(x) min(x, na.rm = TRUE)),
  ymax = scale * apply(boot_mat, 2, function(x) max(x, na.rm = TRUE))
)

# Plot using ggplot2
c1 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
  geom_line(linewidth = 1) +
  ylim(scale * c(ymin, ymax)) +
  #xlim(0, 7e6/10e3) +
  scale_x_continuous(n.breaks = 7, limits=c(0, 7e6/10e3)) +
  labs(
    x = "Distance to nearest Cla-element\n(in 10 kb windows)",
    y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
  ) +
  mytheme
c1

# CHROMOSOME 2
# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr2")) {
  
  df <- subset(dfm_mean, chromosome == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df[which(df$mean_dist_to_Cla >= bin_starts[i] & df$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin[i]=NA
    }
  }
}


# Boootstrap (resample rows from df with replacement)
n_boot=100
boot_mat=matrix(nrow=n_boot, ncol=length(bin_starts))

for (b in 1:n_boot){
  print(b)
  newrows=sample(1:dim(df)[1], size=dim(df)[1], replace=T)
  df_boot=df[newrows,]
  mean_rho_per_dist_bin_boot=NULL
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df_boot[which(df_boot$mean_dist_to_Cla >= bin_starts[i] & df_boot$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin_boot[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin_boot[i]=NA
    }
  }
  boot_mat[b,]=mean_rho_per_dist_bin_boot
}		



scale=1e3
myylab=expression(paste("Mean ", rho, "/bp x 10"^-3, sep=""))

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
plot(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, type="n", ylim=scale*c(ymin, ymax), xlab="", ylab="")
title(xlab="Distance to nearest Cla-element\n(in 10 kb windows)", line=3.5)
title(ylab=myylab, line=2.5)

# Add ribbon showing range of bootstrap min and max rho values per bin
ribbon_col="#c8c8c8"
ribbon_xs=c(1:dim(boot_mat)[2], dim(boot_mat)[2]:1)
ribbon_ys=c(apply(boot_mat, 2, function(x) max(x, na.rm = TRUE)), rev(apply(boot_mat, 2, function(x) min(x, na.rm = TRUE))))
polygon(ribbon_xs, scale*ribbon_ys, col=ribbon_col, border=NA)

# Add the line for the actual mean rho per distance bin
lines(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, lwd=2)

# AS GGPLOT (with help of ChatGPT):
# Combine mean_rho_per_dist_bin and boot_mat into a data frame
plot_data <- data.frame(
  x = 1:length(mean_rho_per_dist_bin),
  y = scale * mean_rho_per_dist_bin,
  ymin = scale * apply(boot_mat, 2, function(x) min(x, na.rm = TRUE)),
  ymax = scale * apply(boot_mat, 2, function(x) max(x, na.rm = TRUE))
)

# Plot using ggplot2
c2 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
  geom_line(lwd = 1) +
  ylim(scale * c(ymin, ymax)) +
  #xlim(0, 7e6/10e3) +
  scale_x_continuous(n.breaks = 7, limits=c(0, 7e6/10e3)) +
  labs(
    x = "Distance to nearest Cla-element\n(in 10 kb windows)",
    y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
  ) +
  mytheme

c2
# CHROMOSOME 3
# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr3")) {
  
  df <- subset(dfm_mean, chromosome == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df[which(df$mean_dist_to_Cla >= bin_starts[i] & df$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin[i]=NA
    }
  }
}


# Boootstrap (resample rows from df with replacement)
n_boot=100
boot_mat=matrix(nrow=n_boot, ncol=length(bin_starts))

for (b in 1:n_boot){
  print(b)
  newrows=sample(1:dim(df)[1], size=dim(df)[1], replace=T)
  df_boot=df[newrows,]
  mean_rho_per_dist_bin_boot=NULL
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df_boot[which(df_boot$mean_dist_to_Cla >= bin_starts[i] & df_boot$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin_boot[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin_boot[i]=NA
    }
  }
  boot_mat[b,]=mean_rho_per_dist_bin_boot
}	



scale=1e3
myylab=expression(paste("Mean ", rho, "/bp x 10"^-3, sep=""))

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
plot(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, type="n", ylim=scale*c(ymin, ymax), xlab="", ylab="")
title(xlab="Distance to nearest Cla-element\n(in 10 kb windows)", line=3.5)
title(ylab=myylab, line=2.5)

# Add ribbon showing range of bootstrap min and max rho values per bin
ribbon_col="#c8c8c8"
ribbon_xs=c(1:dim(boot_mat)[2], dim(boot_mat)[2]:1)
ribbon_ys=c(apply(boot_mat, 2, function(x) max(x, na.rm = TRUE)), rev(apply(boot_mat, 2, function(x) min(x, na.rm = TRUE))))
polygon(ribbon_xs, scale*ribbon_ys, col=ribbon_col, border=NA)

# Add the line for the actual mean rho per distance bin
lines(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, lwd=2)

# AS GGPLOT:
# Combine mean_rho_per_dist_bin and boot_mat into a data frame
plot_data <- data.frame(
  x = 1:length(mean_rho_per_dist_bin),
  y = scale * mean_rho_per_dist_bin,
  ymin = scale * apply(boot_mat, 2, function(x) min(x, na.rm = TRUE)),
  ymax = scale * apply(boot_mat, 2, function(x) max(x, na.rm = TRUE))
)

# Plot using ggplot2
c3 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
  geom_line(lwd = 1) +
  ylim(scale * c(ymin, ymax)) +
  #xlim(0, 7e6/10e3) +
  scale_x_continuous(n.breaks = 7, limits=c(0, 7e6/10e3)) +
  labs(
    x = "Distance to nearest Cla-element\n(in 10 kb windows)",
    y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
  ) +
  mytheme

c3

# CHROMOSOME 4
# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr4")) {
  
  df <- subset(dfm_mean, chromosome == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df[which(df$mean_dist_to_Cla >= bin_starts[i] & df$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin[i]=NA
    }
  }
}


# Boootstrap (resample rows from df with replacement)
n_boot=100
boot_mat=matrix(nrow=n_boot, ncol=length(bin_starts))

for (b in 1:n_boot){
  print(b)
  newrows=sample(1:dim(df)[1], size=dim(df)[1], replace=T)
  df_boot=df[newrows,]
  mean_rho_per_dist_bin_boot=NULL
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df_boot[which(df_boot$mean_dist_to_Cla >= bin_starts[i] & df_boot$mean_dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$mean_rho*n_sites
      mean_rho_per_dist_bin_boot[i]=sum(total_rho)/sum(n_sites)
    } 
    else {
      mean_rho_per_dist_bin_boot[i]=NA
    }
  }
  boot_mat[b,]=mean_rho_per_dist_bin_boot
}	



scale=1e3
myylab=expression(paste("Mean ", rho, "/bp x 10"^-3, sep=""))

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
plot(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, type="n", ylim=scale*c(ymin, ymax), xlab="", ylab="")
title(xlab="Distance to nearest Cla-element\n(in 10 kb windows)", line=3.5)
title(ylab=myylab, line=2.5)

# Add ribbon showing range of bootstrap min and max rho values per bin
ribbon_col="#c8c8c8"
ribbon_xs=c(1:dim(boot_mat)[2], dim(boot_mat)[2]:1)
ribbon_ys=c(apply(boot_mat, 2, function(x) max(x, na.rm = TRUE)), rev(apply(boot_mat, 2, function(x) min(x, na.rm = TRUE))))
polygon(ribbon_xs, scale*ribbon_ys, col=ribbon_col, border=NA)

# Add the line for the actual mean rho per distance bin
lines(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, lwd=2)

# AS GGPLOT:
# Combine mean_rho_per_dist_bin and boot_mat into a data frame
plot_data <- data.frame(
  x = 1:length(mean_rho_per_dist_bin),
  y = scale * mean_rho_per_dist_bin,
  ymin = scale * apply(boot_mat, 2, function(x) min(x, na.rm = TRUE)),
  ymax = scale * apply(boot_mat, 2, function(x) max(x, na.rm = TRUE))
)

# Plot using ggplot2
c4 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
  geom_line(lwd = 1) +
  ylim(scale * c(ymin, ymax)) +
  #xlim(0, 7e6/10e3) +
  scale_x_continuous(n.breaks = 7, limits=c(0, 7e6/10e3)) +
  labs(
    x = "Distance to nearest Cla-element\n(in 10 kb windows)",
    y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
  ) +
  mytheme

c4
# ARRANGE ALL PLOTS TOGETHER
library(cowplot)

plot_grid(
  c1,
  c2,
  c3,
  c4,
  nrow = 2,
  ncol = 2,
  align = 'hv', axis = 'lrb'
)

#-
library(patchwork)

theme_blankx <- theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      #title=element_blank()
                      )

theme_blanky <- theme(axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      #title=element_blank()
                      )
                      

theme_blankxy <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       #title=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())

c1_blankx <- c1 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) + theme_blankx + labs(title = "Chromosome 1")

c2_blankxy <- c2 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) + theme_blankxy + labs(title = "Chromosome 2")

c3_noblank <- c3 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) + labs(title = "Chromosome 3")

c4_blanky <- c4 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) + theme_blanky + labs(title = "Chromosome 4")


p_merged <- c1_blankx + c2_blankxy + c3_noblank + c4_blanky +
  patchwork::plot_layout(ncol = 2, heights = c(1,1))

p_merged

save_plot("../../../Plots/DecayPlots_round7_less500bp.png", p_merged, nrow = 2, ncol = 2, base_asp = 1.9, dpi = 500, bg = "white")
save_plot("../../../Plots/DecayPlots_round7_less500bp.svg", p_merged, nrow = 2, ncol = 2, base_asp = 1.9, dpi = 500, bg = "white")


# How many entires are left in file?

dfm_mean %>%
  group_by(chromosome) %>%
  summarize(n_entries = n(),
            n_gt7e6 = sum(mean_dist_to_Cla >= 7e6),
            p_gt7e6  = n_gt7e6 / n_entries)
# 


# save chr1
c1_save <- c1 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) + labs(title = "Chromosome 1")

save_plot("../../../Plots/DecayPlots_round7_less500bp_max7e6_CHR1.png", c1_save, scale = 1.5, nrow = 1, ncol = 1, base_asp = 1.6, dpi = 500, bg = "white")
save_plot("../../../Plots/DecayPlots_round7_less500bp_max7e6_CHR1.svg", c1_save, scale = 1.5, nrow = 1, ncol = 1, base_asp = 1.6, bg = "white")


