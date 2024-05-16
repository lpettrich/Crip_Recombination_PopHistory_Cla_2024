###########################################################
#                                                         #
#           Histogram Cla-element cluster size            #
#                   Script by Laura Pettrich              #
#                           January 2024                  #
###########################################################

## Clean environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
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

common_path = "/home/alle/recombination-map/MELT/MELT-run2/Cla1/bed-files"

# Read in bed files
## MF
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "^MF_Chr.*_clustersize\\.bed$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_1

data_lst_1 = lapply(files_to_read_1, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## MG
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "^MG_Chr.*_clustersize\\.bed$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_2

data_lst_2 = lapply(files_to_read_2, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## NMF
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "^NMF_Chr.*_clustersize\\.bed$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_3

data_lst_3 = lapply(files_to_read_3, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## SI
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "^SI_Chr.*_clustersize\\.bed$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_4

data_lst_4 = lapply(files_to_read_4, read.table, header = FALSE)  # read all the matching files

# Read in bedtools closest results
## SS
files_to_read_5 = list.files(
  path = common_path,        # directory to search within
  pattern = "^SS_Chr.*_clustersize\\.bed$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_5

data_lst_5 = lapply(files_to_read_5, read.table, header = FALSE)  # read all the matching files

# V6 is autocorrected to TRUE instead of keeping the "T" (= Thymin) in data_lst_2 - 5
# Define data types and read in the data again
# Assuming data_lst_1 is a list of data frames
# data_lst_1 has the correct data types
data_types <- sapply(data_lst_1, function(df) sapply(df, class))
data_types

data_lst_2 = lapply(files_to_read_2, read.table, header = FALSE, colClasses = data_types[1:14])
data_lst_3 = lapply(files_to_read_3, read.table, header = FALSE, colClasses = data_types[1:14])
data_lst_4 = lapply(files_to_read_4, read.table, header = FALSE, colClasses = data_types[1:14])
data_lst_5 = lapply(files_to_read_5, read.table, header = FALSE, colClasses = data_types[1:14])
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

#-----------------------------------------------------------------

# Bind data frames together
dfm <- rbind(df1, df2, df3, df4, df5)
head(dfm)
str(dfm)

# We just need V1:3 and pop
dfm <- dfm %>%
  select(V1:3, pop)

# Rename columns in data frame
# https://stackoverflow.com/questions/28648513/changing-column-names-in-a-list-of-data-frames-in-r
new_colnames <- c("chromosome", "start", "end", "Pop_Chr")

dfm <- dfm %>%
  rename_with(~ new_colnames)

# Add coulmn with cluster size
dfm <- dfm %>%
  mutate(size = end - start)
head(dfm)

# Size is only 1 because it is a vcf file converted to bed
ggplot(dfm, aes(x = start, y = end)) + 
  geom_point() +
  facet_wrap(vars(Pop_Chr), ncol = 4)




# How large are the insertions?
# define theme
mytheme <- theme(axis.text.x = element_text(size = 14, color = "black"),
                 #axis.text.y = element_text(size = 16, color = "black"),
                 #axis.title.y = element_text(size = 18,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 18,face = "bold", color = "black"),
                 title = element_text(size = 17, color = "black"),
                 text = element_text(size=17, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("white"),
                 panel.grid.minor = element_line("white"),
                 panel.grid.major.x = element_line("lightgrey", linetype = 3),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.background = element_rect(fill="white"),
                 legend.position="bottom", 
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank()) 

h <- ggplot(dfm, aes(x = size)) +  
  geom_histogram(binwidth = 10, color="black", fill="lightgrey") +
  facet_wrap(vars(chromosome), ncol = 2) +
  #geom_vline(aes(xintercept=mean(size)),
  #             color="darkred", linetype="dashed", size=1) +
  mytheme +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 18,face = "bold", color = "black"))

h

ggplot(dfm, aes(x = size)) + 
  geom_histogram(aes(y=after_stat(density)), colour="black", fill="lightgrey", binwidth = 10) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_wrap(vars(chromosome), ncol = 2) +
  mytheme +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 18,face = "bold", color = "black"))

# Summary
head(dfm)

summary(dfm$size)

