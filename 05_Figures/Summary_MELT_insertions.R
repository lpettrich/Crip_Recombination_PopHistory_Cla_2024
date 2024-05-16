## Clean environment
rm(list = ls())

#3 Load libraries
library(tidyverse)
library(dplyr)
#library(chromoMap)
#library(scales)
library(ggplot2)

# Set directory
setwd("C:/Users/laupe/Documents/Uni-Köln/PhD/CRIP/RecombinationMap/MELT_AllPop/1-IndivAnalysis/")
getwd()

#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
# Here you define vectors for every data table you have
# This info will be added a column later

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
common_path = getwd()

# Read in bedtools closest results
# These dataframes will be used for analysis
# You add each population seperately
## MF
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "\\bMF\\d+\\.bwamem\\.sambl\\.sort\\.Cla1\\.tmp\\.bed$", # regex pattern, some explanation below
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
  pattern = "MG\\d+\\.bwamem\\.sambl\\.sort\\.Cla1\\.tmp\\.bed$", # regex pattern, some explanation below
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
  pattern = "NMF\\d+\\.bwamem\\.sambl\\.sort\\.Cla1\\.tmp\\.bed$", # regex pattern, some explanation below
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
  pattern = "SI\\d+\\.bwamem\\.sambl\\.sort\\.Cla1\\.tmp\\.bed$", # regex pattern, some explanation below
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
  pattern = "SS\\d+\\.bwamem\\.sambl\\.sort\\.Cla1\\.tmp\\.bed$", # regex pattern, some explanation below
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
names(data_lst_1) <- ind_prefix[1:4]
names(data_lst_2) <- ind_prefix[5:8]
names(data_lst_3) <- ind_prefix[9:12]
names(data_lst_4) <- ind_prefix[13:16]
names(data_lst_5) <- ind_prefix[17:20]

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

# Rename columns in data frame
# https://stackoverflow.com/questions/28648513/changing-column-names-in-a-list-of-data-frames-in-r
new_colnames <- c("chromosome", "start", "end", "insert", 
                  "score", "strand", "ind", "pop")

dfm <- dfm %>%
  rename_with(~ new_colnames)
head(dfm)

# Find unique and shared insertions
# Find inserts occurring in only one population
inserts_one_pop <- dfm %>% 
  group_by(insert) %>% 
  filter(n_distinct(pop) == 1)

# Find inserts occurring in all populations
inserts_all_pop <- dfm %>% 
  group_by(insert) %>% 
  filter(n_distinct(pop) == n_distinct(dfm$pop))

# Output
print("Inserts occurring in only one population:")
print(inserts_one_pop)

print("Inserts occurring in all populations:")
print(inserts_all_pop)

# Count the number of populations for each insert
pop_counts <- dfm %>% 
  group_by(insert) %>% 
  summarise(num_populations = n_distinct(pop))

# Output
print("Count of populations for each insert:")
print(pop_counts)

allpop <- pop_counts %>%
  filter(num_populations == 5) 
