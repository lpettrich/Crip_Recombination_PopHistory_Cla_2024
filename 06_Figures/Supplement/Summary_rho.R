# Summary stats on rho
#----------------------------

rm(list = ls())
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

# Load libraries
library(tidyverse)
library(scales)
library(cowplot)

# Set directory
setwd("/home/alle/recombination-map")
getwd()


#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
ind <- c("MF1", "MF2", "MF3", "MF4", "MG2", "MG3", "MG4", "MG5", "NMF1", "NMF2", "NMF3", "NMF4", 
         "SI1", "SI2", "SI3", "SI4", "SS1", "SS2", "SS3", "SS4")
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

common_path = "/home/alle/recombination-map/ismc-output/ismc_newwnd"

# Chr1
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr1_ismc.rho.10kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_1 = lapply(files_to_read_1, read.table, header = TRUE)  # read all the matching files

# Chr2
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr2_ismc.rho.10kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_2 = lapply(files_to_read_2, read.table, header = TRUE)  # read all the matching files

# Chr3
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr3_ismc.rho.10kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
data_lst_3 = lapply(files_to_read_3, read.table, header = TRUE)  # read all the matching files

# Chr4
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "*Chr4_ismc.rho.10kb.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE,          # return the full path
)
data_lst_4 = lapply(files_to_read_4, read.table, header = TRUE)  # read all the matching files

#-------------------------------------------------------------------
# Rename columns
#-------------------------------------------------------------------
# Rename data list names
names(data_lst_1) <- ind
names(data_lst_2) <- ind
names(data_lst_3) <- ind
names(data_lst_4) <- ind

# Check if renaming worked fine
#Make empty list, otherwise loop won't work
check_rename <- list()

for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_1[i]), colnames(data_lst_1[[i]][4]))
}
check_rename


for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_2[i]), colnames(data_lst_2[[i]][4]))
}
check_rename

for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_3[i]), colnames(data_lst_3[[i]][4]))
}
check_rename

for (i in 1:20){
  check_rename[i] <- identical(names(data_lst_4[i]), colnames(data_lst_4[[i]][4]))
}
check_rename

# Rename columns in every data frame
# https://stackoverflow.com/questions/28648513/changing-column-names-in-a-list-of-data-frames-in-r
new_colnames <- c("chromosome", "start", "end", "rho")

data_lst_1 <- lapply(data_lst_1, rename_with, ~ new_colnames)
data_lst_2 <- lapply(data_lst_2, rename_with, ~ new_colnames)
data_lst_3 <- lapply(data_lst_3, rename_with, ~ new_colnames)
data_lst_4 <- lapply(data_lst_4, rename_with, ~ new_colnames)

# Change names in chromosome column to the correct chromosome
for (i in 1:20){
  data_lst_1[[i]][1] <- "Chr1"
}

for (i in 1:20){
  data_lst_2[[i]][1] <- "Chr2"
}

for (i in 1:20){
  data_lst_3[[i]][1] <- "Chr3"
}

for (i in 1:20){
  data_lst_4[[i]][1] <- "Chr4"
}

# Add column with individual name
for (i in 1:20){
  data_lst_1[[i]]$individual <- names(data_lst_1[i])
}

for (i in 1:20){
  data_lst_2[[i]]$individual <- names(data_lst_2[i])
}

for (i in 1:20){
  data_lst_3[[i]]$individual <- names(data_lst_3[i])
}

for (i in 1:20){
  data_lst_4[[i]]$individual <- names(data_lst_4[i])
}

# Merge dataframes into one big dataframe (df1 = chr1, df2 = chr2, ...)
df1 <- bind_rows(data_lst_1)
df2 <- bind_rows(data_lst_2)
df3 <- bind_rows(data_lst_3)
df4 <- bind_rows(data_lst_4)

df <-  rbind(df1,df2,df3,df4)

unique(df1$chromosome)
unique(df2$chromosome)
unique(df3$chromosome)
unique(df4$chromosome)

# Remove incomplete windows
df <- df %>%
  mutate(wnd_size = end-start)

unique(df$wnd_size) # window size should be only 10,000

window_size = 10e3

df <- df %>% 
  filter(wnd_size == window_size)

unique(df$wnd_size)
unique(df$chromosome)

# also for the separate df
df1 <- df1 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

df2 <- df2 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

df3 <- df3 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

df4 <- df4 %>%
  mutate(wnd_size = end-start) %>% 
  filter(wnd_size == window_size)

unique(df4$wnd_size)

#----------------------------------------------------------
# Calculate mean and confidence interval
## Chr1
head(df1)
mf_mean1 <- df1 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))

head(mf_mean1)

mg_mean1 <- df1 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(mg_mean1)

nmf_mean1 <- df1 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(nmf_mean1)

si_mean1 <- df1 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(si_mean1)

ss_mean1 <- df1 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(ss_mean1)

## Chr2
head(df2)
mf_mean2 <- df2 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(mf_mean2)

mg_mean2 <- df2 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(mg_mean2)

nmf_mean2 <- df2 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(nmf_mean2)

si_mean2 <- df2 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(si_mean2)

ss_mean2 <- df2 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(ss_mean2)

## Chr3
head(df3)
mf_mean3 <- df3 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(mf_mean3)

mg_mean3 <- df3 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(mg_mean3)

nmf_mean3 <- df3 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(nmf_mean3)

si_mean3 <- df3 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(si_mean3)

ss_mean3 <- df3 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(ss_mean3)

## Chr4
head(df4)
mf_mean4 <- df4 %>% 
  filter(individual == "MF1" | individual == "MF2" | 
           individual == "MF3" | individual == "MF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(mf_mean4)

mg_mean4 <- df4 %>% 
  filter(individual == "MG2" | individual == "MG3" | 
           individual == "MG4" | individual == "MG5") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(mg_mean4)

nmf_mean4 <- df4 %>% 
  filter(individual == "NMF1" | individual == "NMF2" | 
           individual == "NMF3" | individual == "NMF4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(nmf_mean4)

si_mean4 <- df4 %>% 
  filter(individual == "SI1" | individual == "SI2" | 
           individual == "SI3" | individual == "SI4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(si_mean4)

ss_mean4 <- df4 %>% 
  filter(individual == "SS1" | individual == "SS2" | 
           individual == "SS3" | individual == "SS4") %>%
  group_by(chromosome, start, end)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(ss_mean4)



#----------------------------------------------------------
# SUMMARIES
#----------------------------------------------------------
# MF
summary(mf_mean1)
summary(mf_mean2)
summary(mf_mean3)
summary(mf_mean4)

# MG
summary(mg_mean1)
summary(mg_mean2)
summary(mg_mean3)
summary(mg_mean4)

# NMF
summary(nmf_mean1)
summary(nmf_mean2)
summary(nmf_mean3)
summary(nmf_mean4)

# SI
summary(si_mean1)
summary(si_mean2)
summary(si_mean3)
summary(si_mean4)

# SS
summary(ss_mean1)
summary(ss_mean2)
summary(ss_mean3)
summary(ss_mean4)


#---------------------------------------------------------------------
# Global mean value of rho for each chromosome of each population
#---------------------------------------------------------------------

str(df)

df_mean <- df %>%
  mutate(population = str_remove(individual, "\\d+$"))

df_mean <- df_mean %>% 
  group_by(population, chromosome)  %>% 
  summarise(
    N = n(),
    mean_rho = mean(rho),
    median_rho = median(rho),
    Q25 = quantile(rho, 0.25),
    Q75 = quantile(rho, 0.75),
    SE.low = mean_rho - (sd(rho)/sqrt(N)),
    SE.high = mean_rho + (sd(rho)/sqrt(N))
  )
head(df_mean)

# Add total mean per pop
df_totalpop <- df %>%
  mutate(population = str_remove(individual, "\\d+$")) 

df_totalpop <- df_totalpop %>%
  group_by(population)  %>% 
  summarise(
    N = n(),
    mean_rho = mean(rho),
    median_rho = median(rho),
    Q25 = quantile(rho, 0.25),
    Q75 = quantile(rho, 0.75),
    SE.low = mean_rho - (sd(rho)/sqrt(N)),
    SE.high = mean_rho + (sd(rho)/sqrt(N))
  )
head(df_totalpop)


# Add total mean per chr
df_totalchr <- df %>%
  mutate(population = str_remove(individual, "\\d+$")) 

df_totalchr <- df_totalchr %>%
  group_by(chromosome)  %>% 
  summarise(
    N = n(),
    mean_rho = mean(rho),
    median_rho = median(rho),
    Q25 = quantile(rho, 0.25),
    Q75 = quantile(rho, 0.75),
    SE.low = mean_rho - (sd(rho)/sqrt(N)),
    SE.high = mean_rho + (sd(rho)/sqrt(N))
  )

# Total mean
df_total <- df %>%
  mutate(population = str_remove(individual, "\\d+$")) 

df_total <- df_total %>%
  summarise(
    N = n(),
    mean_rho = mean(rho),
    median_rho = median(rho),
    Q25 = quantile(rho, 0.25),
    Q75 = quantile(rho, 0.75),
    SE.low = mean_rho - (sd(rho)/sqrt(N)),
    SE.high = mean_rho + (sd(rho)/sqrt(N))
  )
head(df_total)


# Bind them together
df_totalchr$population <- "total"
df_totalpop$chromosome <- "total"
df_total$population <- "total"
df_total$chromosome <- "total"

df_mean <- rbind(df_mean, df_totalchr, df_totalpop, df_total)

# Save csv
write.csv(df_mean, "/home/alle/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/global_mean_rho_population-wise_chromosome-wise.csv", row.names = FALSE)


# Not included
# Group by population and calculate the Pearson correlation between start position and rho
df_stats <- df %>%
  mutate(population = str_remove(individual, "\\d+$")) 


# Function to compute correlation and p-value
get_correlation <- function(x, y) {
  test <- cor.test(x, y, method = "pearson")  # You can change to "spearman" if needed
  tibble(correlation = test$estimate,
         p_value = test$p.value)
}

# Apply to each population
correlations <- df_stats %>%
  group_by(population) %>%
  summarise(result = list(get_correlation(start, rho))) %>%
  unnest(cols = c(result))

# View the table with correlation and p-values
print(correlations)




