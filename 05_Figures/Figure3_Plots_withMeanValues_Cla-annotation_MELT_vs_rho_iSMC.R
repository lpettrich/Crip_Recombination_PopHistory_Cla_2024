
###########################################################
#                                                         #
# Comparison recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           August 2023                   #
###########################################################

## Clean environment
rm(list = ls())
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

# Load libraries
library(tidyverse)
library(dplyr)
library(scales)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)
library(gridExtra)

# Set directory
setwd("/home/alle/recombination-map")
getwd()

#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
ind <- c("MF1", "MF2", "MF3", "MF4", "MG2", "MG3", "MG4", "MG5", "NMF1", "NMF2", "NMF3", "NMF4", 
         "SI1", "SI2", "SI3", "SI4", "SS1", "SS2", "SS3", "SS4")
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

common_path = "/home/alle/recombination-map/ismc-output/ismc/"

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

# Remove incomplete windows
df <- df %>%
  mutate(wnd_size = end-start)

unique(df$wnd_size) # window size should be only 10,000

window_size = 10e3

df <- df %>% 
  filter(wnd_size == window_size)

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

pop <-  c("MF", 
          "NMF", 
          "MG", 
          "SI", 
          "SS")

#my_palette <- c("#355BA4",
#                "#8f4e8a",
#                "#287f97",
#                 "#b26712",
#                "#AB3232")

my_palette <- c("#4C83EB",
                "#CD70C6",
                "#3AB6D8",
                "#FF941A",
                "#AB3232")

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

# Chr1 with log10 transform
p1_mf <- mf_mean1 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[1],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[1]) +
  scale_x_continuous(limits=c(1,last(df1$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p1_mf

p1_nmf <- nmf_mean1 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[2],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[2]) +
  scale_x_continuous(limits=c(1,last(df1$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "NMF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p1_nmf

p1_mg <- mg_mean1 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[3],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[3]) +
  scale_x_continuous(limits=c(1,last(df1$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MG population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p1_mg

p1_si <- si_mean1 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[4],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[4]) +
  scale_x_continuous(limits=c(1,last(df1$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SI population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p1_si

p1_ss <- ss_mean1 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[5],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[5]) +
  scale_x_continuous(limits=c(1,last(df1$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SS population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p1_ss


# Chr2 with log10 transform
p2_mf <- mf_mean2 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[1],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[1]) +
  scale_x_continuous(limits=c(1,last(df2$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p2_mf

p2_nmf <- nmf_mean2 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[2],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[2]) +
  scale_x_continuous(limits=c(1,last(df2$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "NMF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p2_nmf

p2_mg <- mg_mean2 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[3],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[3]) +
  scale_x_continuous(limits=c(1,last(df2$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MG population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p2_mg

p2_si <- si_mean2 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[4],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[4]) +
  scale_x_continuous(limits=c(1,last(df2$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SI population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p2_si

p2_ss <- ss_mean2 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[5],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[5]) +
  scale_x_continuous(limits=c(1,last(df2$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SS population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p2_ss

# Chr3 with log10 transform
p3_mf <- mf_mean3 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[1],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[1]) +
  scale_x_continuous(limits=c(1,last(df3$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p3_mf

p3_nmf <- nmf_mean3 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[2],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[2]) +
  scale_x_continuous(limits=c(1,last(df3$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "NMF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p3_nmf

p3_mg <- mg_mean3 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[3],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[3]) +
  scale_x_continuous(limits=c(1,last(df3$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MG population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p3_mg

p3_si <- si_mean3 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[4],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[4]) +
  scale_x_continuous(limits=c(1,last(df3$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SI population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p3_si

p3_ss <- ss_mean3 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[5],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[5]) +
  scale_x_continuous(limits=c(1,last(df3$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SS population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p3_ss

# Chr4 with log10 transform
p4_mf <- mf_mean4 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[1],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[1]) +
  scale_x_continuous(limits=c(1,last(df4$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p4_mf

p4_nmf <- nmf_mean4 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[2],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[2]) +
  scale_x_continuous(limits=c(1,last(df4$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "NMF population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p4_nmf

p4_mg <- mg_mean4 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[3],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[3]) +
  scale_x_continuous(limits=c(1,last(df4$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "MG population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p4_mg

p4_si <- si_mean4 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[4],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[4]) +
  scale_x_continuous(limits=c(1,last(df4$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SI population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p4_si

p4_ss <- ss_mean4 %>%
  ggplot(aes(y = mean_rho, x = start)) +  
  geom_ribbon(aes(ymin = SE.low,
                  ymax = SE.high),
              fill = my_palette[5],
              alpha = 0.5) +
  geom_line(size=0.5, colour = my_palette[5]) +
  scale_x_continuous(limits=c(1,last(df4$end)), labels = comma, expand = c(0.125,0.125)) +
  scale_y_continuous(limits = c(min(df3$rho), max(df4$rho)), trans="log10", labels = comma) +
  labs(x = "Position (bp)", y = ~rho, title = "SS population") +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  annotation_logticks(sides="l")

p4_ss

#----------------------------------------------------------------------------------------
# CLA-ELEMENT
##############################################
##    2. Read in Cla-element location file  ##
##        Population-wise annotation        ##
##############################################

# Bed file of extracted Cla-element locations extracted from RepeatMasker gff and now converted to csv-format
# Bed-file is composed of these columns (from https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html):

#  RepeatMasker annotation field |	BED column index 	| BED field
#---------------------------------------------------------------------
# Query sequence |	1 |	chromosome
# Query start |	2 |	start
# Query end |	3 |	stop
# Repeat name |	4  |	id
# Smith-Waterman score |	5 |	score
# Strand |	6 |	strand
# Percentage, substitutions |	7 	 
# Percentage, deleted bases |	8 	 
# Percentage, inserted bases |	9 	 
# Bases in query, past match |	10 	 
# Repeat class |	11 	 
# Bases in complement of the repeat consensus sequence |	12 	 
# Match start |	13 	 
# Match end |	14 	 
# Unique ID  |	15 	 
# Higher-scoring match (optional) |	16


Cla_colnames = c("chromosome", "start", "end", "repeatname")

# Read in csv file
files_to_read_5 = list.files(
  path = "/home/alle/recombination-map/MELT/MELT-run2/Cla1/",        # directory to search within
  pattern = "*.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE,          # return the full path
)


files_to_read_5

data_lst_5 = lapply(files_to_read_5, read.table, header = FALSE, na.string = "NA", sep = "\t", col.names = Cla_colnames)  # read all the matching files

# Convert to data frame
claMF <- bind_rows(data_lst_5[[1]])

claMG <- bind_rows(data_lst_5[[2]])

claNMF <- bind_rows(data_lst_5[[3]])

claSI <- bind_rows(data_lst_5[[4]])

claSS <- bind_rows(data_lst_5[[5]])

# Get length and add population
claMF <- claMF %>% 
  mutate(length = end-start,
         pop = "MF")

claMG <- claMG %>% 
  mutate(length = end-start,
         pop = "MG")

claNMF <- claNMF %>% 
  mutate(length = end-start,
         pop = "NMF")

claSI <- claSI %>% 
  mutate(length = end-start,
         pop = "SI")

claSS <- claSS %>% 
  mutate(length = end-start,
         pop = "SS")




#############################################
##      5. PLOTTING OF RHO AND CLA         ##
##############################################
# CHROMOMAP NOT ABLE TO PLOT CATEGRORICAL AND NUMERICAL DATA
# USE GGPLOT INSTEAD
# https://stackoverflow.com/questions/64103097/best-approach-to-visualise-presence-absence-of-events-in-multiple-groups

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
                 legend.background = element_rect(fill="white"),
                 legend.position="bottom") 


# Chromosome data
# https://stackoverflow.com/questions/55594521/make-a-plot-using-one-column-as-start-and-another-as-stop-for-boxes
#read-in chromosome and scaffold data 
chr.data <- read.table("/home/alle/recombination-map/MELT/crip4.0_genome_chr_length.csv", sep = ";", header = TRUE)
class(chr.data)
head(chr.data)

chr.data <- chr.data %>% 
  rename("chromosome" = "Name",
         "start"= "Start",
         "end"= "End")
str(chr.data)
chr.data

# PLOT IT
target <- c("Chr1", "Chr2", "Chr3", "Chr4")

chromosomes <- chr.data %>%
  filter(chromosome %in% target) %>%
  select(chromosome, end) %>%
  mutate(start = 0) %>%
  mutate(end = end/1e6)

chromosomes

###################################### 
###################################### to find conflicting commands
#library(conflicted)

d <- data.frame(a = 1:10, b = 1:10)
select(d, a)
#######################################

# PLOT FOR MF
genes <- claMF %>%
  filter(chromosome %in% target) %>%
  rename(geneStart = start,
         geneEnd = end) %>%
  mutate(type = "Cla1",
         geneStart = geneStart/1e6,
         geneEnd = geneEnd/1e6)

genes

ggplot(chromosomes) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(genes, 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) +
  facet_grid(chromosome~.) 

# Save one plot per chromosome
claMFchr1 <- ggplot(subset(chromosomes, chromosome %in% "Chr1")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr1"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claMFchr2 <- ggplot(subset(chromosomes, chromosome %in% "Chr2")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr2"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claMFchr3 <- ggplot(subset(chromosomes, chromosome %in% "Chr3")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr3"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claMFchr4 <- ggplot(subset(chromosomes, chromosome %in% "Chr4")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr4"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 



# PLOT FOR MG
genes <- claMG %>%
  filter(chromosome %in% target) %>%
  rename(geneStart = start,
         geneEnd = end) %>%
  mutate(type = "Cla1",
         geneStart = geneStart/1e6,
         geneEnd = geneEnd/1e6)

genes

ggplot(chromosomes) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(genes, 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) +
  facet_grid(chromosome~.) 

# Save one plot per chromosome
claMGchr1 <- ggplot(subset(chromosomes, chromosome %in% "Chr1")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr1"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claMGchr2 <- ggplot(subset(chromosomes, chromosome %in% "Chr2")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr2"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claMGchr3 <- ggplot(subset(chromosomes, chromosome %in% "Chr3")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr3"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claMGchr4 <- ggplot(subset(chromosomes, chromosome %in% "Chr4")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr4"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

# PLOT FOR NMF
genes <- claNMF %>%
  filter(chromosome %in% target) %>%
  rename(geneStart = start,
         geneEnd = end) %>%
  mutate(type = "Cla1",
         geneStart = geneStart/1e6,
         geneEnd = geneEnd/1e6)

genes

ggplot(chromosomes) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(genes, 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) +
  facet_grid(chromosome~.) 

# Save one plot per chromosome
claNMFchr1 <- ggplot(subset(chromosomes, chromosome %in% "Chr1")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr1"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claNMFchr2 <- ggplot(subset(chromosomes, chromosome %in% "Chr2")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr2"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claNMFchr3 <- ggplot(subset(chromosomes, chromosome %in% "Chr3")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr3"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claNMFchr4 <- ggplot(subset(chromosomes, chromosome %in% "Chr4")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr4"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

# PLOT FOR SI
genes <- claSI %>%
  filter(chromosome %in% target) %>%
  rename(geneStart = start,
         geneEnd = end) %>%
  mutate(type = "Cla1",
         geneStart = geneStart/1e6,
         geneEnd = geneEnd/1e6)

genes

ggplot(chromosomes) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(genes, 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) +
  facet_grid(chromosome~.) 
# Save one plot per chromosome
claSIchr1 <- ggplot(subset(chromosomes, chromosome %in% "Chr1")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr1"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claSIchr2 <- ggplot(subset(chromosomes, chromosome %in% "Chr2")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr2"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claSIchr3 <- ggplot(subset(chromosomes, chromosome %in% "Chr3")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr3"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claSIchr4 <- ggplot(subset(chromosomes, chromosome %in% "Chr4")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr4"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

# PLOT FOR SS
genes <- claSS %>%
  filter(chromosome %in% target) %>%
  rename(geneStart = start,
         geneEnd = end) %>%
  mutate(type = "Cla1",
         geneStart = geneStart/1e6,
         geneEnd = geneEnd/1e6)

genes

ggplot(chromosomes) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(genes, 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) +
  facet_grid(chromosome~.) 
# Save one plot per chromosome
claSSchr1 <- ggplot(subset(chromosomes, chromosome %in% "Chr1")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr1"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claSSchr2 <- ggplot(subset(chromosomes, chromosome %in% "Chr2")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr2"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claSSchr3 <- ggplot(subset(chromosomes, chromosome %in% "Chr3")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr3"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 

claSSchr4 <- ggplot(subset(chromosomes, chromosome %in% "Chr4")) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(subset(genes, chromosome %in% "Chr4"), 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels = comma) 



#----------------------------------------------------------------------------------------
# PUT PLOTS ON SAME PAGE
theme_blankx <- theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      title=element_blank())

theme_blankxy <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       title=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())


##Chr1
p1_mg <- p1_mg + 
  labs(title = "Chromosome 1") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(size = 11, color = "black")) 

claMGchr1 <- claMGchr1  + theme_blankx

p1_nmf <- p1_nmf + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claNMFchr1 <- claNMFchr1 + theme_blankx 

p1_mf <- p1_mf + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claMFchr1 <- claMFchr1  + theme_blankx 

p1_si <- p1_si + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claSIchr1 <- claSIchr1 + theme_blankx 

p1_ss <- p1_ss + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claSSchr1  <- claSSchr1 + theme(axis.text.x = element_text(size = 11, color = "black"))


p_chr1 <- p1_mg + claMGchr1 + p1_nmf + claNMFchr1 + p1_mf + claMFchr1 +  p1_si + claSIchr1 + p1_ss + claSSchr1 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))

##Chr2
p2_mg <- p2_mg + 
  labs(title = "Chromosome 2") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

claMGchr2 <- claMGchr2  + theme_blankx

p2_nmf <- p2_nmf + theme_blankxy

claNMFchr2 <- claNMFchr2  + theme_blankx

p2_mf <- p2_mf + theme_blankxy

claMFchr2 <- claMFchr2  + theme_blankx

p2_si <- p2_si +theme_blankxy

claSIchr2 <- claSIchr2  + theme_blankx

p2_ss <- p2_ss + theme_blankxy

claSSchr2  <- claSSchr2 + theme(axis.text.x = element_text(size = 11, color = "black")) 


p_chr2 <- p2_mg + claMGchr2 + p2_nmf + claNMFchr2 +  p2_mf + claMFchr2 + p2_si + claSIchr2 + p2_ss + claSSchr2 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))

##Chr3
p3_mg <- p3_mg + 
  labs(title = "Chromosome 3") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

claMGchr3 <- claMGchr3  + theme_blankx

p3_nmf <- p3_nmf + theme_blankxy

claNMFchr3 <- claNMFchr3  + theme_blankx

p3_mf <- p3_mf + theme_blankxy

claMFchr3 <- claMFchr3  + theme_blankx

p3_si <- p3_si + theme_blankxy

claSIchr3 <- claSIchr3  + theme_blankx

p3_ss <- p3_ss + theme_blankxy

claSSchr3 <- claSSchr3 + theme(axis.text.x = element_text(size = 11, color = "black"))  


p_chr3 <- p3_mg + claMGchr3 + p3_nmf + claNMFchr3 + p3_mf + claMFchr3 + p3_si + claSIchr3 + p3_ss + claSSchr3 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))

##Chr4
p4_mg <- p4_mg + 
  labs(title = "Chromosome 4") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

claMGchr4 <- claMGchr4  + theme_blankx

p4_nmf <- p4_nmf + theme_blankxy 

claNMFchr4 <- claNMFchr4  + theme_blankx

p4_mf <- p4_mf + theme_blankxy 

claMFchr4 <- claMFchr4  + theme_blankx

p4_si <- p4_si + theme_blankxy 

claSIchr4 <- claSIchr4  + theme_blankx

p4_ss <- p4_ss + theme_blankxy 

claSSchr4  <- claSSchr4 + theme(axis.text.x = element_text(size = 11, color = "black")) 


p_chr4 <- p4_mg + claMGchr4 + p4_nmf + claNMFchr4 + p4_mf + claMFchr4 + p4_si + claSIchr4 + p4_ss + claSSchr4 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))
#+ plot_layout(guides = "collect")




# Wrap merged plots to one element
a <- wrap_elements(p_chr1)
b <- wrap_elements(p_chr2)
c <- wrap_elements(p_chr3)
d <- wrap_elements(p_chr4)

# Patchwork graph with shared y-axis

layout <- "
AAAAABBBBCCCCDDDD"



allplots <- a + b + c + d + plot_layout(design = layout) 

#?plot_layout

gt <- patchwork::patchworkGrob(allplots)
gridExtra::grid.arrange(gt, 
                        left = textGrob("\u03c1/10 kb", gp = gpar(fontsize = 15, fontype = "bold"),
                                        just = "centre", rot = 90), 
                        bottom = textGrob("Position (Mb)", gp = gpar(fontsize = 15, fontype = "bold"),
                                          vjust = -1, just = "centre"))


#-------------------------------------------------
# Different layout 
#----------------------------------------------------------------------------------------
# PUT PLOTS ON SAME PAGE
theme_blankx <- theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      title=element_blank())

theme_blankxy <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       title=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())


##Chr1
p1_mg <- p1_mg + 
  labs(title = "Chromosome 1") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(size = 11, color = "black")) 

claMGchr1 <- claMGchr1  + theme_blankx

p1_nmf <- p1_nmf + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claNMFchr1 <- claNMFchr1 + theme_blankx 

p1_mf <- p1_mf + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claMFchr1 <- claMFchr1  + theme_blankx 

p1_si <- p1_si + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claSIchr1 <- claSIchr1 + theme_blankx 

p1_ss <- p1_ss + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claSSchr1  <- claSSchr1 + theme(axis.text.x = element_text(size = 11, color = "black"))


p_chr1 <- p1_mg + claMGchr1 + p1_nmf + claNMFchr1 + p1_mf + claMFchr1 + p1_si + claSIchr1 + p1_ss + claSSchr1 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))

##Chr2
p2_mg <- p2_mg + 
  labs(title = "Chromosome 2") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

claMGchr2 <- claMGchr2  + theme_blankx

p2_nmf <- p2_nmf + theme_blankxy

claNMFchr2 <- claNMFchr2  + theme_blankx

p2_mf <- p2_mf + theme_blankxy

claMFchr2 <- claMFchr2  + theme_blankx

p2_si <- p2_si +theme_blankxy

claSIchr2 <- claSIchr2  + theme_blankx

p2_ss <- p2_ss + theme_blankxy

claSSchr2  <- claSSchr2 + theme(axis.text.x = element_text(size = 11, color = "black")) 


p_chr2 <- p2_mg + claMGchr2 + p2_nmf + claNMFchr2 + p2_mf + claMFchr2 + p2_si + claSIchr2 + p2_ss + claSSchr2 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))

##Chr3
p3_mg <- p3_mg + 
  labs(title = "Chromosome 3") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(size = 11, color = "black")) 

claMGchr3 <- claMGchr3  + theme_blankx

p3_nmf <- p3_nmf + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claNMFchr3 <- claNMFchr3 + theme_blankx 

p3_mf <- p3_mf + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claMFchr3 <- claMFchr3  + theme_blankx 

p3_si <- p3_si + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claSIchr3 <- claSIchr3 + theme_blankx 

p3_ss <- p3_ss + theme_blankx  + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 11, color = "black"))

claSSchr3  <- claSSchr3 + theme(axis.text.x = element_text(size = 11, color = "black"))


p_chr3 <- p3_mg + claMGchr3 + p3_nmf + claNMFchr3 + p3_mf + claMFchr3 + p3_si + claSIchr3 + p3_ss + claSSchr3 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))

##Chr4
p4_mg <- p4_mg + 
  labs(title = "Chromosome 4") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

claMGchr4 <- claMGchr4  + theme_blankx

p4_nmf <- p4_nmf + theme_blankxy 

claNMFchr4 <- claNMFchr4  + theme_blankx

p4_mf <- p4_mf + theme_blankxy 

claMFchr4 <- claMFchr4  + theme_blankx

p4_si <- p4_si + theme_blankxy 

claSIchr4 <- claSIchr4  + theme_blankx

p4_ss <- p4_ss + theme_blankxy 

claSSchr4  <- claSSchr4 + theme(axis.text.x = element_text(size = 11, color = "black")) 


p_chr4 <- p4_mg + claMGchr4 + p4_nmf + claNMFchr4 + p4_mf + claMFchr4 + p4_si + claSIchr4 + p4_ss + claSSchr4 + 
  patchwork::plot_layout(ncol = 1, heights = c(4, 1, 4, 1, 4, 1, 4, 1, 4, 1))
#+ plot_layout(guides = "collect")




# Wrap merged plots to one element
a <- wrap_elements(p_chr1)
b <- wrap_elements(p_chr2)
c <- wrap_elements(p_chr3)
d <- wrap_elements(p_chr4)

# Patchwork graph with shared y-axis

layout <- "
AAAAABBBBB
CCCCCDDDDD"



allplots <- a + b + c + d + plot_layout(design = layout) 
#+ plot_annotation(tag_levels = "A")

#?plot_layout

gt <- patchwork::patchworkGrob(allplots)
plot_final <- gridExtra::grid.arrange(gt, 
                                      left = textGrob("\u03c1", gp = gpar(fontsize = 15, fontface = "bold"),
                                                      just = "centre", rot = 90), 
                                      bottom = textGrob("Position (Mb)", gp = gpar(fontsize = 15, fontface = "bold"),
                                                        vjust = -1, just = "centre"))
#plot_final <- wrap_elements(plot_final)

# Get legend
legend_df <- as.data.frame(cbind(pop,my_palette,c(1,2,3,4,5)))
legend_df

legend_df$pop <- factor(legend_df$pop ,                 # Relevel group factor
                        levels = c("MG", "NMF", "MF", "SI", "SS"))
legend_df$my_palette <- c("#3AB6D8", "#CD70C6", "#4C83EB", "#FF941A", "#AB3232")
legend_df$names <- factor(c("Hesse (MG)", "Lorraine (NMF)", "Rhône-Alpes (MF)", "Piemont (SI)", "Andalusia (SS)"))
legend_df
legend_df$names <- factor(legend_df$names,
                          levels = c("Hesse (MG)", "Lorraine (NMF)", "Rhône-Alpes (MF)", "Piemont (SI)", "Andalusia (SS)"))

lp <- ggplot(legend_df, aes(y=legend_df$names, x=V3, color=legend_df$names)) +
  geom_line(size=0.5) +
  scale_color_manual(values=legend_df$my_palette) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, title= "Population")) +
  mytheme + 
  theme(legend.title = element_text(size = 14, #face = "bold", 
                                    color = "black"))

legend_b <- get_legend(lp + theme(legend.position="bottom"))


pf1 <- wrap_elements(plot_final) / legend_b + plot_layout(heights = c(10,.0001))


pf2 <- plot_grid(plot_final,
                 legend_b,
                 align = 'h', axis = 't',
                 ncol= 1,
                 rel_heights = c(10,.5),
                 rel_widths = c(10,.75)
)

pf2
#save_plot("rec-landscape_and_Cla-element.png", pf1, nrow = 8, ncol = 4, base_asp = 1.8, dpi = 500, bg = "white", scale = 0.6)
save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/png-files/rec-landscape_and_Cla-element_v3.png", pf2, nrow = 8, ncol = 4, base_asp = 1.9, dpi = 1000, bg = "white", scale = 0.6)
save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/svg-files/rec-landscape_and_Cla-element_v3.svg", pf2, nrow = 8, ncol = 4, base_asp = 1.9, dpi = 1000, bg = "white", scale = 0.6)

# Global mean value of rho for each chromosome of each population
str(df)

population <- sub("\\d+$", "", names(df$individual))

df_mean <- df %>%
  mutate(population = str_remove(individual, "\\d+$"))

df_mean <- df_mean %>% 
  group_by(population, chromosome)  %>% 
  summarise(N = n(),
            mean_rho = mean(rho),
            SE.low = mean_rho - (sd(rho)/sqrt(N)),
            SE.high = mean_rho + (sd(rho)/sqrt(N)))
head(df_mean)


write.csv(df_mean, "/home/alle/sciebo/RecombinationLandscape_CRIP/Scripts/global_mean_rho_population-wise_chromosome-wise.csv", row.names = FALSE)

