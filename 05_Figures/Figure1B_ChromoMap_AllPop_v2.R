## Clean environment
rm(list = ls())

#3 Load libraries
library(tidyverse)
library(dplyr)
#library(chromoMap)
library(scales)
library(ggplot2)

# Set directory
setwd("/home/alle/recombination-map")
getwd()

#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
ind <- c("MF1", "MF2", "MF3", "MF4", "MG2", "MG3", "MG4", "MG5", "NMF1", "NMF2", "NMF3", "NMF4", 
         "SI1", "SI2", "SI3", "SI4", "SS1", "SS2", "SS3", "SS4")
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

#########################################################
##    2. Read in Cla-element location file             ##
##        Total annotation with every population       ##
#########################################################

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
files_to_read = list.files(
  path = "/home/alle/recombination-map/MELT/MELT-run2/Cla1/AllPop/",        # directory to search within
  pattern = "*.bedgraph", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE,          # return the full path
)

files_to_read

data_lst = lapply(files_to_read, read.table, header = FALSE, na.string = "NA", sep = "\t", col.names = Cla_colnames)  # read all the matching files

# Convert to data frame
claAll <- bind_rows(data_lst[[1]])

# Get length 
claAll <- claAll %>% 
  mutate(length = end-start)

##############################################
##           4. Chromomap                   ##
##############################################
#read-in chromosome and scaffold data 
chr.data <- read.table("/home/alle/recombination-map/MELT/crip4.0_genome_chr_length.csv", sep = ";", header = TRUE)
class(chr.data)
head(chr.data)

chr.data <- chr.data %>% 
  rename("chromosome" = "Name",
         "chr_start"= "Start",
         "chr_end"= "End")
str(chr.data)
chr.data

# create chromomap with MF Cla annotation
claAll <- claAll %>%
  select(repeatname, chromosome, start, end)
str(claAll)
unique(claAll$chromosome)


# With ggplot
# https://stackoverflow.com/questions/64103097/best-approach-to-visualise-presence-absence-of-events-in-multiple-groups
# https://stackoverflow.com/questions/55594521/make-a-plot-using-one-column-as-start-and-another-as-stop-for-boxes

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

# define chromosomes of interest
target <- c("Chr1", "Chr2", "Chr3", "Chr4")

# modify chromosome data frame to filter and start at 0 and get end in Mb
chromosomes <- chr.data %>%
  filter(chromosome %in% target) %>%
  select(chromosome, chr_end) %>%
  mutate(chr_start = 0) %>%
  mutate(chr_end = chr_end/1e6)

chromosomes

# filter Cla element data frame and mutate to get position in Mb
genes <- claAll %>%
  filter(chromosome %in% target) %>%
  rename(geneStart = start,
         geneEnd = end) %>%
  mutate(type = "Cla1",
         geneStart = geneStart/1e6,
         geneEnd = geneEnd/1e6)

genes

# plot it with facet_grid by chromosome
p <- ggplot(chromosomes) +
  geom_rect(mapping=aes(xmin=chr_start, xmax=chr_end, ymin=0, ymax=1), 
            fill="lightblue") +
  geom_rect(genes, 
            mapping=aes(xmin=geneStart, xmax=geneEnd, 
                        ymin=0, ymax=1, color=type, fill = "type"), color = "darkred", fill = "darkred") +
  mytheme + 
  labs(x = "Position (Mb)") +
  scale_x_continuous(labels = comma, n.breaks = 6, position = "top") +
  facet_grid(chromosome~.) 

p + theme(axis.text.x = element_text(size = 14, color = "black", angle = -90, hjust = 1))

# How large are the insertions?
size_cla <- genes %>%
  mutate(geneEnd = geneEnd * 10^6,
         geneStart = geneStart * 10^6) %>%
  mutate(size = geneEnd - geneStart)
  
head(size_cla)

mean_size <- size_cla %>%
  group_by(chromosome) %>%
  reframe(m_size = mean(size))
mean_size

h <- ggplot(size_cla, aes(x = size)) +  
  geom_histogram(binwidth = 50, color="black", fill="lightgrey") +
  facet_wrap(vars(chromosome), ncol = 2) +
  #geom_vline(aes(xintercept=mean(size)),
  #             color="darkred", linetype="dashed", size=1) +
  mytheme +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 18,face = "bold", color = "black"))

h


#calaculate optimal bindwidth based on the Freedman-Diaconis rule
bw <- 2 * IQR(size_cla$size) / length(size_cla$size)^(1/3) #https://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram

h <- size_cla %>%
  ggplot(aes(x = size)) + 
  geom_histogram(aes(y=after_stat(density), x = size), colour="black", fill="lightgrey", binwidth = bw) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_wrap(vars(chromosome), ncol = 2) +
  labs(x = "Cluster length (bp)", y = "Density") +
  mytheme +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 18,face = "bold", color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 18,face = "bold", color = "black"),
        axis.ticks.x = element_line(color = "black"))

h

#breaks <- quantile(size_cla$size/1000,seq(0,1,by=0.05))
bw <- 2 * IQR(size_cla$size/1000) / length(size_cla$size/1000)^(1/3) #https://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram

# define theme
mytheme <- theme(axis.text.x = element_text(size = 14, color = "black"),
                 #axis.text.y = element_text(size = 16, color = "black"),
                 #axis.title.y = element_text(size = 18,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 16,face = "bold", color = "black"),
                 title = element_text(size = 16, color = "black"),
                 text = element_text(size=14, color = "black"),
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
h <- size_cla %>%
  mutate(size = size/1000) %>%
  ggplot(aes(x = size)) + 
  geom_histogram(aes(y=after_stat(density)), colour="black", fill="lightgrey", binwidth = bw) +
  geom_density(alpha=.2, fill="#FF6666") +
  #geom_vline(aes(xintercept = mean(size)), linetype = "dashed") +
  facet_wrap(vars(chromosome), ncol = 2) +
  labs(x = "Cluster length (kb)", y = "Density") +
  mytheme +
  theme(axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16,face = "bold", color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16,face = "bold", color = "black"),
        axis.ticks.x = element_line(color = "black"))
# explanation why density higher than 1
h

library(cowplot)
save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/svg-files/DensityPlot_v1.svg", h, ncol = 1.5, nrow = 2, base_asp = 2, dpi = 500, bg = "white")
save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/svg-files/DensityPlot_v1.png", h, ncol = 1.5, nrow = 2, base_asp = 2, dpi = 500, bg = "white")

# Get summary
summary(size_cla$size)
