###########################################################
# #
# Statistics recombination spots and Cla-element location #
#   Script by Laura Pettrich  #
#   REVISED   #
# Juli 2024   #
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
dfm <- dfm %>%
  mutate(chromosome_arm = paste(chromosome, chr_arms, sep = "_"))

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


#---------------------------------------------------
# Regression done according to Condor-paper
#---------------------------------------------------
#window_size <- 100e3


# CHROMOSOME 1
# ARM 1

# define max size
max_dist <- max(dfm$dist_to_Cla[dfm$chromosome_arm == "Chr1_NA"])

# Generate bins according to max-dist and window size
# Generate bins 
bin_size=100e3
bin_starts=seq(0, max_dist-bin_size, by=100e3)
bin_ends=bin_starts+bin_size


# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr1_NA")) {
  
  df <- subset(dfm, chromosome_arm == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
# Get windows within distance bin
tt=df[which(df$dist_to_Cla >= bin_starts[i] & df$dist_to_Cla < bin_ends[i]),]
# Get mean rho (weighted)
if (dim(tt)[1] > 0){
  n_sites=tt$end-tt$start
  total_rho=tt$rho*n_sites
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
tt=df_boot[which(df_boot$dist_to_Cla >= bin_starts[i] & df_boot$dist_to_Cla < bin_ends[i]),]
# Get mean rho (weighted)
if (dim(tt)[1] > 0){
  n_sites=tt$end-tt$start
  total_rho=tt$rho*n_sites
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
title(xlab="Distance to nearest Cla-element\n(in 100 kb windows)", line=3.5)
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
c1.1 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
  geom_line(linewidth = 1) +
  ylim(scale * c(ymin, ymax)) +
  #xlim(0, 7e6/10e3) +
  scale_x_continuous(n.breaks = 7, limits=c(0, max_dist/100e3)) +
  labs(
x = "Distance to nearest Cla-element\n(in 100 kb windows)",
y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
  ) +
  mytheme
c1.1















# CHROMOSOME 2
# ARM 1

# define max size
max_dist <- max(dfm$dist_to_Cla[dfm$chromosome_arm == "Chr2_NA"])

# Generate bins according to max-dist and window size
# Generate bins 
bin_size=100e3
bin_starts=seq(0, max_dist-bin_size, by=100e3)
bin_ends=bin_starts+bin_size


# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr2_NA")) {
  
  df <- subset(dfm, chromosome_arm == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
# Get windows within distance bin
tt=df[which(df$dist_to_Cla >= bin_starts[i] & df$dist_to_Cla < bin_ends[i]),]
# Get mean rho (weighted)
if (dim(tt)[1] > 0){
  n_sites=tt$end-tt$start
  total_rho=tt$rho*n_sites
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
tt=df_boot[which(df_boot$dist_to_Cla >= bin_starts[i] & df_boot$dist_to_Cla < bin_ends[i]),]
# Get mean rho (weighted)
if (dim(tt)[1] > 0){
  n_sites=tt$end-tt$start
  total_rho=tt$rho*n_sites
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
title(xlab="Distance to nearest Cla-element\n(in 100 kb windows)", line=3.5)
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
  c2.1 <- ggplot(plot_data, aes(x = x, y = y)) +
geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
geom_line(linewidth = 1) +
ylim(scale * c(ymin, ymax)) +
#xlim(0, 7e6/100e3) +
scale_x_continuous(n.breaks = 7, limits=c(0, max_dist/100e3)) +
labs(
  x = "Distance to nearest Cla-element\n(in 100 kb windows)",
  y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
) +
mytheme
  c2.1


  




# CHROMOSOME 3
# ARM 1

# define max size
max_dist <- max(dfm$dist_to_Cla[dfm$chromosome_arm == "Chr3_NA"])

# Generate bins according to max-dist and window size
# Generate bins 
bin_size=100e3
bin_starts=seq(0, max_dist-bin_size, by=100e3)
bin_ends=bin_starts+bin_size


# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr3_NA")) {
  
  df <- subset(dfm, chromosome_arm == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
    # Get windows within distance bin
    tt=df[which(df$dist_to_Cla >= bin_starts[i] & df$dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$rho*n_sites
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
    tt=df_boot[which(df_boot$dist_to_Cla >= bin_starts[i] & df_boot$dist_to_Cla < bin_ends[i]),]
    # Get mean rho (weighted)
    if (dim(tt)[1] > 0){
      n_sites=tt$end-tt$start
      total_rho=tt$rho*n_sites
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
title(xlab="Distance to nearest Cla-element\n(in 100 kb windows)", line=3.5)
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
  c3.1 <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
    geom_line(linewidth = 1) +
    ylim(scale * c(ymin, ymax)) +
    #xlim(0, 7e6/100e3) +
    scale_x_continuous(n.breaks = 7, limits=c(0, max_dist/100e3)) +
    labs(
      x = "Distance to nearest Cla-element\n(in 100 kb windows)",
      y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
    ) +
    mytheme
  c3.1
  








    
# CHROMOSOME 4
# define max size
max_dist <- max(dfm$dist_to_Cla[dfm$chromosome_arm == "Chr4_NA"])

# Generate bins according to max-dist and window size
# Generate bins 
bin_size=100e3
bin_starts=seq(0, max_dist-bin_size, by=100e3)
bin_ends=bin_starts+bin_size


# Calculate mean rho per distance bin
mean_rho_per_dist_bin <- NULL

for (group_value in c("Chr4_NA")) {
  
  df <- subset(dfm, chromosome_arm == group_value)
  df <- df
  
  for (i in 1:length(bin_starts)){
# Get windows within distance bin
tt=df[which(df$dist_to_Cla >= bin_starts[i] & df$dist_to_Cla < bin_ends[i]),]
# Get mean rho (weighted)
if (dim(tt)[1] > 0){
  n_sites=tt$end-tt$start
  total_rho=tt$rho*n_sites
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
tt=df_boot[which(df_boot$dist_to_Cla >= bin_starts[i] & df_boot$dist_to_Cla < bin_ends[i]),]
# Get mean rho (weighted)
if (dim(tt)[1] > 0){
  n_sites=tt$end-tt$start
  total_rho=tt$rho*n_sites
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
title(xlab="Distance to nearest Cla-element\n(in 100 kb windows)", line=3.5)
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
      c4 <- ggplot(plot_data, aes(x = x, y = y)) +
        geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#c8c8c8") +
        geom_line(linewidth = 1) +
        ylim(scale * c(ymin, ymax)) +
        #xlim(0, 7e6/100e3) +
        scale_x_continuous(n.breaks = 7, limits=c(0, max_dist/100e3)) +
        labs(
          x = "Distance to nearest Cla-element\n(in 100 kb windows)",
          y = expression(bold(paste("Mean ",  "\u03c1", "/bp x 10"^"-3", sep = "")))
        ) +
        mytheme
      c4

# ARRANGE ALL PLOTS TOGETHER
library(cowplot)

plot_grid(
  c1.1,
  c2.1,
  c3.1,
  #c1.2,
  #c2.2,
  #c3.2,
  nrow = 1,
  ncol = 3,
  align = 'hv', axis = 'lrb'
)

#-
library(patchwork)

theme_blankx <- theme(axis.title.x=element_blank(),
  #axis.text.x=element_blank(),
  #axis.ticks.x=element_blank(),
  #title=element_blank()
)

theme_blanky <- theme(axis.title.y=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  #title=element_blank()
)


theme_blankxy <- theme(axis.title.x=element_blank(),
   #axis.text.x=element_blank(),
   #axis.ticks.x=element_blank(),
   #title=element_blank(),
   axis.title.y=element_blank(),
   axis.text.y=element_blank(),
   axis.ticks.y=element_blank())

c1.1_blankx <- c1.1 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) +
  theme_blankx + labs(title = "Chromosome 1, arm 1") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", sum(dfm$chromosome_arm == "Chr1_arm1")), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black")


c2.1_blankxy <- c2.1 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) +
  theme_blankxy + labs(title = "Chromosome 2, arm 1") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", sum(dfm$chromosome_arm == "Chr2_arm1")), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black")


c3.1_blankxy <- c3.1 + scale_y_continuous(n.breaks = 5, limits=c(0, 80)) +
  theme_blankxy + labs(title = "Chromosome 3, arm 1") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", sum(dfm$chromosome_arm == "Chr3_arm1")), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black")






p_merged <- c1.1_blankx + c2.1_blankxy + c3.1_blankxy + 
  #c1.2_noblank + c2.2_blanky + c3.2_blanky +
  patchwork::plot_layout(ncol = 3, heights = c(1))

p_merged

save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_all_more500bp.png", p_merged, nrow = 2, ncol = 3, base_asp = 1.9, dpi = 500, bg = "white")
save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_all_more500bp.svg", p_merged, nrow = 2, ncol = 3, base_asp = 1.9, dpi = 500, bg = "white")


# chr4
c4 <- c4 + scale_y_continuous(n.breaks = 5, limits=c(0, 200)) + 
  #labs(title = "Chromosome 4") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", sum(dfm$chromosome_arm == "Chr4_NA")), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black")


#save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_arms_100kb_chr4_all_more500bp.png", c4, base_asp = 1.9, dpi = 500, bg = "white")
#save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_arms_100kb_chr4_all_more500bp.svg", c4, base_asp = 1.9, dpi = 500, bg = "white")


# chr1.1
c1.1_n <- c1.1 + scale_y_continuous(n.breaks = 5, limits=c(-10, 80)) + 
  #labs(title = "Chromosome 4") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", sum(dfm$chromosome_arm == "Chr1_NA")), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black")


save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_chr1.1_all_more500bp.png", c1.1_n, base_asp = 1.9, dpi = 500, bg = "white")
save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_chr1.1_all_more500bp.svg", c1.1_n, base_asp = 1.9, dpi = 500, bg = "white")


# chr2.1
c2.1_n <- c2.1 + scale_y_continuous(n.breaks = 5, limits=c(-10, 80)) + 
  #labs(title = "Chromosome 4") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", sum(dfm$chromosome_arm == "Chr2_NA")), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black")


save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_chr2.1_all_more500bp.png", c2.1_n, base_asp = 1.9, dpi = 500, bg = "white")
save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_chr2.1_all_more500bp.svg", c2.1_n, base_asp = 1.9, dpi = 500, bg = "white")


# chr3.1
c3.1_n <- c3.1 + scale_y_continuous(n.breaks = 5, limits=c(-10, 80)) + 
  #labs(title = "Chromosome 4") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", sum(dfm$chromosome_arm == "Chr3_NA")), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black")


save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_chr3.1_all_more500bp.png", c3.1_n, base_asp = 1.9, dpi = 500, bg = "white")
save_plot("~/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/DecayPlots_run12_centromeres_100kb_chr3.1_all_more500bp.svg", c3.1_n, base_asp = 1.9, dpi = 500, bg = "white")


