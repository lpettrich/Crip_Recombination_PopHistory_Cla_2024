# Clean environment
rm(list = ls())

# Load libraries
 library(ggplot2)
 library(reshape2)
 library(RColorBrewer) 
 library(sf)
 library(rnaturalearth)
 library(rnaturalearthdata)
 library(rgeos)
 library(maptools)
 library(rgdal)
 library(raster)
 library(scales)
 library(ggpubr)
 library(egg)
 library(cowplot)
 library(gridGraphics)
 library(grid)
 library(tidyverse)

# Set directory
setwd("~/sciebo/final-plots/")
getwd()

#-------------------------------------------------------------------------------------------------------------
############
# 1. MSMC2 #
############

# A) EFFECTIVE POPULATION SIZE
#-----------------------------
# Read in data
# Combined Cross-Coalescence; mean values per population
# Read in data
mf <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/Mean-values-msmc2-MF.csv", header = TRUE, sep = ";")
mg <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/Mean-values-msmc2-MG.csv", header = TRUE, sep = ";")
nmf <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/Mean-values-msmc2-NMF.csv", header = TRUE, sep = ";")
si <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/Mean-values-msmc2-SI.csv", header = TRUE, sep = ";")
ss <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/Mean-values-msmc2-SS.csv", header = TRUE, sep = ";")

# Check data
str(mf)

# determine mutation rate and generation time
mu <- 4.27*10^-9    #from Waldvogel & Pfenninger 2021; mutation rate
gen_mf <- 1/9.07    # generation time = time/generations
gen_mg <-  1/7.85 
gen_nmf <- 1/7.7
gen_si <- 1/10.57
gen_ss <- 1/14.86

# remove outliers in data
mf <- mf[6:31,]
mg <- mg[6:31,]
nmf <- nmf[6:31,]
si <- si[6:31,]
ss <- ss[6:31,]

# merge data frames
mdata <- rbind(mf[,29:32],mg[,29:32],nmf[,29:32],si[,29:32],ss[,29:32])

# add column with population names
# mdata <- data.frame(mdata,Population = rep(c("MF (Rhône-Alpes, FRA)","MG (Hessen, GER)","NMF (Lorraine, FRA)","SI (Piemonte, ITA)","SS (Andalucia, ESP)"),
#                                           times=c(nrow(mf),nrow(mg),nrow(nmf),nrow(si),nrow(ss))))

mdata <- data.frame(mdata,Population = rep(c("MF","MG","NMF","SI","SS"),
                                           times=c(nrow(mf),nrow(mg),nrow(nmf),nrow(si),nrow(ss))))

# A.1) YEARS AGO     
# define colour palette
my_palette1 <- "#3AB6D8"
my_palette2 <- "#CD70C6"
my_palette3 <- "#4C83EB"
my_palette4 <- "#FF941A"
my_palette5 <- "#AB3232"
my_palette <- c(my_palette1, my_palette2, my_palette3, my_palette4, my_palette5)

mdata$Population <- factor(mdata$Population ,                 # Relevel group factor
                     levels = c("MG", "NMF", "MF", "SI", "SS"))

# plot Ne against time
scientific_10 <- function(x) {   parse(text=gsub("e\\+*", " %*% 10^", 
                                                 scales::scientific_format()(x))) }

ggp <- ggplot(mdata, aes(x = mean_years, y = mean_Ne, colour = Population)) + 
  geom_step(size=0.75)  +
  scale_x_continuous(limits=c(10^2,10^6), trans="log10", labels = scientific_10) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = my_palette)  +
  labs(x = "Years ago", y = "Effective \npopulation size") + 
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14,face = "bold", color = "black"),
        axis.title.x = element_text(size = 14,face = "bold", color = "black"),
        title = element_text(size = 15, color = "black"),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line("lightgrey"),
        panel.grid.minor = element_line("white"),
        panel.border = element_rect("black", fill = NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.position="none") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
ggp



pNeYears <-  ggp


# B) CROSS-COALESCENCE RATE
#---------------------------
# Read in data
# Combined Cross-Coalescence
mfmg <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_MF1-MG5_msmc2.final.txt", header = TRUE)
mfnmf <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_MF1-NMF4_msmc2.final.txt", header = TRUE)
mfsi <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_MF1-SI4_msmc2.final.txt", header = TRUE)
mfss <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_MF1-SS4_msmc2.final.txt", header = TRUE)
mgnmf <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_MG2-NMF4_msmc2.final.txt", header = TRUE)
mgsi <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_MG2-SI4_msmc2.final.txt", header = TRUE)
mgss <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_MG2-SS4_msmc2.final.txt", header = TRUE)
nmfsi <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_NMF1-SI4_msmc2.final.txt", header = TRUE)
nmfss <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_NMF1-SS4_msmc2.final.txt", header = TRUE)
siss <- read.table("~/sciebo/MasterThesis/MSMC2/cross-coal/combined-crosscoal-run2/combined_SI1-SS4_msmc2.final.txt", header = TRUE)


# mutation rate and generation time defined in previous step

# to get relative gene flow, you can compute the relative cross-coalescence rate: 2 * lambda01 / (lambda00 + lambda11)

str(mfmg)
mfmg <- mfmg %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

mfnmf <- mfnmf %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

mfsi <- mfsi %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

mfss <- mfss %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

mgnmf <- mgnmf %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

mgsi <- mgsi %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

mgss <- mgss %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

nmfsi <- nmfsi %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

nmfss <- nmfss %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))

siss <- siss %>% 
  mutate(rel.cc = (2*lambda_01)/(lambda_00+lambda_11))



# add column with pop.size, generations and years
mfmg <- mfmg %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),    # generations ago
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_mf+gen_mg)/2)))   # multiple by generation time to get years

mfnmf <- mfnmf %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),    
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_mf+gen_nmf)/2)))   

mfsi <- mfsi %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),    
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_mf+gen_si)/2)))   

mfss <- mfss %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),    
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_mf+gen_ss)/2)))   

mgnmf <- mgnmf %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_mg+gen_nmf)/2)))

mgsi <- mgsi %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_mg+gen_si)/2)))

mgss <- mgss %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_mg+gen_ss)/2)))

nmfsi <- nmfsi %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_nmf+gen_si)/2)))

nmfss <- nmfss %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_nmf+gen_ss)/2)))

siss <- siss %>% 
  mutate(generations = (((left_time_boundary+right_time_boundary)/2)/mu),
         years = (((left_time_boundary+right_time_boundary)/2)/mu*((gen_si+gen_ss)/2)))



# remove outliers in data
mfmg <- mfmg[6:31,]
mfnmf <- mfnmf[6:31,]
mfsi <- mfsi[6:31,]
mfss <- mfss[6:31,]
mgnmf <- mgnmf[6:31,]
mgsi <- mgsi[6:31,]
mgss <- mgss[6:31,]
nmfsi <- nmfsi[6:31,]
nmfss <- nmfss[6:31,]
siss <- siss[6:31,]

# merge data frames
mdata <- rbind(mfmg,mfnmf,mfsi,mfss,mgnmf,mgsi,mgss,nmfsi,nmfss,siss)

# add column with population names
mdata <- data.frame(mdata,Population = rep(c("MF-MG","MF-NMF","MF-SI","MF-SS",
                                             "MG-NMF","MG-SI","MG-SS",
                                             "NMF-SI","NMF-SS","SI-SS"),
                                           times=c(nrow(mfmg),nrow(mfnmf),nrow(mfsi),nrow(mfss),nrow(mgnmf),nrow(mgsi),nrow(mgss),
                                                   nrow(nmfsi),nrow(nmfss),nrow(siss))))


# B.1) CROSS-COALESCENCE - YEARS AGO
# plot cross-coal against time
# scale already defined in previous step

# define palette
my_palette <- colorRampPalette(c("firebrick4", "indianred2", "orange", "moccasin", "mediumpurple3", "skyblue2", "royalblue"))(10)
my_palette <- rev(my_palette)

# plot
ggp2 <- ggplot(mdata, aes(x = years, y = rel.cc, colour = Population)) + 
  geom_step(size=0.75)  +
  geom_hline(yintercept = 0.5, colour="gray27", linetype = "dashed") +
  scale_x_continuous(limits=c(10^2,10^6), trans="log10", labels = scientific_10) +
  scale_y_continuous(limits = c(0,1.2), labels = comma, breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = my_palette)  +
  labs(x = "Years ago", y = "Relative \ncross-coalescence \nrate") + 
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14,face = "bold", color = "black"),
        axis.title.x = element_text(size = 14,face = "bold", color = "black"),
        title = element_text(size = 15, color = "black"),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line("lightgrey"),
        panel.grid.minor = element_line("white"),
        panel.border = element_rect("black", fill = NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_rect("white"),
        legend.position="top") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE, title= "Population\npair"))
ggp2
#---------------------------------------------------------------------------------------------
# It's better to use annual mean temperature https://chelsa-climate.org/chelsa-trace21k/
# bio1 = Annual Mean Temperature [°C]
# Read in tif 
setwd("~/Documents/03_Master/Thesis/ClimateData/CHELSA/bio1")
rlist=list.files(getwd(), pattern=".tif", full.names=FALSE)
for(i in rlist) { assign(unlist(strsplit(i, "[.]"))[1], raster(i)) } 

#kBP1 <- raster("CHELSA_TraCE21k_bio01_10_V1.0.tif")
#kBP2 <- raster("CHELSA_TraCE21k_bio01_0_V1.0.tif")
#kBP3 <- raster("CHELSA_TraCE21k_bio01_-10_V1.0.tif")
#kBP4 <- raster("CHELSA_TraCE21k_bio01_-20_V1.0.tif")
#kBP5 <- raster("CHELSA_TraCE21k_bio01_-30_V1.0.tif")
#kBP6 <- raster("CHELSA_TraCE21k_bio01_-40_V1.0.tif")
#kBP7 <- raster("CHELSA_TraCE21k_bio01_-50_V1.0.tif")
#kBP8 <- raster("CHELSA_TraCE21k_bio01_-60_V1.0.tif")
#kBP9 <- raster("CHELSA_TraCE21k_bio01_-70_V1.0.tif")
#kBP10 <- raster("CHELSA_TraCE21k_bio01_-80_V1.0.tif")
#kBP11 <- raster("CHELSA_TraCE21k_bio01_-90_V1.0.tif")
#kBP12 <- raster("CHELSA_TraCE21k_bio01_-100_V1.0.tif")
#kBP13 <- raster("CHELSA_TraCE21k_bio01_-110_V1.0.tif")
#kBP14 <- raster("CHELSA_TraCE21k_bio01_-120_V1.0.tif")
#kBP15 <- raster("CHELSA_TraCE21k_bio01_-130_V1.0.tif")
#kBP16 <- raster("CHELSA_TraCE21k_bio01_-140_V1.0.tif")
#kBP17 <- raster("CHELSA_TraCE21k_bio01_-150_V1.0.tif")
#kBP18 <- raster("CHELSA_TraCE21k_bio01_-160_V1.0.tif")
#kBP19 <- raster("CHELSA_TraCE21k_bio01_-170_V1.0.tif")
#kBP20 <- raster("CHELSA_TraCE21k_bio01_-180_V1.0.tif")
#kBP21 <- raster("CHELSA_TraCE21k_bio01_-190_V1.0.tif")
kBP22 <- raster("CHELSA_TraCE21k_bio01_-200_V1.0.tif")

GDALinfo("CHELSA_TraCE21k_bio01_-200_V1.0.tif")


# Create a data.frame with sample site coordinates
site <- c("MG","NMF", "MF", "SI", "SS")
lon <- c(9.0819270 , 6.2156670, 4.8865000 , 8.3473320, -4.5267980)
lat <- c(50.1680610, 49.1765430, 45.8616760, 45.4036180, 37.399080)
samples <- data.frame(site, lon, lat, row.names="site")
samples


plot(kBP22)
points(samples, pch = 16)

# Extract data from tif for your sites
ID <- seq(-200, 10, by=10)
d <- c(paste("CHELSA_TraCE21k_bio01_", ID, "_V1", sep = ""))


temp.data <- samples 
temp.data$kBP1 <- raster::extract(`CHELSA_TraCE21k_bio01_10_V1`, samples)
temp.data$kBP2 <- raster::extract(`CHELSA_TraCE21k_bio01_0_V1`, samples)
temp.data$kBP3 <- raster::extract(`CHELSA_TraCE21k_bio01_-10_V1`, samples)
temp.data$kBP4 <- raster::extract(`CHELSA_TraCE21k_bio01_-20_V1`, samples)
temp.data$kBP5 <- raster::extract(`CHELSA_TraCE21k_bio01_-30_V1`, samples)
temp.data$kBP6 <- raster::extract(`CHELSA_TraCE21k_bio01_-40_V1`, samples)
temp.data$kBP7 <- raster::extract(`CHELSA_TraCE21k_bio01_-50_V1`, samples)
temp.data$kBP8 <- raster::extract(`CHELSA_TraCE21k_bio01_-60_V1`, samples)
temp.data$kBP9 <- raster::extract(`CHELSA_TraCE21k_bio01_-70_V1`, samples)
temp.data$kBP10 <- raster::extract(`CHELSA_TraCE21k_bio01_-80_V1`, samples)
temp.data$kBP11 <- raster::extract(`CHELSA_TraCE21k_bio01_-90_V1`, samples)
temp.data$kBP12 <- raster::extract(`CHELSA_TraCE21k_bio01_-100_V1`, samples)
temp.data$kBP13 <- raster::extract(`CHELSA_TraCE21k_bio01_-110_V1`, samples)
temp.data$kBP14 <- raster::extract(`CHELSA_TraCE21k_bio01_-120_V1`, samples) 
temp.data$kBP15 <- raster::extract(`CHELSA_TraCE21k_bio01_-130_V1`, samples)
temp.data$kBP16 <- raster::extract(`CHELSA_TraCE21k_bio01_-140_V1`, samples)
temp.data$kBP17 <- raster::extract(`CHELSA_TraCE21k_bio01_-150_V1`, samples)
temp.data$kBP18 <- raster::extract(`CHELSA_TraCE21k_bio01_-160_V1`, samples)
temp.data$kBP19 <- raster::extract(`CHELSA_TraCE21k_bio01_-170_V1`, samples)
temp.data$kBP20 <- raster::extract(`CHELSA_TraCE21k_bio01_-180_V1`, samples)
temp.data$kBP21 <- raster::extract(`CHELSA_TraCE21k_bio01_-190_V1`, samples)
temp.data$kBP22 <- raster::extract(`CHELSA_TraCE21k_bio01_-200_V1`, samples)

fix("temp.data")

# new data fram with only temperature
temp <- as.data.frame(t(temp.data[,-c(1:2)]))
temp$kBP <- c(1:22)



d <- melt(temp, id.vars="kBP")

d$variable <- factor(d$variable ,                 # Relevel group factor
                     levels = c("MG", "NMF", "MF", "SI", "SS"))

scientific_10 <- function(x) {   parse(text=gsub("e\\+*", " %*% 10^", 
                                                 scales::scientific_format()(x))) }

d$kBPtrans <- d$kBP*1000

my_palette1 <- "#3AB6D8"
my_palette3 <- "#4C83EB"
my_palette2 <- "#CD70C6"
my_palette4 <- "#FF941A"
my_palette5 <- "#AB3232"
my_palette <- c(my_palette1, my_palette2, my_palette3, my_palette4, my_palette5)

# Everything on the same plot
p <- ggplot(d, aes(kBPtrans,value, col=variable, group = variable)) + 
  geom_point() + geom_line() +
  scale_x_continuous(limits=c(10^2,10^6), trans="log10", labels = scientific_10) +
  scale_y_continuous(limits = c(-7, 20), n.breaks = 7) +
  scale_color_manual(values = my_palette)  +
  labs(x = "Years ago", y = "Temperature (°C)") + 
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14,face = "bold", color = "black"),
        axis.title.x = element_text(size = 14,face = "bold", color = "black"),
        title = element_text(size = 15, color = "black"),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line("lightgrey"),
        panel.grid.minor = element_line("white"),
        panel.border = element_rect("black", fill = NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_rect("white"),
        legend.position="bottom") + 
  guides(color = guide_legend(nrow = 1, byrow = TRUE, title= "Population"))

p

# SAVE PLOT
#svg("../../CHELSA-Climate-kBP22-1_line.longaxis.svg", width=10,height=6)
#p
#dev.off()


#-----------------------------------------------------------------------------------------------
# Create a map of Europe with sampling locations
# Crop raster to only have the coordinates of Europe
plot(`CHELSA_TraCE21k_bio01_-200_V1`)

europe <- c(-15, 20, 30, 60)
kBP22 <- crop(`CHELSA_TraCE21k_bio01_-200_V1`, europe)
plot(kBP22)

kBP10 <- crop(`CHELSA_TraCE21k_bio01_-80_V1`, europe)
plot(kBP10)

kBP1 <- crop(`CHELSA_TraCE21k_bio01_-10_V1`, europe)
plot(kBP1)

# Load map to have a map with land mass underneath
m <- raster("../../GRAY_50M_SR_W/GRAY_50M_SR_W.tif") # https://www.naturalearthdata.com/downloads/50m-raster-data/50m-gray-earth/
plot(m)

# Crop map to coordinates of Europe
m <- crop(m,europe)
plot(m)

# Load political map
data(wrld_simpl)
wrld_simpl
wrld_simpl_c <- crop(wrld_simpl,europe)
wrld_simpl_c


my_palette <- c("#3AB6D8", "#CD70C6", "#4C83EB","#FF941A", "#AB3232")
my_palette2 <- brewer.pal(n = 9, name = "Greys")
my_palette3 <- brewer.pal(n = 11, name = "RdYlBu")[c(11:8,6,5:1)]


# PLOT kBP22
#-------------------------------------------------
#sizing window to a good shape for the world
dev.new(width=10,height=5)
#so that maps extends to edge of window
oldpar <- par(mai=c(0,0,.5,0.1))
plot(m,  
     col = my_palette2,
     main = "22k-BP",
     axes = FALSE,
     legend = FALSE
)
plot(kBP22, alpha = 0.7, add = TRUE, col = my_palette3, breaks = c(-15,-10,-5,0,5,10,15,20,25))
points(samples, pch = 21, col = "black", bg = my_palette)
#plot(wrld_simpl_c,add=TRUE)
#dev.off()
par(oldpar) #reset graphics settings

# Save base R plots as objects
## grab the scene as a grid object

grid.echo()
a <- grid.grab()

## draw it, changes optional
grid.newpage()
a <- editGrob(a)
grid.draw(a)

# PLOT kBP10
#---------------------------------------------
#sizing window to a good shape for the world
dev.new(width=10,height=5)
#so that maps extends to edge of window
oldpar <- par(mai=c(0,0,.5,0.1))
plot(m,  
     col = my_palette2,
     main = "10k-BP",
     axes = FALSE,
     legend = FALSE,
)
plot(kBP10, alpha = 0.7, add = TRUE, col = my_palette3, breaks = c(-15,-10,-5,0,5,10,15,20,25))
points(samples, pch = 21, col = "black", bg = my_palette)
#plot(wrld_simpl_c,add=TRUE)
#dev.off()
par(oldpar) #reset graphics settings


# Save base R plots as objects
## grab the scene as a grid object

grid.echo()
b <- grid.grab()

## draw it, changes optional
grid.newpage()
b <- editGrob(b)
grid.draw(b)

# PLOT kBP1
#--------------------------------------------
#sizing window to a good shape for the world
dev.new(width=10,height=5)
#so that maps extends to edge of window
oldpar <- par(mai=c(0,0,.5,0.1))
plot(m,  
     col = my_palette2,
     main = "1k-BP",
     axes = FALSE,
     legend = FALSE,
)
plot(kBP1, alpha = 0.7, add = TRUE, col = my_palette3, breaks = c(-15,-10,-5,0,5,10,15,20,25))
points(samples, pch = 21, col = "black", bg = my_palette)
#plot(wrld_simpl_c,add=TRUE)
#dev.off()
par(oldpar) #reset graphics settings



# Save base R plots as objects
## grab the scene as a grid object

grid.echo()
c <- grid.grab()

## draw it, changes optional
grid.newpage()
c <- editGrob(c)
grid.draw(c)


dev.off()


#------------------------------------------------------------------------------------
#Merge plots

p <- p + theme(legend.position="none")
legend_a <- get_legend(ggp2 + theme(legend.position="top", 
                                    legend.text=element_text(size = 14, color = "black")))
legend_b <- get_legend(ggp + theme(legend.position="bottom", 
                                   legend.text=element_text(size = 14, color = "black")))

p <- p +
  scale_y_continuous(limits = c(-7, 18), n.breaks = 7) +
  scale_x_continuous(limits=c(3*10^2,10^6), trans="log10", labels = scientific_10)  +
  geom_vline(xintercept = 1*10^3, colour="gray27", linetype = "dashed") +
  geom_vline(xintercept = 10*10^3, colour="gray27", linetype = "dashed") +
  geom_vline(xintercept = 22*10^3, colour="gray27", linetype = "dashed") 

ggp <- ggp +
  scale_x_continuous(limits=c(3*10^2,10^6), trans="log10", labels = scientific_10) +
  geom_vline(xintercept = 1*10^3, colour="gray27", linetype = "dashed") +
  geom_vline(xintercept = 10*10^3, colour="gray27", linetype = "dashed") +
  geom_vline(xintercept = 22*10^3, colour="gray27", linetype = "dashed") 

ggp2 <- ggp2 +
  scale_x_continuous(limits=c(3*10^2,10^6), trans="log10", labels = scientific_10)  +
  geom_vline(xintercept = 1*10^3, colour="gray27", linetype = "dashed") +
  geom_vline(xintercept = 10*10^3, colour="gray27", linetype = "dashed") +
  geom_vline(xintercept = 22*10^3, colour="gray27", linetype = "dashed") 

pcol <- NULL
pcol <- plot_grid( legend_a,
                   ggp2 + theme(legend.position="none"),
                   ggp + theme(legend.position="none"),
                   p + theme(legend.position="none"),
                   legend_b,
                   labels = c("", "A", "B", "C", ""),
                   align = 'hv', axis = 'lrb',
                   hjust = -.8,
                   vjust = -.5,
                   ncol= 1,
                   label_size = 16,
                   rel_heights = c(0.7, 1, 1, 1, 0.7)
)

pcol

#pcol <- plot_grid(legend_a, pcol, legend_b, align = 'hv', axis = 'lrb', ncol = 1, rel_heights = c(0.1, 1, .05))

pcol

#svg("~/sciebo/MasterThesis/final-plots/MSMC2-Temperature-combined-20221029_LRUA22.svg", height = 10, width = 8)
#pcol
#dev.off()

#svg("../../MSMC2-Temperature-combined2_20221029_LRUA22.svg", width = 8, height = 10)
gA <- ggplotGrob(ggp2)
gB <- ggplotGrob(ggp)
gC <- ggplotGrob(p)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB, gC))
dev.off()

#library(patchwork)
#wrap_elements(a) + wrap_elements(b) + wrap_elements(c) 

# SAVE FOR MANUSCIPT
bottom_row <- plot_grid(c, 
                        b, 
                        a, 
                        labels = c("D", "E", "F"), 
                        hjust = -.8,
                        vjust = -.5,
                        label_size = 16, ncol = 3
)  + theme(plot.margin = unit(c(20,20,10,10), "points"))

bottom_row


theme_set(theme_cowplot(font_size=16))
all_plots <- 
  plot_grid(pcol,
           bottom_row,
           align = 'hv', axis = 'lrb',
           ncol= 1,
           rel_heights = c(6, 3),
           rel_widths = c(1, 1),
           scale = c(1,1) 
           ) 


plot(all_plots)

#library(patchwork)
#pcol + wrap_elements(bottom_row) + plot_layout(heights = c(5,2))

#gridExtra::grid.arrange(pcol, bottom_row, padding = unit(0.5, "line"), nrow = 2, 
#                        widths = c(1),
#                        heights = c(5,1.5))

#ggsave("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/png-files/MSMC2-climate-maps-CRIP-Years-combined_Manuscript_20230911.png", 
#       all_plots, height = 1.1*10, width = 1.5*10, bg = "white")

save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/MSMC2-climate-maps-CRIP-Years-combined_Manuscript_20230911.svg", 
          all_plots, nrow = 7, ncol = 4, base_asp = 1.7, dpi = 500, bg = "white", scale = 0.57)

save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/MSMC2-climate-maps-CRIP-Years-combined_Manuscript_20230911.png", 
          all_plots, nrow = 7, ncol = 4, base_asp = 1.7, dpi = 500, bg = "white", scale = 0.57)

# SAVE SEPERATE AND MERGE WITH INKSCAPE
pcol
save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/MSMC2-CRIP-Years_Manuscript_20240227.svg", 
          pcol, nrow = 3, ncol = 2, base_asp = 1.45, dpi = 500, bg = "white", scale = 1)



bottom_row
save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/Plots/Cimate-Maps_Manuscript_20240227.svg", 
          bottom_row, nrow = 2, ncol = 3, base_asp = 1.9, dpi = 500, bg = "white", scale = 0.8)
