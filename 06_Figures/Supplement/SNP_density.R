# SNP density

# Clean environment
rm(list = ls())

# Load necessary libraries
library(vcfR)
library(dplyr)
library(ggplot2)
library(scales)



#---------------------------------------------------------------------------------------------------------
#   BASED ON VCF
#---------------------------------------------------------------------------------------------------------
# Set directory
setwd("~/Documents/03_Master/Thesis/CRIP/phased-vcfs/")

# Step 1: Read the VCF file
# Replace 'your_file.vcf.gz' with the path to your VCF file
vcf <- read.vcfR("MF1_Chr1_phased_merged.vcf.gz")  # Can also read .vcf.gz files

# Step 2: Extract the necessary data from the VCF file
# Extract fixed data (CHROM, POS)
vcf_data <- as.data.frame(vcf@fix)

# Step 3: Create a data frame with chromosome and position for SNPs
snps <- vcf_data %>%
  filter(REF != "." & ALT != ".") %>%  # Ensure it's a valid SNP
  select(CHROM, POS) %>%                 # Keep only chromosome and position
  mutate(POS = as.numeric(POS))        # Ensure position is numeric

# Step 4: Define the correct order of chromosomes
goodChrOrder <- c(paste("Chr", c(1:4), sep = ""))
snps$CHROM <- factor(snps$CHROM, levels = goodChrOrder)

# Step 5: Plot the SNP density using ggplot2
snpDensity <- ggplot(snps, aes(x = POS)) + 
  geom_histogram(binwidth = 10e3, color = "brown4", fill = "brown3") +  # Adjust bin width as needed
  facet_wrap(~ CHROM, ncol = 2, scales = "free_x") +  # Separate plots for each chromosome
  ggtitle("Chromosome-wise SNP Distribution") + 
  xlab("Genomic location (bp)") + 
  ylab("SNP density") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

# Step 6: Customize the y-axis text size for better visibility
snpDensity + theme(axis.text.y = element_text(size = 6))



#---------------------------------------------------------------------------------------------------------
#   BASED ON MULTIHETSEP
#---------------------------------------------------------------------------------------------------------
# Set directory
setwd("~/Documents/03_Master/Thesis/CRIP/multihetsep-files/")

# Step 1: Read the multihetsep file
# Replace 'your_file.tabular' with the path to your multihetsep file
snps1 <- read.table("multihetsep_MF1_MF2_MF3_MF4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                   comment.char = "#")
snps2 <- read.table("multihetsep_MF1_MF2_MF3_MF4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps3 <- read.table("multihetsep_MF1_MF2_MF3_MF4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps4 <- read.table("multihetsep_MF1_MF2_MF3_MF4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps5 <- read.table("multihetsep_MG2_MG3_MG4_MG5_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps6 <- read.table("multihetsep_MG2_MG3_MG4_MG5_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps7 <- read.table("multihetsep_MG2_MG3_MG4_MG5_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps8 <- read.table("multihetsep_MG2_MG3_MG4_MG5_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps9 <- read.table("multihetsep_NMF1_NMF2_NMF3_NMF4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps10 <- read.table("multihetsep_NMF1_NMF2_NMF3_NMF4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps11 <- read.table("multihetsep_NMF1_NMF2_NMF3_NMF4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps12 <- read.table("multihetsep_NMF1_NMF2_NMF3_NMF4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps13 <- read.table("multihetsep_SI1_SI2_SI3_SI4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps14 <- read.table("multihetsep_SI1_SI2_SI3_SI4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps15 <- read.table("multihetsep_SI1_SI2_SI3_SI4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps16 <- read.table("multihetsep_SI1_SI2_SI3_SI4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps17 <- read.table("multihetsep_SS1_SS2_SS3_SS4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps18 <- read.table("multihetsep_SS1_SS2_SS3_SS4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps19 <- read.table("multihetsep_SS1_SS2_SS3_SS4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps20 <- read.table("multihetsep_SS1_SS2_SS3_SS4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

# Put them together
mf <- rbind(snps1,snps2,snps3,snps4)
mg <- rbind(snps5, snps6,snps7,snps8)
nmf <- rbind(snps9,snps10,snps11,snps12)
si <- rbind(snps13, snps14,snps15,snps16)
ss <- rbind(snps17, snps18,snps19,snps20)

# Add population label
mf$pop <- "MF"
mg$pop <- "MG"
nmf$pop <- "NMF"
si$pop <- "SI"
ss$pop <- "SS"

# Bind them all
snps <- rbind(mf,mg,nmf,si,ss)

# Step 2: Assign column names
colnames(snps) <- c("chr", "start", "id", "refallele", "pop")

# Step 3: Filter out rows with invalid SNP data (if necessary)
# Here, we are keeping all rows as there is no filtering criterion mentioned
# If needed, you can add conditions to filter SNPs based on specific criteria

# Step 4: Define the correct order of chromosomes and populations
goodChrOrder <- c(paste("Chr", c(1:4), sep = ""))
snps$chr <- factor(snps$chr, levels = goodChrOrder)

goodPopOrder <- c("MG", "NMF", "MF", "SI", "SS")
snps$pop <- factor(snps$pop, levels = goodPopOrder)

# Step 5: Plot the SNP density using ggplot2
snpDensity <- ggplot(snps, aes(x = start/1e6,  fill = pop)) + 
  geom_histogram(binwidth = 50e3/1e6) +  # Adjust bin width as needed
  facet_grid(pop ~ chr, scales = "free_x") +  # Separate plots for each chromosome
  #ggtitle("Chromosome-wise SNP distribution based on multihetsep files") + 
  xlab("Genomic location (Mb)") + 
  ylab("SNP density") + 
  labs(fill = "Population") +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16,face = "bold", color = "black"),
        axis.title.x = element_text(size = 16,face = "bold", color = "black"),
        title = element_text(size = 16, color = "black"),
        text = element_text(size = 16, color = "black"),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line("white"),
        panel.grid.minor = element_line("white"),
        panel.grid.major.x = element_line("lightgrey", linetype = 3),
        panel.grid.major.y = element_line("lightgrey", linetype = 3),
        panel.border = element_rect("black", fill = NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.position="bottom",
        #strip.text.y = element_text(angle = 0),
        strip.text.y = element_blank())


my_palette <- c("MF" = "#4C83EB",   # Assign colors to populations
                "NMF" = "#CD70C6",
                "MG" = "#3AB6D8",
                "SI" = "#FF941A",
                "SS" = "#AB3232")

snpDensity <- snpDensity + scale_fill_manual(values = my_palette) + scale_x_continuous(label = comma)

snpDensity

cowplot::save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/SNP_density.svg", snpDensity, nrow = 5, ncol = 4, base_asp = 1.6, dpi = 500, bg = "white", scale = 0.8)
cowplot::save_plot("/home/alle/sciebo/RecombinationLandscape_CRIP/01_Revison/NewResults/Plots/SNP_density.png", snpDensity, nrow = 5, ncol = 4, base_asp = 1.6, dpi = 500, bg = "white", scale = 0.8)

