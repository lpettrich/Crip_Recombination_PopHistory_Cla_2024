## 1. Run msmc2 cross-coalescence

### Example with MF1-MG5

    # Create cross-coalesence of populations 
    # Cheops1
    module purge

    cd /projects/ag-waldvogel/pophistory/CRIP/msmc2/cross-coalescene

    # MF (within1) - MG (within2)
    /home/lpettric/bin/msmc2/build/release/msmc2 --skipAmbiguous -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0,1,2,3,4,5,6,7 -o ./run2/within1_msmc2_MF1-MG5 ../multihetsep-Chr1/multihetsep_MF1-MG5_joined_Chr1.txt ../multihetsep-Chr2/multihetsep_MF1-MG5_joined_Chr2.txt ../multihetsep-Chr3/multihetsep_MF1-MG5_joined_Chr3.txt ../multihetsep-Chr4/multihetsep_MF1-MG5_joined_Chr4.txt &
    wait
    /home/lpettric/bin/msmc2/build/release/msmc2 --skipAmbiguous -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 8,9,10,11,12,13,14,15 -o ./run2/within2_msmc2_MF1-MG5 ../multihetsep-Chr1/multihetsep_MF1-MG5_joined_Chr1.txt ../multihetsep-Chr2/multihetsep_MF1-MG5_joined_Chr2.txt ../multihetsep-Chr3/multihetsep_MF1-MG5_joined_Chr3.txt ../multihetsep-Chr4/multihetsep_MF1-MG5_joined_Chr4.txt &
    wait
    /home/lpettric/bin/msmc2/build/release/msmc2 --skipAmbiguous -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0-8,0-9,0-10,0-11,0-12,0-13,0-14,0-15,1-8,1-9,1-10,1-11,1-12,1-13,1-14,1-15,2-8,2-9,2-10,2-11,2-12,2-13,2-14,2-15,3-8,3-9,3-10,3-11,3-12,3-13,3-14,3-15,4-8,4-9,4-10,4-11,4-12,4-13,4-14,4-15,5-8,5-9,5-10,5-11,5-12,5-13,5-14,5-15,6-8,6-9,6-10,6-11,6-12,6-13,6-14,6-15,7-8,7-9,7-10,7-11,7-12,7-13,7-14,7-15 -o ./run2/across_msmc2_MF1-MG5 ../multihetsep-Chr1/multihetsep_MF1-MG5_joined_Chr1.txt ../multihetsep-Chr2/multihetsep_MF1-MG5_joined_Chr2.txt ../multihetsep-Chr3/multihetsep_MF1-MG5_joined_Chr3.txt ../multihetsep-Chr4/multihetsep_MF1-MG5_joined_Chr4.txt

### Create combined cross-coalescence

    #Cheops0
    module purge
    module load python/3.4.3
    cd /projects/ag-waldvogel/pophistory/CRIP/msmc2/cross-coalescene/run2
    # MF-MG
    /home/lpettric/bin/msmc-tools/combineCrossCoal.py across_msmc2_MF1-MG5.final.txt within1_msmc2_MF1-MG5.final.txt within2_msmc2_MF1-MG5.final.txt > combinedMF1-MG5_msmc2.final.txt

### Account for uncertainities

a)  Form mean value for estimates of every population

b)  Convert lambda and left-time-boundary to real time data 

         time <- ((data$left_time_boundary+data$right_time_boundary)/2)/mu*gen) # years ago
         time2 <- ((data$left_time_boundary+data$right_time_boundary)/2)/mu) # generations ago
         pop.size <- (1/data$lambda)/(2*mu)

         mu = 4.27*10^-9
         gen = generation time = time/generations = 1 year / 10.1 generations 
         # 10.1 = mean value generations per year from Waldvogel et al. 2018 supplement
         # use genertaion time of every population!!!


         ## IN ANALYSIS USE GENERATION TIME OF EVERY POPULATION seperate! ##


c)  Remove values with unrealistic lambda (jump in value) => first 5 values and last --> remove it constant across populations

d)  Estimate mean haplotype length (MHL) and time to most recent ancestor (tMRCA) per population

    #Cheops0 module purge module load samtools/1.13

    cd
    /projects/ag-waldvogel/pophistory/CRIP/msmc2/callable-sites/snp-call/

    # repeat variant calling without splitting data per chromosome, indels should be included to get correct number of records

    while read y; do bcftools mpileup -q 30 -Q 20 -C 50 -r Chr1,Chr2,Chr3,Chr4 -f /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta /projects/ag-waldvogel/pophistory/CRIP/bam-files/$y”.bwamem.sort.q30.rmd.bam” | bcftools call -c | bgzip -c > $y”_allChr_inclIndel.vcf.gz” ;done < /projects/ag-waldvogel/pophistory/CRIP/bam-files/list-crip & 
    wait 
    for f in *_allChr_inclIndel.vcf.gz; do bcftools index $f; done &
    wait 
    while read x; do bcftools view -M 2 -O z -o ./biallelic-snp/$x"_allChr_inclIndel_biallelic.vcf.gz" $x"_allChr_inclIndel.vcf.gz" ;done < /projects/ag-waldvogel/pophistory/CRIP/bam-files/list-crip & 
    wait 
    cd /projects/ag-waldvogel/pophistory/CRIP/msmc2/callable-sites/snp-call/biallelic-snp/ & 
    wait 
    for f in *_allChr_inclIndel_biallelic.vcf.gz; do bcftools index $f; done & 
    wait 
    while read x; do bcftools stats $x"_allChr_inclIndel_biallelic.vcf.gz" >  $x"_allChr_inclIndel_biallelic.vcf.stats" ;done < /projects/ag-waldvogel/pophistory/CRIP/bam-files/list-crip

Biallelic filtering not neccessary =\> total number of SNPs -
multiallelic SNPs = diallelic SNPs

      total SNPs - multiallelic SNPs = diallelic SNPs
      diallelic SNPs/number of records = Heterozygosity
      Mean Heterozygosity*100000 = X heterzygote position per 100000 bases
      100000 bases/X heterzygote position = 1 het per X bases (i.e. einfacher Dreisatz)
      X bases * SER = MHL
      1/(2*r*MHL) = tMRCA
