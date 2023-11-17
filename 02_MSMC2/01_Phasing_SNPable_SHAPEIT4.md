Phased haplotype data of single chromosomes as input for MSMC2 if using more than 2 haplotypes!

## 1. Parts of genome to be analysed

### Get names of chromosomes

    grep '^>' Chironomus_riparius_genome_010921.fasta | sed 's/>//' > chromosome-nr.txt

I will only use the chromosomes not the scaffolds

### Index bam-files

    while read f; do samtools index $f".bwamem.sort.q30.rmd.bam" > $f".bwamem.sort.q30.rmd.bam.bai" ; done < list-crip

### Extract chromosome information + index

a)  Extract chromosome

    samtools faidx Chironomus_riparius_genome_010921.fasta Chr1 > Chr1_Chironomus_riparius_genome_010921.fasta 
    samtools faidx Chironomus_riparius_genome_010921.fasta Chr2 > Chr2_Chironomus_riparius_genome_010921.fasta 
    samtools faidx Chironomus_riparius_genome_010921.fasta Chr3 > Chr3_Chironomus_riparius_genome_010921.fasta 
    samtools faidx Chironomus_riparius_genome_010921.fasta Chr4 > Chr4_Chironomus_riparius_genome_010921.fasta

b)  Index

    /home/lpettric/bin/bwa/bwa index Chr1_Chironomus_riparius_genome_010921.fasta
    /home/lpettric/bin/bwa/bwa index Chr2_Chironomus_riparius_genome_010921.fasta
    /home/lpettric/bin/bwa/bwa index Chr3_Chironomus_riparius_genome_010921.fasta
    /home/lpettric/bin/bwa/bwa index Chr4_Chironomus_riparius_genome_010921.fasta

## 2. Create mappability mask per chromosome using SNPable

a)  Extract overlapping 145mers subsequences as artificial reads from
    Chr1

I have the same data as Ann-Marie so I choose the same length

    /home/lpettric/bin/seqbility-20091110/splitfa Chr1_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000
    /home/lpettric/bin/seqbility-20091110/splitfa Chr2_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000
    /home/lpettric/bin/seqbility-20091110/splitfa Chr3_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000
    /home/lpettric/bin/seqbility-20091110/splitfa Chr4_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000


    cat xaa xab xac xad xae xaf xag > Chr1_145splits.fa
    cat xaa xab xac xad xae xaf > Chr2_145splits.fa
    cat xaa xab xac xad xae xaf > Chr3_145splits.fa
    cat xaa xab > Chr4_145splits.fa

b)  Map artificial reads back to chromosomes

    /home/lpettric/bin/bwa/bwa aln -R 1000000 -O 3 -E 3 Chr1/Chr1_Chironomus_riparius_genome_010921.fasta Chr1/145mer-subsequences/Chr1_145splits.fa > Chr1/145mer-subsequences/Chr1_145splits_bwaaln.sai

c)  Convert sai to sam

    /home/lpettric/bin/bwa/bwa samse Chr1/Chr1_Chironomus_riparius_genome_010921.fasta Chr1/145mer-subsequences/Chr1_145splits_bwaaln.sai Chr1/145mer-subsequences/Chr1_145splits.fa > Chr1_145splits_bwaaln.sam

    gzip Chr1_145splits_bwaaln.sam

d)  Generate rawMask

    gzip -dc Chr1_145splits_bwaaln.sam.gz | /home/lpettric/bin/seqbility-20091110/gen_raw_mask.pl > rawMask_Chr1_145.fa

e)  Generate the final mask

    gen_mask -l 145 -r 0.5 rawMask_Chr1_145.fa > mask_Chr1_145_50.fa

**length 145bp** **stringency 0.5**

f)  convert final-masks to .bed using makeMappabilitMask.py change paths
    of input and output

## 3. Variant calling and phasing

a)  Try script from msmc-tools (bamCaller.py)

Get coverage statistics per chromosome

    while read f; do samtools depth -r Chr1 $f".bwamem.sort.q30.rmd.bam" | awk '{sum += $3} END {print sum / NR}' > $f".Chr1.cov" ; done < list-crip

repeat for every chromosome

run bamCaller.py

    samtools mpileup -q 30 -Q 20 -C 50 -u -r <chr> -f <ref.fa> <bam> | bcftools call -c -V indels | /home/lpettric/bin/msmc-tools/bamCaller.py <mean_cov> <out_mask.bed.gz> | gzip -c > <out.vcf.gz>

    # samtools:
    # q = Minimum mapping quality for an alignment to be used
    # Q = Minimum base quality for a base to be considered.
    # C = Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50. 
    # u = uncompressed
    # r = Only generate pileup in region. Requires the BAM files to be indexed. If used in conjunction with -l then considers the intersection of the two requests.
    # f = fasta-ref
    # bcftools:
    # c = consensus-caller
    # V = skip-variants snps|indels

Get summary list of coverage per chromosome per bam

    awk '{print $0 "\t" FILENAME}' *.Chr1.cov > Chr1.summary

Final command: Instead of samtools use bcftools because mpielup migrated
to it

    module purge
    module load samtools/1.13
    module load python/3.4.3

    cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

    # Chr1
    while read -r x y; do bcftools mpileup -q 30 -Q 20 -C 50 -r Chr1 -f /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta $y".bwamem.sort.q30.rmd.bam" | bcftools call -c -V indels | /home/lpettric/bin/msmc-tools/bamCaller.py $x $y"_mask.bed.gz" | gzip -c > $y".vcf.gz" ;done < ./mean-coverage-chromosome/Chr1.summary

--\> changed script, so that it is directly bgzipped

    # Chr2
    while read -r x y; do bcftools mpileup -q 30 -Q 20 -C 50 -r Chr2 -f /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta $y".bwamem.sort.q30.rmd.bam" | bcftools call -c -V indels | /home/lpettric/bin/msmc-tools/bamCaller.py $x $y"_Chr2_mask.bed.gz" | bgzip -c > $y"_Chr2.vcf.gz" ;done < ./mean-coverage-chromosome/Chr2.summary

merge vcf.files and zip and index them (because no reference panel)

    bgzip *.vcf --> bcftools only index if they are bgzipped

normal gzip wrong! need to change to bgzip

command for already existing files of Chr1:

    for f in *vcf*; do zcat $f | bgzip -c > $f".bgz" ; done

check type of vcf with:

    htsfile file.vcf.gz

index

    for f in *vcf.gz; do bcftools index $f; done

you get csi index

    bcftools merge --print-header *.vcf.gz    # to get header info to see if there are recurring headers
    bcftools merge -O z -o merged_Chr1.vcf.gz  *.vcf.gz       

b)  No reference panel, that's why we merged the vcf, now we filter to
    only have monoallelic and biallelic SNPs

    bcftools view -M 2 -O z -o merged_biallelic_Chr1.vcf.gz
    merged_Chr1.vcf.gz bcftools index merged_biallelic_Chr1.vcf.gz

c)  Phasing with SHAPEIT

Shapeit4.2 was modified by Peter Heger to remove AVX2 dependency

start shapeit main run (need to load boost and samtools)

    /home/lpettric/bin/shapeit4/bin/shapeit4.2 -I merged_biallelic_Chr1.vcf.gz -O ./phased/merged_biallelic_Chr1_phased.vcf.gz --sequencing --region Chr1 --log  /phased/shapeit_Chr1.log

after phasing all individuals together, seperate them again

    while read f; do bcftools view -s $f -O z -o $f"_Chr1_phased.vcf.gz" merged_biallelic_Chr1_phased.vcf.gz ; done < ../../../bam-files/list-crip

vcf.gz files indexed

    for f in *.vcf.gz; do bcftools index $f; done

e)  Correct for missed genotypes

These files now contain the pashed alleles for each individual. However,
the pre-phased files might contain more information that should not be
lost.

merging phased and unphased vcfs, keeping all unphased sites from the
original vcf, but replacing the phased ones

use --force-samples because phased and unphase have same headers

    while read f; do bcftools merge --force-samples ../$f"_Chr1.vcf.gz" $f"_Chr1_phased.vcf.gz" | awk '
    BEGIN {OFS="\t"}
    $0 ~ /^##/ {print}
    $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
    $0 !~ /^#/ {
    if(substr($11, 1, 3) != "./.")
    $10 = $11
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
    }' | bcftools view -O z > $f"_Chr1_phased_merged.vcf.gz" ; done < ../../../bam-files/list-crip
