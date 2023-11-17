
## 1. Create input-files (multihetsep)

    while read a b c d; do /home/lpettric/bin/msmc-tools//generate_multihetsep.py --mask=../$a"_Chr1_mask.bed.gz" \
                          --mask=..$b"_Chr1_mask.bed.gz" \
                          --mask=..$c"_Chr1_mask.bed.gz" \
                          --mask=..$d"_Chr1_mask.bed.gz" \
                          --mask=/projects/ag-waldvogel/pophistory/CRIP/masking/final-mask/mask_Chr1_145_50.bed.gz \
                          $a"_Chr1_phased_merged.vcf.gz" $b"_Chr1_phased_merged.vcf.gz" \
                          $c"_Chr1_phased_merged.vcf.gz" $d"_Chr1_phased_merged.vcf.gz" > "multihetsep_"$a"_"$b"_"$c"_"$d"_Chr1.txt"; done < ../../../msmc2/list-populations

## 2. Create input file for cross-coalescence

Create multihetsep with 16 haplotypes

        #Cheops0
        module purge
        module load python/3.4.3
        cd /projects/ag-waldvogel/pophistory/CRIP/phasing/Chr1/phased
        # 4 inidividuals per population
        while read a b c d e f g h; do /home/lpettric/bin/msmc-tools//generate_multihetsep.py --mask=../$a"_Chr1_mask.bed.gz" \
                          --mask=../$b"_Chr1_mask.bed.gz" \
                          --mask=../$c"_Chr1_mask.bed.gz" \
                          --mask=../$d"_Chr1_mask.bed.gz" \
                          --mask=../$e"_Chr1_mask.bed.gz" \
                          --mask=../$f"_Chr1_mask.bed.gz" \
                          --mask=../$g"_Chr1_mask.bed.gz" \
                          --mask=../$h"_Chr1_mask.bed.gz" \
                          --mask=/projects/ag-waldvogel/pophistory/CRIP/masking/final-mask/mask_Chr1_145_50.bed.gz \
                          $a"_Chr1_phased_merged.vcf.gz" $b"_Chr1_phased_merged.vcf.gz" \
                          $c"_Chr1_phased_merged.vcf.gz" \
                          $d"_Chr1_phased_merged.vcf.gz" $e"_Chr1_phased_merged.vcf.gz" $f"_Chr1_phased_merged.vcf.gz" \
                          $g"_Chr1_phased_merged.vcf.gz" \
                          $h"_Chr1_phased_merged.vcf.gz"> /projects/ag-waldvogel/pophistory/CRIP/msmc2/multihetsep-Chr1/"multihetsep_"$a"-"$h"_joined_Chr1.txt"; done < /projects/ag-waldvogel/pophistory/CRIP/msmc2/list-populations-cc

For list-populations-cc see below:

list-populations-cc:

MF1 MF2 MF3 MF4 MG2 MG3 MG4 MG5

MF1 MF2 MF3 MF4 NMF1 NMF2 NMF3 NMF4

MF1 MF2 MF3 MF4 SI1 SI2 SI3 SI4

MF1 MF2 MF3 MF4 SS1 SS2 SS3 SS4

MG2 MG3 MG4 MG5 NMF1 NMF2 NMF3 NMF4

MG2 MG3 MG4 MG5 SI1 SI2 SI3 SI4

MG2 MG3 MG4 MG5 SS1 SS2 SS3 SS4

NMF1 NMF2 NMF3 NMF4 SI1 SI2 SI3 SI4

NMF1 NMF2 NMF3 NMF4 SS1 SS2 SS3 SS4

SI1 SI2 SI3 SI4 SS1 SS2 SS3 SS4
