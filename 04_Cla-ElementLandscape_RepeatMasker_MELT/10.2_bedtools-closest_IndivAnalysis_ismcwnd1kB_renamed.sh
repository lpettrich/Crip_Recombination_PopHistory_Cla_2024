#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80GB
#SBATCH --time=00:01:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=bedtools_closest
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

# Load modules
module purge
module load bedtools/2.31.0 
module load R/4.3.1_system

OUTDIR=/projects/ag-waldvogel/pophistory/CRIP/comparison_Cla_rho/01-bedtools_closest/run11

INPDIR=/projects/ag-waldvogel/pophistory/CRIP/comparison_Cla_rho/02-bedtools_intersect

RHODIR=/projects/ag-waldvogel/pophistory/CRIP/ismc

IND=("MF1" "MF2" "MF3" "MF4" \
       "MG2" "MG3" "MG4" "MG5" \
       "NMF1" "NMF2" "NMF3" "NMF4" \
       "SI1" "SI2" "SI3" "SI4" \
       "SS1" "SS2" "SS3" "SS4")

CHR=("Chr1" "Chr2" "Chr3" "Chr4")

POP=("MF" "MG" "NMF" "SI" "SS")


# RHO: CONVERT BEDGRAPH TO BED AND SORT
### https://www.biostars.org/p/10141/
# Has been already done in 1.6

#for ind in "${IND[@]}"; do
#  for chr in "${CHR[@]}"; do    
#    cd ${RHODIR}/"${ind}"/"${chr}"
#    cat "${ind}"_"${chr}"_ismc.rho.1kb.bedgraph | grep -v '^chrom' \
#    | sed "s/chr1/${chr}/g" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.bed
#    sort -k1,1 -k2,2n ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.bed > ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.sorted.bed
#  done
#done

wait

echo "Converting bedgraph to bed done"

wait

# NO MEAN RHO/10kb FOR EACH POP NEED If USING IndivAnalysis


# CLA: USE IndivAnalysis of AllPop Run
# SPLIT CLA-BED BY CHROMOSOME AND SORT
# Already done during renaming

echo "Split Cla-bed by chromosome and sort it done"

wait

# BEDTOOLS CLOSEST
### with pophistory samples mapped on Crip4.0
### for each individual and each chromosome 
### Cla regions compared to rho (1 kb)

cd ${OUTDIR}

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do
    bedtools closest -d -t first -a ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.sorted.bed \
    -b ${INPDIR}/"${ind}"_"${chr}".bwamem.sambl.sort.Cla1.tmp_renamed.final.sorted.bed > \
    "${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.sorted_distToNearestCla1.bed
  done
done

wait

echo "Bedtools closest finished"

wait

# FILTER FOR CLUSTERSIZE MORE OR EQUAL THAN 500 bp

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do
    cd ${OUTDIR}
    awk '$7 - $6 >= 500' "${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.sorted_distToNearestCla1.bed > "${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.sorted_distToNearestCla1_more500bp.bed
  done
done

# FILTER FOR CLUSTERSIZE LESS THAN 500 bp

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do
    cd ${OUTDIR}
    awk '$7 - $6 < 500' "${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.sorted_distToNearestCla1.bed > "${ind}"_"${chr}"_ismc.rho.1kb.bedgraph_renamed.sorted_distToNearestCla1_less500bp.bed
  done
done

