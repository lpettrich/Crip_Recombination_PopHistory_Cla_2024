#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --time=00:00:10
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=bedtools_intersect_rename
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

# Load modules
module purge
module load bedtools/2.31.0 

OUTDIR=/projects/ag-waldvogel/pophistory/CRIP/comparison_Cla_rho/02-bedtools_intersect

BEDA=/projects/ag-waldvogel/pophistory/CRIP/MELT-analysis/Cla1/AllPop/2-GroupAnalysis

BEDB=/projects/ag-waldvogel/pophistory/CRIP/MELT-analysis/Cla1/AllPop/1-IndivAnalysis

IND=("MF1" "MF2" "MF3" "MF4" \
       "MG2" "MG3" "MG4" "MG5" \
       "NMF1" "NMF2" "NMF3" "NMF4" \
       "SI1" "SI2" "SI3" "SI4" \
       "SS1" "SS2" "SS3" "SS4")

CHR=("Chr1" "Chr2" "Chr3" "Chr4")

POP=("MF" "MG" "NMF" "SI" "SS")

# CLA: USE IndivAnalysis of AllPop Run

# RENAME INSERTIONS BASED ON GROUP ANALYSIS
for ind in "${IND[@]}"; do
    cd ${OUTDIR}
    # Identify overlapping regions
    bedtools intersect -a ${BEDA}/Cla1.master.bed -b ${BEDB}/"${ind}".bwamem.sambl.sort.Cla1.tmp.bed -wa -wb > "${ind}"_overlapping_regions.bed   
    
    # Rename overlapping regions
    awk 'BEGIN {OFS="\t"} {if ($4 != $10) print $7, $8, $9, $4, $11, $12; else print $7, $8, $9, $10, $11, $12}' "${ind}"_overlapping_regions.bed > "${ind}".bwamem.sambl.sort.Cla1.tmp_renamed.bed

    # Merge renamed regions with non-overlapping regions
    cat "${ind}".bwamem.sambl.sort.Cla1.tmp_renamed.bed  <(bedtools intersect -v -a ${BEDB}/"${ind}".bwamem.sambl.sort.Cla1.tmp.bed -b ${BEDA}/Cla1.master.bed) > "${ind}".bwamem.sambl.sort.Cla1.tmp_renamed.final.bed
    
    # Get a list of non-overlapping regions (should be empty)
    bedtools intersect -v -a ${BEDB}/"${ind}".bwamem.sambl.sort.Cla1.tmp.bed -b ${BEDA}/Cla1.master.bed > "${ind}"_non_overlapping_regions.bed

done

wait

echo "Rename insertions done"

# SPLIT CLA-BED BY CHROMOSOME AND SORT
for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do
    cd ${OUTDIR}
    cat "${ind}".bwamem.sambl.sort.Cla1.tmp_renamed.final.bed | grep "${chr}" | sort -k1,1 -k2,2n > "${ind}"_"${chr}".bwamem.sambl.sort.Cla1.tmp_renamed.final.sorted.bed       
  done
done

wait

echo "Split Cla-bed by chromosome and sort it done"

wait
