#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80GB
#SBATCH --time=00:10:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=bedtools_closest
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

module purge
module load bedtools/2.31.0 
module load R/4.3.1_system

export R_LIBS_USER=/home/lpettric/R/x86_64-pc-linux-gnu-library/4.3

OUTDIR=/projects/ag-waldvogel/pophistory/CRIP/comparison_Cla_rho/01-bedtools_closest

CLADIR=/projects/ag-waldvogel/pophistory/CRIP/MELT-analysis/Cla1

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

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do    
    cd ${RHODIR}/"${ind}"/"${chr}"
    cat "${ind}"_"${chr}"_ismc.rho.10kb.bedgraph | grep -v '^chrom' \
    | sed "s/chr1/${chr}/g" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.10kb.bedgraph_renamed.bed
    sort -k1,1 -k2,2n ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.10kb.bedgraph_renamed.bed > ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed
  done
done

wait

echo "Converting bedgraph to bed done"

wait

# GET MEAN RHO/10kb FOR EACH POP
### bedtools merge: -c=column, -o=operator
### All bed-files are moved to shared directory

cd ${RHODIR}/bed-files

# Merge for MF
for chr in "${CHR[@]}"; do
    cat ${IND[0]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[1]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    ${IND[2]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[3]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    | sort -k1,1 -k2,2n > ${POP[0]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.bed
done

wait

echo "Merge Rho/10kb for ${IND[0]} ${IND[1]} ${IND[2]} ${IND[3]} of ${POP[0]}"

# Merge for MG
for chr in "${CHR[@]}"; do
    cat ${IND[4]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[5]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    ${IND[6]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[7]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    | sort -k1,1 -k2,2n > ${POP[1]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.bed
done

wait

echo "Merge Rho/10kb for ${IND[4]} ${IND[5]} ${IND[6]} ${IND[7]} of ${POP[1]}"

# Merge for NMF
for chr in "${CHR[@]}"; do
    cat ${IND[8]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[9]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    ${IND[10]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[11]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    | sort -k1,1 -k2,2n > ${POP[2]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.bed
done

wait

echo "Merge Rho/10kb for ${IND[8]} ${IND[9]} ${IND[10]} ${IND[11]} of ${POP[2]}"

# Mean for SI
for chr in "${CHR[@]}"; do
    cat ${IND[12]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[13]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    ${IND[14]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[15]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    | sort -k1,1 -k2,2n > ${POP[3]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.bed
done

wait

echo "Merge Rho/10kb for ${IND[12]} ${IND[13]} ${IND[14]} ${IND[15]} of ${POP[3]}"

# Merge for SS
for chr in "${CHR[@]}"; do
    cat ${IND[16]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[17]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    ${IND[18]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed ${IND[19]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.bed \
    | sort -k1,1 -k2,2n > ${POP[4]}_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.bed
done

wait

echo "Merge Rho/10kb for ${IND[16]} ${IND[17]} ${IND[18]} ${IND[19]} of ${POP[4]}"

wait

# Use R to get mean

R --vanilla -f Merged-bed-files-mean.R

wait

echo "Callaculating means done"
# Saved as "${pop}"_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.mean.bed

wait

# CLA: REMOVE SITES WITH HOM REF ALLELES (0/0)[ac0] AND NO-CALL (./.)[s25] 
# CONVERT VCF TO BED AND SORT BED
# SPLIT BED BY CHROMOSOME

cd ${CLADIR}

for pop in "${POP[@]}"; do
    cd ${CLADIR}/"${pop}"/4-MakeVCF
    cat Cla1.final_comp.vcf | grep -v 'ac0' > ${CLADIR}/bed-files/"${pop}"_Cla1.final_comp_noac0.vcf
    cd ${CLADIR}/bed-files
    cat "${pop}"_Cla1.final_comp_noac0.vcf | grep -v 's25' > "${pop}"_Cla1.final_comp_noac0_nos25.vcf
done

wait 

# https://www.biostars.org/p/106249/
module load miniconda/py38_4.12.0
conda activate bedops_env

for pop in "${POP[@]}"; do
    cd ${CLADIR}/bed-files
    vcf2bed --do-not-sort < "${pop}"_Cla1.final_comp_noac0_nos25.vcf > "${pop}"_Cla1.final_comp_noac0_nos25.bed
    sort -k1,1 -k2,2n "${pop}"_Cla1.final_comp_noac0_nos25.bed > "${pop}"_Cla1.final_comp_noac0_nos25.sorted.bed
done

wait

conda deactivate

echo "Filtering vcf, converting to bed and sorting done"

wait

# SPLIT CLA-BED BY CHROMOSOME
for pop in "${POP[@]}"; do
  for chr in "${CHR[@]}"; do
    cd ${CLADIR}/bed-files
    cat "${pop}"_Cla1.final_comp_noac0_nos25.sorted.bed | grep "${chr}" | sort -k1,1 -k2,2n > "${pop}"_"${chr}"_Cla1.final_comp_noac0_nos25.sorted.bed      
  done
done

wait

echo "Split Cla-bed by chromosome done"

wait

# UPDATE THIRD COLUMN TO REFLECT CLUSTER SIZE
for pop in "${POP[@]}"; do
  for chr in "${CHR[@]}"; do
    cd ${CLADIR}/bed-files
    awk 'BEGIN{OFS=FS="\t"} {if (match($9, /SVLEN=([0-9]+)/)) {len = substr($9, RSTART+6, RLENGTH-6); $3 = $2 + len} print}' \
    "${pop}"_"${chr}"_Cla1.final_comp_noac0_nos25.sorted.bed > "${pop}"_"${chr}"_Cla1.final_comp_noac0_nos25.sorted_clustersize.bed
  done
done

wait

echo "Replacing third column with actual cluster size done"

wait


# BEDTOOLS CLOSEST
### with pophistory samples mapped on Crip4.0
### for each individual and each chromosome 
### Cla regions compared to rho (10 kb)
cd ${OUTDIR}

for pop in "${POP[@]}"; do
  for chr in "${CHR[@]}"; do
    bedtools closest -d -t first -a ${RHODIR}/bed-files/"${pop}"_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.mean.bed \
    -b ${CLADIR}/bed-files/"${pop}"_"${chr}"_Cla1.final_comp_noac0_nos25.sorted_clustersize.bed > \
    "${pop}"_"${chr}"_ismc.rho.10kb.bedgraph_renamed.sorted.merged.mean_distToNearestCla1_clustersize.bed
  done
done

wait

echo "Bedtools closest finished"

wait
