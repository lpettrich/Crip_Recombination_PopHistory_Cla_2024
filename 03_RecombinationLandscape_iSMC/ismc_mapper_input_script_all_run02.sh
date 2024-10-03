#!/bin/bash

# --------------------------------------------------------------------------------------
#                                  iSMC Script - Mapper
# --------------------------------------------------------------------------------------

# Change input information (MF1 already done)
IND=("MF1" "MF2" "MF3" "MF4" "MG2" "MG3" "MG4" "MG5" "NMF1" "NMF2" "NMF3" "NMF4" "SI1" "SI2" "SI3" "SI4" "SS1" "SS2" "SS3" "SS4")
CHR=("Chr1" "Chr2" "Chr3" "Chr4")

# Change directory and create ismc mapper file
for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do    
    cd /home/lpettrich/pophistory/ismc/"${ind}"/"${chr}"
    cat > "${ind}"_"${chr}"_ismc_mapper.bpp <<EOF 
    dataset_label = ${ind}_${chr}_ismc
    bin_sizes = (1000, 100000)
    tab_file_path = ${ind}_${chr}_phased_merged.vcf.gz_ismc.tab
    bin_rate = rho 
    tmrca = true 
EOF
  done
done

