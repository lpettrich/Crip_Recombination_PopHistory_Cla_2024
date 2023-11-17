#!/bin/bash

# --------------------------------------------------------------------------------------
#                                  iSMC Script - Mapper RUN
# --------------------------------------------------------------------------------------

IND=("MF2" "MF3" "MF4" "MG2" "MG3" "MG4" "MG5" "NMF1" "NMF2" "NMF3" "NMF4" "SI1" "SI2" "SI3" "SI4" "SS1" "SS2" "SS3" "SS4")
CHR=("Chr1" "Chr2" "Chr3" "Chr4")

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do    
    cd /home/lpettrich/pophistory/ismc/"${ind}"/"${chr}"
    /home/lpettrich/pophistory/ismc/ismc_mapper params="${ind}"_"${chr}"_ismc_mapper.bpp |& tee "${ind}"_"${chr}"_ismc_mapper.log
  done
done
