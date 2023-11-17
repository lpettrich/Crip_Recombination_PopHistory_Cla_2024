#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=80GB
#SBATCH --time=8:00:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=MELT-for-MF
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

module purge
module load bowtie2/2.4.1
module load samtools/1.13
# cheops1 runs with java1.8_345 without loading module but for some reason it's not working
module load openjdk/1.8.0_202
java -version 

cd /projects/ag-waldvogel/pophistory/CRIP/MELT-analysis

# Melt Analysis 
### with pophistory samples mapped on Crip4.0
### melt-run for Cla1
## using script environment of Malte Petersen: 

# 1) IndivAnalysis

while read list; do bash /projects/ag-waldvogel/pophistory/CRIP/MELT-analysis/MELT-IndivAnalysis.sh parameterfile_MF_Cla1.tsv $list Cla1/MF ; done < list_MF

wait

echo "IndivAnalysis MF Cla1 done"

wait

# 2) GroupAnalysis

bash /projects/ag-waldvogel/pophistory/CRIP/MELT-analysis/MELT-GroupAnalysis.sh parameterfile_MF_Cla1.tsv Cla1/MF 

wait

echo "GroupAnalysis MF Cla1 done"

wait

# 3) Genotyping

while read list; do bash /projects/ag-waldvogel/pophistory/CRIP/MELT-analysis/MELT-Genotype.sh parameterfile_MF_Cla1.tsv $list Cla1/MF ; done < list_MF

wait

echo "Genotyping MF Cla1 done"

wait

# 4) MakeVCF

bash /projects/ag-waldvogel/pophistory/CRIP/MELT-analysis/MELT-MakeVCF.sh parameterfile_MF_Cla1.tsv Cla1/MF 

wait

echo "MakeVCF MF Cla1 done"

wait


