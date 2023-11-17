#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=50GB
#SBATCH --time=20:00:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=repeatmaskerCRIP
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

module purge
module load repeatmasker/4.1.1

LIB=/projects/ag-waldvogel/pophistory/CRIP/repeatmasker/RM-only-with-Cla1-no-Cla-element/TE-lib-CRIP-Vladimir2_MolEcol17_RMhead_merged_rmd_noClasome_noCla-element.fasta
DIR=/projects/ag-waldvogel/pophistory/CRIP/repeatmasker/RM-only-with-Cla1-no-Cla-element
GENOME=/projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta

cd ${DIR}

RepeatMasker -s -xsmall -cutoff 250 -u -gff -pa 10 -lib ${LIB} -dir ${DIR} -engine rmblast ${GENOME}
