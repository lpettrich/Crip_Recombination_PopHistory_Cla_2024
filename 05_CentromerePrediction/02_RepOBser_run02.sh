#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128000M
#SBATCH --cpus-per-task=15
#SBATCH --time=05:00:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=RepOBserver_run02
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

# Load modules
module purge
module load miniconda/py311_23.9.0
conda activate /scratch/lpettric/conda/RepOBserver_env/

# Set directory
DIR=/projects/ag-waldvogel/pophistory/CRIP/centromeres/run02

cd ${DIR}

# Start RepeatOBserverV1
srun Setup_Run_Repeats.sh -i Criparius -f Chironomus_riparius_genome_010921.fasta -h H0 -c 15 -m 128000 -g TRUE
