#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=80GB
#SBATCH --time=20:00:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=preprocessing_MELT
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

DIR=/projects/ag-waldvogel/pophistory/CRIP/bam-files/samblaster
BIN=/home/lpettric/bin

module purge
module load samtools/1.13

cd ${DIR}

# Add MQ tags, convert to bam and sort bam
ls -1 *_bwamem.sam | sed 's/_bwamem.sam//g' > list-crip

while read f; do ${BIN}/samblaster/samblaster --addMateTags -i $f"_bwamem.sam" -o $f"_bwamem.sambl.sam"; done < list-crip |& tee samblaster.log 
 
while read f; do samtools view -b $f"_bwamem.sambl.sam" > $f".bwamem.sambl.bam"; done < list-crip

cat list-crip | ${BIN}/parallel/bin/parallel -j 20 'samtools sort -@ 4 -o {}.bwamem.sambl.sort.bam {}.bwamem.sambl.bam'

wait

echo "MQ tag added and conversion to bam done"

wait

# Gather flagstat statistics

${BIN}/parallel/bin/parallel -j 20 'samtools flagstat {}.bwamem.sambl.sort.bam > {}.bwamem.sambl.sort.flagstat' < list-crip

wait

echo "Flagstat staistics gathered"

wait

# Remove low quality
# remove reads with map. qual. smaller 30
# keep read mapped in proper pair (-f 0x0002)
# remove unmapped read (-F 0x0004) and mate unmapped (-F 0x0008) 
${BIN}/parallel/bin/parallel -j 20 'samtools view -q 30 -f 0x0002 -F 0x0004 -F 0x0008 -b -o {}.bwamem.sambl.sort.q30.bam {}.bwamem.sambl.sort.bam' < list-crip

wait

echo "Low quality alignments removed"

wait

# Remove duplicates
module purge
module load miniconda/py38_4.12.0

conda activate /scratch/lpettric/conda/picard_env

while read f; do picard MarkDuplicates -I $f".bwamem.sambl.sort.q30.bam" -O $f".bwamem.sambl.sort.q30.rmd.bam" -M $f".bwamem.sambl.sort.q30.rmd.stat" -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true; done < list-crip

wait

echo "Duplicates removed"

wait

# Collect mapping statistics
module purge
module load miniconda/py38_4.12.0
conda activate /scratch/lpettric/conda/qualimap_env

# Filtered bam-files
qualimap multi-bamqc -d ./list-qualimap -r --java-mem-size=4G -outformat PDF:HTML -outdir ./qualimap-output

# Unfiltered bam-files
qualimap multi-bamqc -d ./list-qualimap_unfiltered -r --java-mem-size=4G -outformat PDF:HTML -outdir ./qualimap-output_unfiltered

wait

echo "Mapping statistics generated"

wait
