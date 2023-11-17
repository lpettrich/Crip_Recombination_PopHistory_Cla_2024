## 1. First check with fastqc and multiqc

FastQC

    module purge
    module load openjdk/1.8.0_202

     cd /projects/ag-waldvogel/pophistory/CRIP/

    /home/lpettric/bin/FastQC/fastqc --threads 10 -o ./MF1/ ./MF1/MF1_R1.paired.fastq_true.gz ./MF1/MF1_R2.paired.fastq_true.gz &
    wait
    ...

MultiQC

    module purge
    module load miniconda/py38_4.9.2
    eval "$(conda shell.bash hook)"
    conda activate /scratch/lpettric/conda/multiqc_env/

    cd /projects/ag-waldvogel/pophistory/CRIP/

    multiqc --filename multiqc_report_all_trimmed.html -o /projects/ag-waldvogel/pophistory/CRIP/ --file-list /projects/ag-waldvogel/pophistory/CRIP/multiqc-file.list

## 2. Mapping with bwa mem

    module purge

    cd /projects/ag-waldvogel/pophistory/CRIP/

    /home/lpettric/bin/bwa/bwa mem -t 20 -M -R '@RG\tID:MF\tSM:MF1\tPL:ILLUMINA' /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta ./MF1/MF1_R1.paired.fastq_true.gz ./MF1/MF1_R2.paired.fastq_true.gz > ./bam-files/MF1_bwamem.sam &
    wait
    ...

## 3. Sort bam files

    module purge
    module load samtools

    cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

    ls -1 *_bwamem.sam | sed 's/_bwamem.sam//g' > list-crip
    while read f; do samtools view -b $f"_bwamem.sam" > $f".bam" ;done < list-crip
    cat list-crip | /home/lpettric/bin/parallel/bin/parallel -j 20 'samtools sort -@ 4 -o {}.bwamem.sort.bam {}.bam'

"list-crip" contains names of files without file-extension

## 4. Collect flagstat statistics

    module purge
    module load samtools

    cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

    /home/lpettric/bin/parallel/bin/parallel -j 20 'samtools flagstat {}.bwamem.sort.bam > {}.bwamem.sort.flagstat' < list-crip

## 5. Remove low qulaity alignments

    module purge
    module load samtools

    cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

    /home/lpettric/bin/parallel/bin/parallel -j 20 'samtools view -q 30 -f 0x0002 -F 0x0004 -F 0x0008 -b -o {}.bwamem.sort.q30.bam {}.bwamem.sort.bam' < list-crip

## 6. Remove duplicates

    module purge
    module load miniconda/py38_4.9.2
    eval "$(conda shell.bash hook)"
    conda activate /scratch/lpettric/conda/picard_env

    cd /projects/ag-waldvogel/pophistory/CRIP/bam-files/

    picard MarkDuplicates -I MF1.bwamem.sort.q30.bam -O MF1.bwamem.sort.q30.rmd.bam -M MF1.bwamem.sort.q30.rmd.stat -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true &
    wait
    ...

## 7. Collect Mapping statistics

    module purge
    module load miniconda/py38_4.9.2
    eval "$(conda shell.bash hook)"
    conda activate /scratch/lpettric/conda/qualimap_env

    cd /projects/ag-waldvogel/pophistory/CRIP/bam-files/

    qualimap multi-bamqc -d ./list-qualimap -r --java-mem-size=4G -outformat PDF:HTML -outdir ./qualimap-output
