### Fetch Databases

Fetch the NCBI Taxdump

    mkdir -p taxdump; cd taxdump; curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -; cd ..;

Fetch the nt database

    mkdir -p nt wget
    "<ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt>.??.tar.gz" -P nt/ &&\
    for file in nt/\*.tar.gz;\
    do tar xf \$file -C nt && rm \$file;\
    done

Fetch any BUSCO lineages that you plan to use

### BUSCO

    # Download BUSCO for offline usage
    wget -q -O diptera_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/diptera_odb10.2020-08-05.tar.gz" \
           && tar xf diptera_odb10.gz 

    singularity exec -B /scratch/lpettric/blobtools/CRIP/:/busco_wd/ -B /scratch/lpettric/busco_downloads/:/busco/ /scratch/lpettric/singularity/busco_v5.3.2_cv1.sif busco -m genome -i /busco_wd/Chironomus_riparius_genome_010921.fasta -o /busco_wd/CRIP/Chironomus_riparius_genome_010921.busco.diptera_odb10.tsv -l diptera_odb10 -f --offline --download_path /busco/

### Blastx

    # Cheops 1
    module purge
    module load miniconda/py38_4.9.2

    eval "$(conda shell.bash hook)"

    conda activate /scratch/lpettric/conda/blastx_env/

    blastn -db /scratch/lpettric/nt/nt \
           -query /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta \
           -outfmt "6 qseqid staxids bitscore std" \
           -max_target_seqs 10 \
           -max_hsps 1 \
           -evalue 1e-25 \
           -num_threads 32 \
           -out /scratch/lpettric/blobtools/CRIP/Chironomus_riparius_genome_010921.ncbi.blastn.run.out

### Blobtools


    # taxid
    --taxid 315576

    # create dataset
    singularity exec -B /scratch/lpettric/blobtools/CRIP/:/schluppsi -B /scratch/lpettric/blobtools/CRIP/output/:/outolino -B /scratch/lpettric/taxdump/:/taxdumpi /opt/rrzk/software/singularity_images/blobtoolkit_latest.sif blobtools create \
        --fasta /schluppsi/Chironomus_riparius_genome_010921.fasta \
        --meta /schluppsi/Chironomus_riparius_genome_010921.yaml \
        --hits /schluppsi/Chironomus_riparius_genome_010921.ncbi.blastn.run.out \
        --busco /schluppsi/Chironomus_riparius_genome_010921.busco.diptera_odb10.tsv/run_diptera_odb10/full_table.tsv \
        --taxid 315576 \
        --taxdump /taxdumpi \
        --taxrule bestsumorder \
        --cov /schluppsi/MF1.bwamem.sort.q30.rmd.bam \
        /outolino/


    # metadata
    nano Chironomus_riparius_genome_010921.yaml

    assembly:
      alias: Chironomus_riparius_genome_010921
      record_type: chromosome
    taxon:
      name: Chironomus riparius


    # View blobplot
    singularity exec -B /scratch/lpettric/blobtools/CRIP/:/schluppsi -B /scratch/lpettric/blobtools/CRIP/output/output-view:/outolino /opt/rrzk/software/singularity_images/blobtoolkit_latest.sif blobtools view --remote --out /outolino/output --view blob /schluppsi/output
