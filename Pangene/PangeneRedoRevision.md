
## Used this R script to get a tsv file of genes missing in individuals (cnv_artifact_filtered_genes_absent_by_individual.tsv)

## Run pangeneFilteringtoGetPAVList.R


## make text files like this for each individual
```bash
[kelsielopez@boslogin08 blast_validate_cnv]$ head 6433_absent_Dec.txt
GOLPH3L
galGal.LOC427618
PCDHA13
NRN1
SULT1B1
FAM134B
OR9Q1
RACGAP1L
EDPE
chiLan.LOC116781923
```


## Run blast_CNV_6433_only_absent_test_partition_dec.sh for each individual. Probably easier to make this an array for all assemblies 

```bash
[kelsielopez@boslogin08 blast_validate_cnv]$ cat blast_CNV_6433_only_absent_test_partition_dec.sh
#!/bin/bash
#SBATCH -p test
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH -o blast_CNV_6433_only_absent_test_partition_dec_%j.out
#SBATCH -e blast_CNV_6433_only_absent_test_partition_dec_%j.err 
#SBATCH --mem=150G
#SBATCH --mail-type=END

# 6433


GENELIST="6433_absent_Dec.txt"
FASTAIN="/n/netscratch/edwards_lab/Lab/kelsielopez/blast_validate_cnv/protein_gene_colon_transcript.faa"
OUTFA="6433_absent_Dec.faa"

# make a protein fasta of the genes that are missing in that individual from the protein fasta that I used for pangene 

bioawk -c fastx -v genes="$GENELIST" '
    BEGIN{
        while ((getline < genes) > 0) {
            gene[$1]=1
        }
    }
    {
        split($name,a,":")
        if (gene[a[1]]) {print ">"$name"\n"$seq}
    }' $FASTAIN > $OUTFA


reads="/n/holylfs04/LABS/edwards_lab/Lab/klopez/Hemitriccus/6433/00_raw_reads/m64408e_220825_141218.hifi_reads.fq.gz"
reads_fa="/n/holylfs04/LABS/edwards_lab/Lab/klopez/Hemitriccus/6433/00_raw_reads/m64408e_220825_141218.hifi_reads.fasta"
blast_dir="/n/netscratch/edwards_lab/Lab/kelsielopez/blast_validate_cnv"
sample="6433"


cd ${blast_dir}

# turn fastq file into a fasta file 
zcat ${reads} | seqtk seq -A - > ${reads_fa}

# make that the blast database 
makeblastdb -in ${reads_fa} -dbtype nucl -out ${blast_dir}/${sample}_DB

# blast the missing genes against the reads 

GENE_FA="/n/netscratch/edwards_lab/Lab/kelsielopez/blast_validate_cnv/6433_absent_Dec.faa"

for BIRDDB in ${sample}_DB
do
    tblastn -query $GENE_FA -db $BIRDDB -out ${BIRDDB}_vs_genes.blast.DEC -evalue 1e-5 -outfmt 6 -num_threads 32
done

```


## testing e value e-10 and sequence identity 95% 

```bash
cd /n/netscratch/edwards_lab/Lab/kelsielopez/blast_validate_cnv

# start file with header
#printf "Gene\tIndividual\n" > blast_validated_hits.tsv


# function: add hits from one BLAST file
add_hits () {
    local DECFILE="$1"
    local IND="$2"    # individual name to match R 'individuals'

    awk -v IND="$IND" '
        $11 <= 1e-10 && $3 >= 95 {
            # $1 is like chiLan.LOC116784105:rna-XM_032682289.1.2
            split($1, a, ":");
            gene = a[1];
            print gene "\t" IND
        }
    ' "$DECFILE" | sort -u >> blast_validated_hits.tsv
}

# 6371
add_hits 6371_DB_vs_genes.blast.DEC HMRG_6371

# 6386
add_hits 6386_DB_vs_genes.blast.DEC HMRG_6386

# 6388
add_hits 6388_DB_vs_genes.blast.DEC HMRG_6388

# 6431
add_hits 6431_DB_vs_genes.blast.DEC HMRG_6431

# 6433
add_hits 6433_DB_vs_genes.blast.DEC HMRG_6433

```


## then re run this R script (FinalPangene.R) to get final counts 
