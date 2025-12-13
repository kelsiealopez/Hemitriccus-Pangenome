#!/bin/bash
#SBATCH -p test
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH -o blast_CNV_6433_only_absent_%j.out
#SBATCH -e blast_CNV_6433_only_absent.sh_%j.err 
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
