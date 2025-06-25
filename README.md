# Hemitriccus-Pangenome
Code for hemitriccus pangenome project 

## Genomic Feature - SV analysis 

**Get density of genes:**

```bash
############################################################################
############################# Gene Density Pipeline ########################
############################################################################

# Path to the GFF3 annotation file
gff="/n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/hemMar_with_CDS_corrected.gff3"

# Extract gene coordinates from GFF3, convert to BED format (0-based), and save to genes.bed
awk '$3=="gene" {OFS="\t"; print $1, $4-1, $5}' ${gff} > genes.bed

# ----------------------------------------------------------------------------
# Section: Getting Genome/Chromosome Sizes from FASTA
# ----------------------------------------------------------------------------

# (Optional) Index your FASTA file to prepare for downstream steps
samtools faidx your.fasta

# Extract sequence names and lengths from the FASTA index; write to genome.sizes
cut -f1,2 your.fasta.fai > genome.sizes

# ----------------------------------------------------------------------------
# Section: Useful One-Offs
# ----------------------------------------------------------------------------

# Count the number of unique entries (likely gene IDs) in first column of a tab-delimited file
awk '{print $1}' HemMar_10_haps_June_2025_new_faa_fix_PAV.Rtab | sort -u | wc -l

# ----------------------------------------------------------------------------
# Section: Preparing Final Inputs for Windowing & Density Calculation
# ----------------------------------------------------------------------------

# Define path to reference FASTA and its .fai index
REF="/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/hemMar_complete_sorted_prefixed_pggb_subset_JBAT.FINAL.full.soft.mask.fasta"
REF_fai="/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/hemMar_complete_sorted_prefixed_pggb_subset_JBAT.FINAL.full.soft.mask.fasta.fai"

# Index the reference FASTA (creates .fai file if not present)
samtools faidx ${REF}

# Extract reference sequence names and sizes to create the genome size file for bedtools
cut -f1,2 ${REF_fai} > genome.sizes

# Remove 'HemMar#1#' prefix from sequence names in genome.sizes (if present)
sed -i 's/HemMar#1#//g' genome.sizes

# ----------------------------------------------------------------------------
# Section: Windowing and Gene Density Calculation
# ----------------------------------------------------------------------------

# Path to your bedtools binary
bedtools_path="/n/home03/kelsielopez/bedtools2/bin/bedtools"

# Create non-overlapping 10kb windows across each sequence/chromosome
${bedtools_path} makewindows -g genome.sizes -w 10000 > windows.bed

# For each 10kb window, count the number of genes that overlap; write output to gene_density_10kb.bed
${bedtools_path} intersect -a windows.bed -b genes.bed -c > gene_density_10kb.bed
```
