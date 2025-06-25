# Hemitriccus-Pangenome
Code for hemitriccus pangenome project 

## Genomic Feature - SV analysis 

**Get density of genes:**

```bash
############################################################################
############################# Gene Density  ################################
############################################################################

# Path to the GFF3 annotation file
gff="/n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/hemMar_with_CDS_corrected.gff3"

# Extract gene coordinates from GFF3, convert to BED format (0-based), and save to genes.bed
awk '$3=="gene" {OFS="\t"; print $1, $4-1, $5}' ${gff} > genes.bed

# ----------------------------------------------------------------------------
# Getting Genome/Chromosome Sizes from FASTA
# ----------------------------------------------------------------------------

# (Optional) Index your FASTA file to prepare for downstream steps
samtools faidx your.fasta

# Extract sequence names and lengths from the FASTA index; write to genome.sizes
cut -f1,2 your.fasta.fai > genome.sizes

# ----------------------------------------------------------------------------
# Preparing Final Inputs for Windowing & Density Calculation
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
# Windowing and Gene Density Calculation
# ----------------------------------------------------------------------------

# Path to your bedtools binary
bedtools_path="/n/home03/kelsielopez/bedtools2/bin/bedtools"

# Create non-overlapping 10kb windows across each sequence/chromosome
${bedtools_path} makewindows -g genome.sizes -w 10000 > windows.bed

# For each 10kb window, count the number of genes that overlap; write output to gene_density_10kb.bed
${bedtools_path} intersect -a windows.bed -b genes.bed -c > gene_density_10kb.bed

# Count the number of 'genes' in each 200kb genomic window,
# writing the output to gene_density_200kb.raw.bed.
# The fourth column of the output will be the count of overlapping genes per window.
bedtools intersect -a genome.200kb.windows.bed -b genes.bed -c > gene_density_200kb.raw.bed

# (Same as above, using an explicit bedtools path variable)
# Count the number of genes in each 200kb window.
${bedtools_path} intersect -a genome.200kb.windows.bed -b genes.bed -c > gene_density_200kb.raw.bed

# Calculate gene *density* (number of genes per base) for each 200kb window.
# Here, $4 = gene count; divide by 200,000 (window size) to get density.
# Output columns: chrom, start, end, density.
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/200000}' gene_density_200kb.raw.bed > gene_density_200kb.bed

# Repeat process for 10kb windows: count genes in each window.
${bedtools_path} intersect -a genome.10kb.windows.bed -b genes.bed -c > gene_density_10kb.raw.bed

# Calculate gene density (per base) for each 10kb window.
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/10000}' gene_density_10kb.raw.bed > gene_density_10kb.bed
```


**Get average recombionation rate in windows:**

```bash
############################################################################
############################################################################
############################# Recombination ################################
############################################################################
############################################################################

# Below are the full file paths for each ReLERNN output you want to merge
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/example_output/pggb_cleaned_final_biallelic_snp.PREDICT.BSCORRECTED.txt
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_1_first_half/example_output/HemMar_scaffold_1_first_half_nomissing.PREDICT.BSCORRECTED.txt
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_1_second_half/example_output/HemMar_scaffold_1_second_half_nomissing.PREDICT.BSCORRECTED.txt
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_2_first_half/example_output/HemMar_scaffold_2_first_half_nomissing.PREDICT.BSCORRECTED.txt
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_2_second_half/example_output/HemMar_scaffold_2_second_half_nomissing.PREDICT.BSCORRECTED.txt

# -----------------------------------------------
# Count lines in each ReLERNN result file (sanity check, e.g. to confirm no data loss/duplication)
wc -l /n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/example_output/pggb_cleaned_final_biallelic_snp.PREDICT.BSCORRECTED.txt
wc -l /n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_1_first_half/example_output/HemMar_scaffold_1_first_half_nomissing.PREDICT.BSCORRECTED.txt
wc -l /n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_1_second_half/example_output/HemMar_scaffold_1_second_half_nomissing.PREDICT.BSCORRECTED.txt
wc -l /n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_2_first_half/example_output/HemMar_scaffold_2_first_half_nomissing.PREDICT.BSCORRECTED.txt
wc -l /n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_2_second_half/example_output/HemMar_scaffold_2_second_half_nomissing.PREDICT.BSCORRECTED.txt

# -----------------------------------------------
# Change to the working directory for combining results
cd /n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN

# Create a new combined file and add a header line with column names
echo -e "chrom\tstart\tend\tnSites\trecombRate\tCI95LO\tCI95HI" > all_chroms_combined_no_missing.txt

# Combine all ReLERNN results files into one, skipping header of each,
# and cleaning up any strange b'' Python string-formatting artifacts in column 1
awk 'FNR>1 {gsub(/^b'\''|'\''$/, "", $1); print}' \
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/example_output/pggb_cleaned_final_biallelic_snp.PREDICT.BSCORRECTED.txt \
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_1_first_half/example_output/HemMar_scaffold_1_first_half_nomissing.PREDICT.BSCORRECTED.txt \
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_1_second_half/example_output/HemMar_scaffold_1_second_half_nomissing.PREDICT.BSCORRECTED.txt \
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_2_first_half/example_output/HemMar_scaffold_2_first_half_nomissing.PREDICT.BSCORRECTED.txt \
/n/netscratch/edwards_lab/Lab/kelsielopez/ReLERNN/HemMar/scaff_1_2_no_missing/HemMar_scaffold_2_second_half/example_output/HemMar_scaffold_2_second_half_nomissing.PREDICT.BSCORRECTED.txt \
>> all_chroms_combined_no_missing.txt

# -----------------------------------------------
# Make a backup or intermediate copy for further editing
cp all_chroms_combined_no_missing.txt all_chroms_combined_no_missing_noHemMar.txt

# Remove 'HemMar#1#' (prefix) from chromosome names in the combined file (in-place edit)
sed -i 's/HemMar#1#//g' all_chroms_combined_no_missing_noHemMar.txt

# -----------------------------------------------
# Get the number of unique chromosome names in combined file
awk '{print $1}' all_chroms_combined_no_missing_noHemMar.txt | sort -u | wc -l

# -----------------------------------------------
# Make a BED file of recombination rates (chrom, start, end, recombRate)
awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $5}' all_chroms_combined_no_missing_noHemMar.txt > all_chroms_combined_no_missing_noHemMar.recombRate.bed

# Make a BED file of recombination rates (chrom, start, end, recombRate, nSites)
# Includes the number of SNP sites used in each window as the 5th column
awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $4}' all_chroms_combined_no_missing_noHemMar.txt > all_chroms_combined_no_missing_noHemMar.recombRate.nSites.bed


nano average_recomb_rate.py
# below is the python script

import pandas as pd

# Setup
columns = ['chrom','start','end','recombRate','nSites']
df = pd.read_csv('all_chroms_combined_no_missing_noHemMar.recombRate.nSites.bed', sep='\t', header=None, names=columns)

bin_size = 200_000

# Ensure positions are int
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)

# Assign a bin index per row
df['bin_index'] = df['start'] // bin_size

# Group by chrom and bin index, aggregate mean recombination
result = (
    df.groupby(['chrom', 'bin_index'])
      .agg(
          start  = ('start',  lambda x: x.min() // bin_size * bin_size),
          end    = ('start',  lambda x: (x.min() // bin_size + 1) * bin_size),
          recombRate = ('recombRate','mean')
      )
      .reset_index()
)

# Save output
result[['chrom','start','end','recombRate']].to_csv('recomb_rate_200kb_windows.bed', sep='\t', header=False, index=False)
```


**Take average GC conetent in windows:**


```bash
############################################################################
############################################################################
############################# GC CONTENT ###################################
############################################################################
############################################################################
############################################################################

REF="/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/hemMar_complete_sorted_prefixed_pggb_subset_JBAT.FINAL.full.soft.mask.fasta"

#if making a bed file of GC content...

cd /n/netscratch/edwards_lab/Lab/kelsielopez/feature_importance


nano extract_GC_content.py

# below is the python code 

from Bio import SeqIO
import pandas as pd

fasta_path = "/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/hemMar_complete_sorted_prefixed_pggb_subset_JBAT.FINAL.full.soft.mask.fasta"
bed_out = "gc_content.10kb.bed"

window_size = 10000  # <-- Change this as needed!

out = []
for record in SeqIO.parse(fasta_path, "fasta"):
    seq = str(record.seq).upper()
    seqlen = len(seq)
    for start in range(0, seqlen, window_size):
        end = min(start + window_size, seqlen)
        window_seq = seq[start:end]
        gc = window_seq.count('G') + window_seq.count('C')
        atgc = gc + window_seq.count('A') + window_seq.count('T')
        gc_content = gc / atgc if atgc > 0 else 0
        out.append([record.id, start, end, round(gc_content, 4)])

# Save to BED
df = pd.DataFrame(out, columns=['chrom','start','end','gc_content'])
df.to_csv(bed_out, sep='\t', index=False, header=False)


# copy so i don't overwrite the oroginal file 
cp gc_content.10kb.bed gc_content.10kb.no.hemMar.bed

# also remove hemMar prefix
sed -i 's/HemMar#1#//g' gc_content.10kb.no.hemMar.bed

```
