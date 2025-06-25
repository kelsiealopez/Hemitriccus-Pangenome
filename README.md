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


**Take proportions of repeat content in windows:**


```bash

############################################################################
############################# Repeat Analysis ###############################
############################################################################

# ------------------------
# 1. Define broad repeat class mapping (in Python, see below)
# ------------------------

# Mapping of various repeat subtypes (from annotation) into broad repeat classes
# (e.g. LINE, SINE, DNA, LTR, etc.); will be used in a Python script below

# ------------------------
# 2. Explore Repeat BED File
# ------------------------

# Count the unique repeat classes in your reformatted repeat BED file
awk '{print $4}' repeat_overlaps_reformatted_final.bed | sort -u | wc -l

# ------------------------
# 3. Create windows across the genome for analysis
# ------------------------

# Create 10kb non-overlapping windows over the reference genome
${bedtools_path} makewindows -g genome.sizes -w 10000 > genome.10kb.windows.bed

# ------------------------
# 4. Extract broad repeat class for each feature (use last element if comma-separated)
# ------------------------

awk -F'\t' '{
  n=split($4, arr, ",");
  print $1"\t"$2"\t"$3"\t"arr[n]
}' repeat_overlaps_reformatted_final.bed > repeats_majorclass.bed

# List all unique repeat classes in the newly created BED
awk '{print $4}' repeats_majorclass.bed | sort -u

# ------------------------
# 5. Map each repeat subtype to a broad class using Python
# ------------------------

# Create a mapping script for broad repeat categories
# (save as "broad_repeat_category_mapping.py")
cat << EOF > broad_repeat_category_mapping.py
import pandas as pd

# Mapping dictionary for detailed repeat classes to broad classes
mapping = {
    "LINE": "LINE", "LINE.inc": "LINE", "LINE?": "LINE",
    "DNA": "DNA", "DNA.inc": "DNA", "DNA?": "DNA",
    "LTR": "LTR", "LTR.inc": "LTR", "LTR?": "LTR",
    "Low_complexity": "Low_complexity",
    "Other": "Other", "RC": "Other", "RC?": "Other",
    "Retroposon": "Other", "Retroposon?": "Other",
    "Segmental": "Other", "Unknown": "Unknown", "Unknown.inc": "Unknown",
    "Unspecified": "Other", "SINE": "SINE", "SINE?": "SINE",
    "Satellite": "Satellite", "Simple_repeat": "Simple_repeat",
    "rRNA": "Other", "scRNA": "Other", "snRNA": "Other",
    "srpRNA": "Other", "tRNA": "tRNA", "ARTEFACT": "Other"
}

# Read repeats BED and map to broad classes
df = pd.read_csv('repeats_majorclass.bed', sep='\t', header=None, names=['chrom','start','end','repeatclass'])
df['repeatclass_fixed'] = df['repeatclass'].map(mapping)
# If mapping failed, set as "Other"
df['repeatclass_fixed'] = df['repeatclass_fixed'].fillna('Other')
df[['chrom','start','end','repeatclass_fixed']].to_csv('repeats_majorclass_fixed.bed', sep='\t', header=False, index=False)
EOF

# Run the python mapping script (assumes pandas is installed)
python3 broad_repeat_category_mapping.py

# Remove prefix "HemMar#1#" from chromosome names (if present)
sed -i 's/HemMar#1#//g' repeats_majorclass_fixed.bed

# ------------------------
# 6. Compute coverage (fraction of window covered by each repeat type) in 10kb windows
# ------------------------

# For each unique (broad) repeat type, create a file and compute window-wise coverage
for type in $(cut -f4 repeats_majorclass_fixed.bed | sort | uniq); do
  awk -v t="$type" '$4==t' repeats_majorclass_fixed.bed > "repeats_${type}.bed"
  ${bedtools_path} coverage -a genome.10kb.windows.bed -b repeats_${type}.bed \
    | awk -v type="$type" 'BEGIN{OFS="\t"}{print $1,$2,$3,type,$7}' \
    > "density_${type}_10kb.bed"
done

# ------------------------
# 7. Repeat analysis for 200kb windows
# ------------------------

# Create 200kb windows
${bedtools_path} makewindows -g genome.sizes -w 200000 > genome.200kb.windows.bed

# For each repeat class, compute coverage in 200kb windows
for type in $(cut -f4 repeats_majorclass_fixed.bed | sort | uniq); do
  awk -v t="$type" '$4==t' repeats_majorclass_fixed.bed > "repeats_${type}.bed"
  ${bedtools_path} coverage -a genome.200kb.windows.bed -b repeats_${type}.bed \
    | awk -v type="$type" 'BEGIN{OFS="\t"}{print $1,$2,$3,type,$7}' \
    > "density_${type}_200kb.bed"
done

# ------------------------
# 8. Compute coverage for all repeats combined ("total repeats")
# ------------------------

# Combine all repeats to one .bed file
cat repeats_majorclass_fixed.bed > repeats_ALL.bed

# Compute fraction of each 10kb window covered by any repeat
${bedtools_path} coverage -a genome.10kb.windows.bed -b repeats_ALL.bed \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"total_repeats",$7}' > density_total_repeat_10kb.bed

# Compute fraction of each 200kb window covered by any repeat
${bedtools_path} coverage -a genome.200kb.windows.bed -b repeats_ALL.bed \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"total_repeats",$7}' > density_total_repeat_200kb.bed

# outputs are the  fraction (percent or proportion) of bases in that window covered by a particular repeat class 

## trying instead to get the count of each repeat class in the different windwos

# for 10kb windows
for type in $(cut -f4 repeats_majorclass_fixed.bed | sort | uniq); do
  awk -v t="$type" '$4==t' repeats_majorclass_fixed.bed > "repeats_${type}.bed"
  ${bedtools_path} intersect -c -a genome.10kb.windows.bed -b repeats_${type}.bed \
    | awk -v type="$type" 'BEGIN{OFS="\t"}{print $1,$2,$3,type,$NF}' \
    > "count_${type}_10kb.bed"
done

# for 200 kb windows

for type in $(cut -f4 repeats_majorclass_fixed.bed | sort | uniq); do
  awk -v t="$type" '$4==t' repeats_majorclass_fixed.bed > "repeats_${type}.bed"
  ${bedtools_path} intersect -c -a genome.200kb.windows.bed -b repeats_${type}.bed \
    | awk -v type="$type" 'BEGIN{OFS="\t"}{print $1,$2,$3,type,$NF}' \
    > "count_${type}_200kb.bed"
done

# for all repeats combined

${bedtools_path} intersect -c -a genome.10kb.windows.bed -b repeats_ALL.bed \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"total_repeats",$NF}' > count_total_repeat_10kb.bed

${bedtools_path} intersect -c -a genome.200kb.windows.bed -b repeats_ALL.bed \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"total_repeats",$NF}' > count_total_repeat_200kb.bed


```
