# Identifying inversions with Syri Svim PGGB Minigraph



### svim-asm

```bash
#!/bin/bash
#SBATCH -p test
#SBATCH -c 8
#SBATCH -t 0-12:00
#SBATCH -o svim_redo_%j.out
#SBATCH -e svim_redo_%j.err 
#SBATCH --mem=100000
#SBATCH --mail-type=END


# Path to directory that has reference
ref_path="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"
# Reference fasta name
ref_fa="HemMar_out_JBAT.FINAL_renamed_sorted.fa"

# Path to my directory that has all my hap1 and hap2 assemblies
haps_path="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"

# Working directory, or analysis path 
analysis_path="/n/holyscratch01/edwards_lab/Users/kelsielopez/svim_redo"

#Text file that has my sample names
sample_names="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/sample_names.txt"


# Diploid mode
cd ${analysis_path}

cat ${sample_names} |
while read LINE;
do 
######################################################################################################################
minimap2 -a -x asm5 --cs -r2k -t 4 ${ref_path}/${ref_fa} ${haps_path}/${LINE}.hap1.p_ctg.fa > ${LINE}.aln.hap1.sam # align hap 1 to reference
minimap2 -a -x asm5 --cs -r2k -t 4 ${ref_path}/${ref_fa} ${haps_path}/${LINE}.hap2.p_ctg.fa > ${LINE}.aln.hap2.sam # align hap 2 to reference
samtools sort -m4G -@4 -o ${LINE}.aln.hap1.sorted.bam ${LINE}.aln.hap1.sam # sort sam 
samtools sort -m4G -@4 -o ${LINE}.aln.hap2.sorted.bam ${LINE}.aln.hap2.sam # sort sam 
samtools index ${LINE}.aln.hap1.sorted.bam # index
samtools index ${LINE}.aln.hap2.sorted.bam # index
svim-asm diploid ${analysis_path} ${LINE}.aln.hap1.sorted.bam ${LINE}.aln.hap2.sorted.bam ${ref_path}/${ref_fa} # run svim-asm
mv variants.vcf ${LINE}.variants.vcf # rename vcf file so it is not overwritten
######################################################################################################################
done

```



```bash
# Merge files with SURVIVOR

in_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/svim_redo"
file_list="/n/holyscratch01/edwards_lab/Users/kelsielopez/svim_redo/file_list.txt"

/n/holyscratch01/edwards_lab/Users/kelsielopez/SURVIVOR/Debug/SURVIVOR merge ${file_list} 500 1 1 1 0 50 ${in_dir}/merged_SURVIVOR_svim_INV.vcf &

nano file_list.txt

HMRG_6371.variants.INV.vcf
HMRG_6386.variants.INV.vcf
HMRG_6388.variants.INV.vcf
HMRG_6431.variants.INV.vcf
HMRG_6433.variants.INV.vcf


grep -c "incomplete_inversion" HMRG_6371.variants.INV.vcf

# remove all incomplete inversions 
sed -i '/incomplete_inversion/d' HMRG_6371.variants.INV.vcf
sed -i '/incomplete_inversion/d' HMRG_6386.variants.INV.vcf
sed -i '/incomplete_inversion/d' HMRG_6388.variants.INV.vcf
sed -i '/incomplete_inversion/d' HMRG_6431.variants.INV.vcf
sed -i '/incomplete_inversion/d' HMRG_6433.variants.INV.vcf



in_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/svim_redo"
file_list="/n/holyscratch01/edwards_lab/Users/kelsielopez/svim_redo/file_list.txt"

/n/holyscratch01/edwards_lab/Users/kelsielopez/SURVIVOR/Debug/SURVIVOR merge ${file_list} 500 1 1 1 0 50 ${in_dir}/merged_SURVIVOR_svim_INV_complete_invs.vcf &


# add sample name so i know which sample is which when i merge them together 
sed 's/Sample/HMRG_6371/' HMRG_6371.variants.INV.vcf > sample_name_HMRG_6371.variants.INV.complete.inv.vcf
sed 's/Sample/HMRG_6386/' HMRG_6386.variants.INV.vcf > sample_name_HMRG_6386.variants.INV.complete.inv.vcf
sed 's/Sample/HMRG_6388/' HMRG_6388.variants.INV.vcf > sample_name_HMRG_6388.variants.INV.complete.inv.vcf
sed 's/Sample/HMRG_6431/' HMRG_6431.variants.INV.vcf > sample_name_HMRG_6431.variants.INV.complete.inv.vcf
sed 's/Sample/HMRG_6433/' HMRG_6433.variants.INV.vcf > sample_name_HMRG_6433.variants.INV.complete.inv.vcf



nano file_list_sample_names.txt

in_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/svim_redo"
file_list="/n/holyscratch01/edwards_lab/Users/kelsielopez/svim_redo/file_list_sample_names.txt"

/n/holyscratch01/edwards_lab/Users/kelsielopez/SURVIVOR/Debug/SURVIVOR merge ${file_list} 500 1 1 1 0 50 ${in_dir}/merged_SURVIVOR_svim_INV_complete_invs_w_names.vcf &


```



### Syri
```bash
# first ragTag

# First map with Minimap

#!/bin/bash
#SBATCH -p test
#SBATCH -c 8
#SBATCH -t 0-12:00
#SBATCH -o minimap_syri_%j.out
#SBATCH -e minimap_syri_%j.err 
#SBATCH --mem=150000
#SBATCH --mail-type=END

# Set the reference genome file
refgenome="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/HemMar_out_JBAT.FINAL_renamed_sorted.fa"

# Set the directory containing the filtered query genomes
query_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/ragtag_redo/top_scaffs_consistent"

# Set the output directory
out_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/syri_redo"

# Loop through each query genome matching the pattern
for qrygenome in "$query_dir"/filtered_*.ragtag.scaffold.fasta;
do
    # Extract the base name of the query genome file (without path and extension)
    base_name=$(basename "$qrygenome" .ragtag.scaffold.fasta)

    # Define the output SAM file name
    out_sam="${base_name}.sam"

    # Run minimap2 with the reference genome and query genome, saving the output to the SAM file
    minimap2 -ax asm5 --eqx "$refgenome" "$qrygenome" > "$out_dir"/"$out_sam"

    # Print a message indicating the file has been processed
    echo "Processed $qrygenome, output saved to $out_sam"
done



#!/bin/bash
#SBATCH -p test
#SBATCH -c 8
#SBATCH -t 0-12:00
#SBATCH -o syri_%j.out
#SBATCH -e syri_%j.err 
#SBATCH --mem=150000
#SBATCH --mail-type=END

# Set the reference genome file
refgenome="/n/holyscratch01/edwards_lab/Users/kelsielopez/ragtag_redo/top_scaffs_consistent/HemMar_out_JBAT.FINAL_renamed_sorted.fa"

# Set the directory containing the SAM files
sam_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/syri_redo/filtered_sam_files"

# Set the directory containing the query genomes
query_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/ragtag_redo/top_scaffs_consistent"

# Loop through each SAM file matching the pattern
for samfile in "$sam_dir"/filtered_*.sam; 
do
    # Extract the base name of the SAM file (without path and extension)
    base_name=$(basename "$samfile" .sam)

    # Define the corresponding query genome file name
    qrygenome="$query_dir/${base_name}.ragtag.scaffold.fasta"

    # Run syri with the SAM file, reference genome, and query genome
    syri -c "$samfile" -r "$refgenome" -q "$qrygenome" -k -F S --prefix ${base_name}_ref.
    
    # Print a message indicating the file has been processed
    echo "Processed $samfile with $qrygenome"
done

```



```bash
# PGGB INVERSIONS
```


```bash
# Minigraph INVERSIONS
```

```bash

# merge and find unique from all 4 programs
```
