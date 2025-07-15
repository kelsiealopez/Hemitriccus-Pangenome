## Pangenome graph building with pggb 

Resources:
- [pggb Github page](https://github.com/pangenome/pggb)
- [pggb nextflow](https://nf-co.re/pangenome/1.0.0/)

### Run pggb using nextflow nf-core
I could not get `pggb` to work when I downloaded it using conda, singluarity, manually, etc. Nextflow was the only thing to get the graph building step to work for me. I must have some conflict with some dependencies. However, pggb downloaded with conda worked for the initial partitioning steps, but the graph building steps (`wfmash`, `seqwish`, `smoothxg`, and `gfaffix`) did not work for me unless I used nf-core.
```bash
# Install with conda
module load python/3.10.9-fasrc01
conda create --name env_nf nextflow
conda activate env_nf
```

## 1. First we need to partition sequences into communities. 

I followed [this](https://gtpb.github.io/CPANG22/pages/Day2a_Homo_sapiens_pangenome_graphs) tutorial first to partition my sequences. 

Partition contigs by chromosome by mapping each assembly against the scaffolded references:

```bash
source activate python_env1

PGGB_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb"
PATH_SAMPLE_FA_GZ="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"

# You will need to bgzip your input FASTA file before running the pipeline
cd ${PATH_SAMPLE_FA_GZ}

for file in *.p_ctg.fa;
do
bgzip -@ 4 $file
done

# Optionally, you can directly index the bgzipped FASTA file with

for file in *.p_ctg.fa.gz;
do
samtools faidx $file
done

# Using -s 5k -p 85 -N -m -t 8
# i used -p 85 because my outgroup is 14% divergent from Hemitriccus

PATH_REFERENCE_FA_GZ="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/HemMar_prefixed_scaffolds_corrected.fa.gz"

mkdir -p /n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning

cd ${PATH_SAMPLE_FA_GZ}

ls *.p_ctg.fa.gz | while read FASTA; do
  NAME=$(basename $FASTA .p_ctg.fa.gz);
  echo $NAME

  PATH_PAF="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning/$NAME.vs.ref.paf"
  wfmash $PATH_REFERENCE_FA_GZ ${PATH_SAMPLE_FA_GZ}/$FASTA -s 5k -p 85 -N -m -t 8 > $PATH_PAF
done
```

Collect unmapped contigs and remap them in split mode:

```bash
cd ${PATH_SAMPLE_FA_GZ}

ls *.p_ctg.fa.gz | while read FASTA; do
  NAME=$(basename $FASTA .p_ctg.fa.gz);
  echo $NAME

  PATH_PAF="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning/$NAME.vs.ref.paf"
  comm -23 <(cut -f 1 $FASTA.fai | sort) <(cut -f 1 $PATH_PAF | sort) > $NAME.unaligned.txt

  wc -l $NAME.unaligned.txt
  if [[ $(wc -l $NAME.unaligned.txt | cut -f 1 -d\ ) != 0 ]];
  then 
    samtools faidx $FASTA $(tr '\n' ' ' < $NAME.unaligned.txt) > $NAME.unaligned.fa
    samtools faidx $NAME.unaligned.fa

    PATH_NO_SPLIT_PAF="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning/$NAME.vs.ref.no_split.paf"
    wfmash $PATH_REFERENCE_FA_GZ ${PATH_SAMPLE_FA_GZ}/$NAME.unaligned.fa -s 5k -p 85 -m -t 8 > $PATH_NO_SPLIT_PAF
  fi
done
```

Collect our best mapping for each of our attempted split rescues:

```bash
cd /n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning

ls *.vs.ref.no_split.paf | while read PAF; do
  cat $PAF | awk -v OFS='\t' '{ print $1,$11,$0 }' | sort -n -r -k 1,2 | \
    awk -v OFS='\t' '$1 != last { print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15); last = $1; }'
done > rescues.paf
```

Subset by chromosome, including the references:

```bash
cd ${PATH_SAMPLE_FA_GZ}

DIR_PARTITIONING="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning"

mkdir -p parts

( seq 34 ) | while read i; do
  awk '$6 ~ "scaffold_'$i'$"' $(ls $DIR_PARTITIONING/*.vs.ref.paf | \
    grep -v unaligned | sort; echo $DIR_PARTITIONING/rescues.paf) | cut -f 1 | sort | uniq \
    > parts/scaffold_$i.contigs;
done

( seq 34 ) | while read i; do
  echo scaffold_$i
  samtools faidx HemMar_prefixed_scaffolds_corrected.fa.gz HemMar#1#scaffold_$i > parts/scaffold_$i.pan.fa

  ls *.p_ctg.fa.gz | while read FASTA; do
    NAME=$(basename $FASTA *.p_ctg.fa.gz);
    echo scaffold_$i $NAME

    samtools faidx $FASTA $( comm -12 <(cut -f 1 $FASTA.fai | sort) <(sort parts/scaffold_$i.contigs) ) >> parts/scaffold_$i.pan.fa
  done

  bgzip -@ 8 parts/scaffold_$i.pan.fa && samtools faidx parts/scaffold_$i.pan.fa.gz
done
```

## 2. Run `pggb` using nf-core

```bash

# make sure you are in the nextflow environment and not normal working python environment
conda activate env_nf

in_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/parts"
outdir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb"

cd /n/holyscratch01/edwards_lab/Users/kelsielopez/pggb

ls $in_dir/*.pan.fa.gz | while read FASTA; do
  NAME=${FASTA%.pan.fa.gz};
  echo $NAME
  mkdir $NAME
  cd $NAME
    nextflow run nf-core/pangenome -r dev --input $FASTA --n_haplotypes 12 --wfmash_map_pct_id 85 --seqwish_min_match_length 79 --outdir $NAME -profile singularity
  cd ../
  done
```

## 3. Deconstruct graphs into vcf files for each chromosome and remove large variants 

Download `vg deconstruct`:

NEED to run 1.40.0. I had problems down the line when I ran vcf wave when I was using different versions of vg. 

https://github.com/vgteam/vg

Download 1.40.0 from github! 
https://github.com/vgteam/vg/releases?page=2

```bash
# Do NOT be in python env of nf_env
conda deactivate

chmod +x vg

# Add to path  
export PATH=/n/home03/kelsielopez/vg:$PATH
```

Download `vcfbub`:

vcfbub v0.1.0
https://github.com/pangenome/vcfbub

```bash
# make executable and add to path 
chmod u+x vcfbub
```

Run `vg deconstruct` , `vcfbub` , and `vcfwave` for each chromosome 

```bash
# change according to which samples are done so far... I just started this process while the rest of the pangenome graphs were building because that can take a long time, especially for larger chromosomes.

samples=(scaffold_10 scaffold_11 scaffold_12 scaffold_13 scaffold_14 scaffold_15 scaffold_16 scaffold_17 scaffold_18 scaffold_19)

for sample in "${samples[@]}"; 
do

vg_input_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/parts/${sample}/FINAL_GFA"
vg_out_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/deconstruct"

~/vg deconstruct -H '#' ${vg_input_dir}/${sample}.pan.fa.gz.gfaffix.unchop.Ygs.view.gfa -e -a -P "HemMar" > ${vg_out_dir}/${sample}.snarls.vcf

vcfbub_input_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/deconstruct"
vcfbub_out_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfbub"

# vcf file needs to be gzipped  or else vcfbub doesn't work 
gzip ${vcfbub_input_dir}/${sample}.snarls.vcf

~/vcfbub -l 0 -a 100000 -i ${vcfbub_input_dir}/${sample}.snarls.vcf.gz > ${vcfbub_out_dir}/${sample}.vcfbub.snarls.vcf

# Run vcfwave
module load python/3.10.9-fasrc01
source activate python_env1

vcfwave_input_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfbub"
vcfwave_out_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave"

vcfwave -I 1000 -t 1 ${vcfwave_input_dir}/${sample}.vcfbub.snarls.vcf > ${vcfwave_out_dir}/${sample}.vcfwave.vcfbub.snarls.vcf

```

## 4. Merge vcfs, clean, fix vcfs, etc. 

Create a `vcf_fofn.txt` file that has the paths to all the vcf files

```bash 
nano vcf_fofn.txt

# looks like this, etc. 
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_1.vcfwave.vcfbub.snarls.vcf
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_2.vcfwave.vcfbub.snarls.vcf
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_3.vcfwave.vcfbub.snarls.vcf
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_4.vcfwave.vcfbub.snarls.vcf
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_6.vcfwave.vcfbub.snarls.vcf
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_7.vcfwave.vcfbub.snarls.vcf
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_8.vcfwave.vcfbub.snarls.vcf
/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/vcfwave/scaffold_9.vcfwave.vcfbub.snarls.vcf
```

Combining and cleaning steps:

```bash
# I'm going back to python environment because this is how I downloaded most of the programs 
source activate python_env1

# 1. Concatenate vcf files 
bcftools concat -O z -o pggb_merged.vcf.gz -f vcf_fofn.txt

# 2. Fix ploidy 
bcftools +fixploidy pggb_merged.vcf.gz > pggb_merged.fixploidy.vcf

# 3. 
bcftools +fill-tags pggb_merged.fixploidy.vcf -Ob -o pggb_merged_cleanup.bcf -- -t AN,AC,AF

# 4.
PATH_REFERENCE_FA="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/HemMar_prefixed_scaffolds_corrected.fa"

bcftools norm -m+both -f ${PATH_REFERENCE_FA} pggb_merged_cleanup.bcf > pggb_merged_normalized.vcf

# 5. Add tags for allele number, allele count, allele proportion (?) and missing proportion 
bcftools sort pggb_merged_normalized.vcf | bcftools +fill-tags -- -t AN,AC,AF,F_MISSING | bcftools view -o pggb_cleaned_final.vcf.gz -Oz

# 6. Make a vcf of biallelic SNPs
bcftools view -m2 -M2 -v snps -o pggb_cleaned_final_biallelic_snp.vcf.gz -Oz pggb_cleaned_final.vcf.gz 

# 7. Make a vcf of biallelic indels 
bcftools view -m2 -M2 -v indels -o pggb_cleaned_final_biallelic_indels.vcf.gz -Oz pggb_cleaned_final.vcf.gz 
```

## 5. Compute variant table

Without repeat information and genomic overlaps:

```bash
# 1. # VEFL_149044 is my outgroup
bcftools view -s VEFL_149044 -Oz -o pggb_cleaned_outgrouponly.vcf.gz pggb_cleaned_final.vcf.gz

# 2.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' pggb_cleaned_outgrouponly.vcf.gz > aa.tab

# 3.
python3 calc_anc_allele.py

# 4.
bgzip aa.processed.tab 

# 5.
tabix -s1 -b2 -e2 aa.processed.tab.gz

# 6.
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > hdr.txt

# 7.
bcftools annotate -x INFO/CONFLICT,INFO/LV,INFO/PS,INFO/AT,INFO/TYPE,INFO/LEN,INFO/ORIGIN,INFO/NS pggb_cleaned_final.vcf.gz | bcftools annotate -a aa.processed.tab.gz -c CHROM,POS,REF,ALT,.INFO/AA -h hdr.txt -Oz -o pggb_recode_aa.vcf.gz

# 8.
bcftools view -Ou -a -s ^VEFL_149044 pggb_recode_aa.vcf.gz | bcftools annotate -Ou -x INFO/AF,INFO/F_MISSING | bcftools view -c 1 -Ou -Oz -o pggb_ingroup_only.vcf.gz

# 9.
bcftools query -f "%CHROM\t%POS0\t%END\t%TYPE\t%REF\t%ALT\t%INFO/AA\t%INV[\t%GT]\n" pggb_ingroup_only.vcf.gz > pggb_variation.tab

# 10. call variants
python3 compute_var_table_no_regions_no_repeats.py
```

With repeat information and genomic overlaps:

```bash
bcftools view -Ou -a -s ^VEFL_149044 pggb_recode_aa.vcf.gz | bcftools annotate -Ou -x INFO/AF,INFO/F_MISSING | bcftools view -c 1 -Ou | bcftools annotate -a /n/holyscratch01/edwards_lab/Users/kelsielopez/agat/all_regions_blocks_tab.bed -h gr_header.txt -c CHROM,FROM,TO,GR --merge-logic GR:unique -Oz -o pggb_ingroup_only.vcf.gz
bcftools query -f "%CHROM\t%POS0\t%END\t%TYPE\t%INFO/GR\t%REF\t%ALT\t%INFO/AA\t%INV[\t%GT]\n" pggb_ingroup_only.vcf > fixed_pggb_variation.tab
cut -f1,2,3 fixed_pggb_variation.tab | bedtools intersect -a - -b /n/holyscratch01/edwards_lab/Users/kelsielopez/repeats/rep_masker/02_custom_lib/HemMar_prefixed_scaffolds_corrected.repeats.bed -loj | cut -f 1,2,3,7 > pggb_variation_repeat_overlaps.tab
cut -f1,2,3 fixed_pggb_variation.tab | bedtools intersect -a - -b /n/holyscratch01/edwards_lab/Users/kelsielopez/agat/all_regions_blocks_tab.bed -loj | cut -f 1,2,3,7 > fixed_pggb_variation_genomic_overlaps.tab

# going to be overwriting everything
python3 compute_var_table_overlaps_fixed.py
```
