# Minigraph

## Description
Graph building with Minigraph. 

## Table of Contents
- [Github page](https://github.com/lh3/minigraph)

## Installation
Steps to install and set up the project. Include code snippets for installation commands.

```bash
# Install with conda
module load python/3.10.9-fasrc01
source activate python_env1
conda install bioconda::minigraph
```

This is simple to run. We want to input the haplotypes in a random order.
```bash
# Working directory
/n/holyscratch01/edwards_lab/Users/kelsielopez/svim/correct_ref

# Path to directory that has reference
ref_path="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"
# Reference fasta name
#ref_fa="HemMar_prefixed_scaffolds_corrected.fa"

# Path to my directory that has all my hap1 and hap2 assemblies
haps_path="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"

# Working directory, or analysis path 
analysis_path="/n/holyscratch01/edwards_lab/Users/kelsielopez/minigraph"

ref="HemMar_prefixed_scaffolds_corrected"

# in random order
sample10="HMRG_6371.hap2.p_ctg"
sample4="HMRG_6386.hap1.p_ctg"
sample1="HMRG_6386.hap2.p_ctg"
sample9="HMRG_6388.hap1.p_ctg"
sample3="HMRG_6388.hap2.p_ctg"
sample2="HMRG_6431.hap1.p_ctg"
sample11="HMRG_6431.hap2.p_ctg"
sample6="HMRG_6433.hap1.p_ctg"
sample5="HMRG_6433.hap2.p_ctg"
sample7="HMRG_6371.hap1.p_ctg"
sample12="VEFL_149044.hap1.p_ctg"
sample8="VEFL_149044.hap2.p_ctg"

minigraph -cxggs -t20 ${ref_path}/${ref}.fa ${haps_path}/${sample1}.fa > ${analysis_path}/${sample1}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample1}.gfa ${haps_path}/${sample2}.fa > ${analysis_path}/${sample2}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample2}.gfa ${haps_path}/${sample3}.fa > ${analysis_path}/${sample3}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample3}.gfa ${haps_path}/${sample4}.fa > ${analysis_path}/${sample4}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample4}.gfa ${haps_path}/${sample5}.fa > ${analysis_path}/${sample5}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample5}.gfa ${haps_path}/${sample6}.fa > ${analysis_path}/${sample6}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample6}.gfa ${haps_path}/${sample7}.fa > ${analysis_path}/${sample7}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample7}.gfa ${haps_path}/${sample8}.fa > ${analysis_path}/${sample8}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample8}.gfa ${haps_path}/${sample9}.fa > ${analysis_path}/${sample9}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample9}.gfa ${haps_path}/${sample10}.fa > ${analysis_path}/${sample10}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample10}.gfa ${haps_path}/${sample11}.fa > ${analysis_path}/${sample11}.gfa
minigraph -cxggs -t20 ${analysis_path}/${sample11}.gfa ${haps_path}/${sample12}.fa > ${analysis_path}/final_with_outgroup.gfa

# Call variants from bubbles and put them in a bed file format
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample1}.fa > ${sample1}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample2}.fa > ${sample2}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample3}.fa > ${sample3}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample4}.fa > ${sample4}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample5}.fa > ${sample5}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample6}.fa > ${sample6}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample7}.fa > ${sample7}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample8}.fa > ${sample8}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample9}.fa > ${sample9}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample10}.fa > ${sample10}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample11}.fa > ${sample11}.bed
minigraph -cxasm --call -t8 ${analysis_path}/final_with_outgroup.gfa ${haps_path}/${sample12}.fa > ${sample12}.bed

# Merge all bed files 
paste \
HMRG_6371.hap1.p_ctg.bed \
HMRG_6371.hap2.p_ctg.bed \
HMRG_6386.hap1.p_ctg.bed \
HMRG_6386.hap2.p_ctg.bed \
HMRG_6388.hap1.p_ctg.bed \
HMRG_6388.hap2.p_ctg.bed \
HMRG_6431.hap1.p_ctg.bed \
HMRG_6431.hap2.p_ctg.bed \
HMRG_6433.hap1.p_ctg.bed \
HMRG_6433.hap2.p_ctg.bed \
VEFL_149044.hap1.p_ctg.bed \
VEFL_149044.hap2.p_ctg.bed \
| k8 /n/home03/kelsielopez/minigraph/misc/mgutils.js merge - > Minigraph_12haps.bed

# Convert merged BED to VCF

/n/netscratch/edwards_lab/Lab/kelsielopez/pangene/pangene-1.1-bin/bin_x64-linux/k8 \
/n/home03/kelsielopez/minigraph/misc/mgutils-es6.js \
merge2vcf \
-s sample_names.txt \
Minigraph_12haps.bed > Minigraph_12haps.vcf


```
