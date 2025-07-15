

```bash
# Prepare fasta file
```

### Run miniprot for protein alignment to all haplotype assemblies
```bash
#!/bin/bash
#SBATCH -p test
#SBATCH -c 16
#SBATCH -t 0-12:00
#SBATCH -o redo_June_miniprot_%j.out
#SBATCH -e redo_June_miniprot_%j.err 
#SBATCH --mem=200000
#SBATCH --mail-type=END

# just align to 
protein_faa="/n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/combined_headers_deduped.faa"
genome_indir="/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed"

workdir="/n/netscratch/edwards_lab/Lab/kelsielopez/miniprot/HemMar/redo/labeled/June_2025_fix"

sample_names="/n/netscratch/edwards_lab/Lab/kelsielopez/miniprot/hap_sample_names.txt"


cd ${workdir}
cat ${sample_names} | while read LINE; do
    miniprot -t16 -d ${workdir}/${LINE}.mpi ${genome_indir}/${LINE}.p_ctg.fa

    miniprot -Iut16 --gff ${workdir}/${LINE}.mpi ${protein_faa} > ${workdir}/${LINE}_aln.gff

    miniprot --outs=0.97 --no-cs -Iut16 ${workdir}/${LINE}.mpi ${protein_faa} > ${workdir}/${LINE}.paf
done


sample_names="/n/netscratch/edwards_lab/Lab/kelsielopez/miniprot/hap_sample_names_6433.txt"

cd ${workdir}
cat ${sample_names} | while read LINE; do
    miniprot -t16 -d ${workdir}/${LINE}.mpi ${genome_indir}/${LINE}.p_ctg.fa

    miniprot -Iut16 --gff ${workdir}/${LINE}.mpi ${protein_faa} > ${workdir}/${LINE}_aln.gff

    miniprot --outs=0.97 --no-cs -Iut16 ${workdir}/${LINE}.mpi ${protein_faa} > ${workdir}/${LINE}.paf
done


```

### Run pangene for pangene gene graph building and call copy number variation and presence absence variation

```bash
#!/bin/bash
#SBATCH -p test
#SBATCH -c 16
#SBATCH -t 0-12:00
#SBATCH -o pangene_redo_June_2025_%j.out
#SBATCH -e pangene_redo_June_2025_%j.err 
#SBATCH --mem=200000
#SBATCH --mail-type=END

cd /n/netscratch/edwards_lab/Lab/kelsielopez/pangene/redo/labeled/June_2025

indir="/n/netscratch/edwards_lab/Lab/kelsielopez/miniprot/HemMar/redo/labeled/June_2025_fix"

pangene \
${indir}/HMRG_6371.hap1.paf \
${indir}/HMRG_6371.hap2.paf \
${indir}/HMRG_6386.hap1.paf \
${indir}/HMRG_6386.hap2.paf \
${indir}/HMRG_6388.hap1.paf \
${indir}/HMRG_6388.hap2.paf \
${indir}/HMRG_6431.hap1.paf \
${indir}/HMRG_6431.hap2.paf \
${indir}/HMRG_6433.hap1.paf \
${indir}/HMRG_6433.hap2.paf \
> HemMar_10_haps_June_2025_fix.gfa

gfa="HemMar_10_haps_June_2025_fix.gfa"

/n/netscratch/edwards_lab/Lab/kelsielopez/pangene/pangene-1.1-bin/bin_x64-linux/k8 \
/n/netscratch/edwards_lab/Lab/kelsielopez/pangene/pangene-1.1-bin/scripts/pangene.js \
gfa2matrix -c ${gfa} > HemMar_10_haps_June_2025_new_faa_fix_PAV.Rtab

# make it readable in excel
cp HemMar_10_haps_June_2025_new_faa_fix_PAV.Rtab HemMar_10_haps_June_2025_new_faa_fix_PAV.Rtab.xls
```
