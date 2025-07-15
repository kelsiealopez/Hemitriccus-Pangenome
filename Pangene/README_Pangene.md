
### 1. Prepare protein fasta file 
```bash
# Prepare fasta file to look like this 'GENE:transcript'

>GPN1:rna-XM_032682815.1.54
WRRQGEAGRAVLCVFWCWAWPAPGKPPSCRAWPPTCTGSAALRT.SI.TPPCTTCPSPPTSVSGTL.STK
KS.NRSDMG.AQTVE..PLSISLLQGLTR..SSLKKDKMHLSMLLLTHRGKLRYSPGQHQEPS.LRPWLP
LFLQLLSM.WTPLAVLTLSLLCPTCCMPAGSCTRQSYLSL.S.TKLT.LTTALQWNGCRTLRLFRMP.IK
RPPMSVT.LVL.V.CWMNFTVH.RWLVFLRCLAQDWMSFLSSFLKL.MNMRGSIVQNTSA.EKHWRKLKI
NKRESSWNTCGRTWAVCVCRAAHWQDLLMLLQWVPLS.Y.HEELSMKSKKRERVILMTLTMKGLRRVMKN
QPSETLCRTCG.NAKGEATRMNE
>GPN1:rna-XM_032682816.1.54
WRRQGEAGRAVLCVFWCWAWPAPGKPPSCRAWPPTCTGSAALRT.SI.TPPCTTCPSPPTSVSGTL.STK
KS.NRSDMG.AQTVE..PLSISLLQGLTR..SSLKKDKMHLSMLLLTHRGKLRYSPGQHQEPS.LRPWLP

#so basically i want to join the gene name and the transcript name oft hese two files with ':'

#I also have this annotation file which i based these files off of 

740793 /n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/hemMar.toga.merged_protein_output.aa
(python_env1) [kelsielopez@boslogin06 toga]$ head hemMar_with_CDS_corrected.gff3
##gff-version 3
scaffold_1	TOGA	gene	30610	43870	.	-	.	ID=nbis-gene-7993;Name=GPN1,chiLan.LOC116791331
scaffold_1	TOGA	RNA	30610	43606	.	-	.	ID=rna-XM_032682815.1.54;Parent=nbis-gene-7993;Name=GPN1
scaffold_1	TOGA	CDS	30610	30701	.	-	.	ID=rna-XM_032682815.1.54.cds1;Parent=rna-XM_032682815.1.54
scaffold_1	TOGA	exon	30610	30701	.	-	.	ID=rna-XM_032682815.1.54.exon1;Parent=rna-XM_032682815.1.54
scaffold_1	TOGA	CDS	31313	31423	.	-	.	ID=rna-XM_032682815.1.54.cds2;Parent=rna-XM_032682815.1.54
scaffold_1	TOGA	exon	31313	31423	.	-	.	ID=rna-XM_032682815.1.54.exon2;Parent=rna-XM_032682815.1.54
scaffold_1	TOGA	CDS	31903	31993	.	-	.	ID=rna-XM_032682815.1.54.cds3;Parent=rna-XM_032682815.1.54
scaffold_1	TOGA	exon	31903	31993	.	-	.	ID=rna-XM_032682815.1.54.exon3;Parent=rna-XM_032682815.1.54
scaffold_1	TOGA	CDS	32700	32739	.	-	.	ID=rna-XM_032682815.1.54.cds4;Parent=rna-XM_032682815.1.54



# use this python code to just merge the two protein fasta

python3 merge_protein_fastas.py


# need to remove duplicate headers.... because it kept saying it failed to build the index

(python_env1) [kelsielopez@boslogin06 labeled]$ grep '^>' /n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/combined_headers.faa | sort | uniq -c | awk '$1>1'
      2 >CDKL2:XM_015276643.2.-1
      2 >CDKL2:XM_015276644.2.-1
      2 >CDKL2:XM_025150398.1.-1
      2 >CDKL2:XM_025150400.1.-1
      2 >chiLan.LOC116779876:rna-XM_032674383.1.-1
      2 >chiLan.LOC116779928:rna-XM_032674440.1.-1
      2 >chiLan.LOC116781295:rna-XM_032676997.1.-1
      2 >galGal.LOC100859273:XM_025146279.1.-1
      2 >galGal.LOC107050721:XM_025145990.1.-1
      2 >galGal.LOC107050724:XM_025145752.1.-1
      2 >PNKP:rna-XM_032677460.1.-1
      2 >PNKP:rna-XM_032677461.1.-1
      2 >PNKP:rna-XM_032677464.1.-1
      2 >RCHY1:XM_015276333.2.-1
      2 >RTN2:XM_025144514.1.-1
      2 >THUMPD2:IDmodified-rna-33033
      2 >THUMPD2:XM_004935288.3.-1
      2 >THUMPD2:XM_025149411.1.-1


awk '
    /^>/ {
        if (seen[$0]++) next
    }
    { print }
' /n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/combined_headers.faa \
 > /n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/combined_headers_deduped.faa
 
 
```

### 2. Run miniprot for protein alignment to all haplotype assemblies
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

### 3. Run pangene for pangene gene graph building and call CNV and PAV

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
