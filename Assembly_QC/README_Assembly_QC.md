## Use Snakemake workflow for assembly with HifiAsm

#### 1. Takes raw PacBio BAM files, converts them from BAM to FASTQ
#### 2.  Assembles using HifiAsm
#### 3. Converts haplotype assemblies from HifiAsm from GFA to FA
#### 4. Runs QUAST on all output haplotype assemblies

- [Snakemake Github page](https://github.com/harvardinformatics/pacbio_hifi_assembly)
- [HifiAsm Github page](https://github.com/chhylp123/hifiasm)
- [QUAST Github page](https://github.com/ablab/quast)


```bash

#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-24:00
#SBATCH --mem=400000
#SBATCH -n 20
#SBATCH -p intermediate
#SBATCH -J snakemake
#SBATCH --mail-type=END
#SBATCH -o output_hifiasm_6371_%j.out
#SBATCH -e errors_hifiasm_6371_%j.err
#check that snakemake is working by running snakemake --help
snakemake --snakefile Snakefile --cores 20 --use-conda --rerun-incomplete
```
