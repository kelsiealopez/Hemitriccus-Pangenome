## Use Snakemake workflow for assembly with HifiAsm

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
