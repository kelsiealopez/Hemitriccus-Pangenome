#!/bin/bash

#SBATCH -N 1
#SBATCH -t 2-24:00
#SBATCH --mem=400000
#SBATCH -n 20
#SBATCH -p intermediate
#SBATCH -J snakemake
#SBATCH --mail-type=END
#SBATCH -o output_hifiasm_6371_%j.out
[kelsielopez@boslogin07 workflow]$ wc -l snakemake_hifiasm_HMRG_6371.sh
13 snakemake_hifiasm_HMRG_6371.sh
[kelsielopez@boslogin07 workflow]$ 
[kelsielopez@boslogin07 workflow]$ 
[kelsielopez@boslogin07 workflow]$ 
[kelsielopez@boslogin07 workflow]$ head snakemake_hifiasm_HMRG_6371.sh -n 15
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
