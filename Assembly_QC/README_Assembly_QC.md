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


#### 5. Run BUSCO on all assemblies

```bash
#!/bin/bash
#SBATCH -p test,shared
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH --mem=100000
#SBATCH --mail-type=END
#SBATCH -o hap_assemblies_busco_%A_%a.out
#SBATCH -e hap_assemblies_busco_%A_%a.err

# List of input files
input_files=(
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6371.hap1.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6371.hap2.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6386.hap1.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6386.hap2.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6388.hap1.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6388.hap2.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6431.hap1.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6431.hap2.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6433.hap1.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/HMRG_6433.hap2.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/VEFL_149044.hap1.p_ctg.fa"
"/n/netscratch/edwards_lab/Lab/kelsielopez/hap_assemblies/prefixed/VEFL_149044.hap2.p_ctg.fa"
)

# Get file for current array task
input_file="${input_files[$SLURM_ARRAY_TASK_ID]}"

# Extract the base name without extension for output folder
base_name=$(basename "$input_file" .p_ctg.fa)

# Directory for the BUSCO executable
BUSCO_DIR="/n/netscratch/edwards_lab/Lab/kelsielopez/busco-5.8.3/bin"

# Run BUSCO using aves_odb10 database
$BUSCO_DIR/busco -i "$input_file" -l aves_odb10 -o "${base_name}_busco_aves" -m genome -f

```
