## Estimating the recombination rate using ReLERNN


```bash
#!/bin/bash

#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH --mem=100000
#SBATCH -c 48
#SBATCH --gres=gpu:4
#SBATCH -p gpu_test
#SBATCH -J ReLERNN_HemMar_gpu
#SBATCH --mail-type=END
#SBATCH -o ReLERNN_HemMar_gpu_%j.out
#SBATCH -e ReLERNN_HemMar_gpu_%j.err

module load python/3.10.9-fasrc01
conda activate recomb
module load python/3.10.9-fasrc01 cuda/11.8.0-fasrc01 cudnn/8.9.2.26_cuda11-fasrc01
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/TensorRT-8.6.1.6/lib
export KERAS_BACKEND="tensorflow"
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/n/sw/helmod-rocky8/apps/Core/cuda/11.8.0-fasrc01/cuda


SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="6.25e-9"
URTR="1"
DIR="./example_output/"
VCF="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb_redo/merging/pggb_cleaned_final_biallelic_snp.vcf"
GENOME="./genome.bed"

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --upperRhoThetaRatio ${URTR} \
    --nTrain 13000 \
    --nVali 2000 \
    --nTest 100 \
    --seed ${SEED}

# Train network
${TRAIN} \
    --projectDir ${DIR} \
    --nCPU 48 \
    --seed ${SEED}

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED}

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${DIR} \
    --nSlice 2 \
    --nReps 2 \
    --seed ${SEED}

```
