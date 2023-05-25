#!/bin/bash
#SBATCH --gpus 1
#SBATCH -t 4320
#SBATCH -A berzelius-2022-218
##SBATCH -A berzelius-2022-216

#change the path to your the location of the repo
#the script that is submitted through sbatch is actually copied to a temporary place
#so it is not possible to use the location of this script to get the location of the repo...

AF_PATH='/proj/wallner-b/apps/AF_multitemplate'
module load Anaconda/2021.05-nsc1
#module load buildenv-gcccuda/11.2-8.3.1-bare
module load buildenv-gcccuda/11.4-8.3.1-bare
#module load gcc/system
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/sse/manual/CUDA/11.4.2_470.57.02/lib64/
#/software/sse/manual/CUDA/11.2.1_460.32.03/lib64/:/proj/wallner/cuda11.2_cudnn8//lib64/
#env
export XLA_FLAGS=--xla_gpu_cuda_data_dir=$CUDA_PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=10.0

echo $XLA_FLAGS
#conda activate /proj/wallner-b/users/x_bjowa/.conda/envs/alphafold/
#conda activate /proj/wallner-b/users/x_bjowa/.conda/envs/alphafold2/
conda activate /proj/beyondfold/apps/.conda/envs/af_server/

which python
nvidia-smi
date
echo Running CMD: python $AF_PATH/run_alphafold.py $@
python $AF_PATH/run_alphafold.py $@
date



