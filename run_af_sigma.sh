#!/bin/bash
#SBATCH -n 1 -c 8
###SBATCH -N 1
##SBATCH --gpus-per-task=1#snic2021-5-229 snic2021-5-373 #liu-compute-2020-10
#SBATCH -A liu-compute-2022-22
#SBATCH -t 960


AF_PATH='/proj/wallner/apps/AF_multitemplate'
PATH=$PATH:/proj/wallner/apps/hhsuite/bin/:/proj/wallner/apps/hmmer-3.2.1/bin/:/proj/wallner/apps/kalign/src/
module load Python/3.7.0-anaconda-5.3.0-extras-nsc1
#module load buildenv-gcccuda/.11.1-9.3.0-bare
module load buildenv-gcccuda/11.4-9.3.0-bare
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/sse/manual/CUDA/11.2.1_460.32.03/lib64/:/proj/wallner/cuda11.2_cudnn8//lib64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/sse/manual/CUDA/11.4.2_470.57.02/lib64/:/proj/wallner/cuda11.2_cudnn8//lib64/
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=3.0

source activate /proj/wallner/users/x_bjowa/.conda/envs/alphafold/

which python
nvidia-smi
date
echo Running CMD: python $AF_PATH/run_alphafold.py $@
python $AF_PATH/run_alphafold.py $@
date



 
