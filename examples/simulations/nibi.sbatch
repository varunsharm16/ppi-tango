#!/bin/bash
# Alliance Canada (nibi) job for the demo simulations: sbatch nibi.sbatch <simulate.py arguments>
#SBATCH --account=def-yourpi         # your own allocation
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.log
export PATH=$(echo $PATH | tr ":" "\n" | grep -v miniconda3 | paste -sd:)
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.6 openmm/8.4.0 python/3.11 mpi4py
source ~/scratch/ppi-tango/env/bin/activate   # a virtualenv with openmm, openmmtools, pymbar, mdtraj
python simulate.py "$@"
