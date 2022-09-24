#!/bin/bash

#SBATCH --job-name=1a_bal_k10
#SBATCH --time=04:00:00
#SBATCH --array=1-200
#SBATCH --ntasks=1
#SBATCH --output=Output/%x-%j.out
#SBATCH --error=Errors/arrayJob_%A_%a.err
#SBATCH --mail-user=nayan.saxena@mail.utoronto.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load gcc/7.3.0 r/4.0.2

export R_LIBS_USER=$HOME/Rpackages:$R_LIBS_USER

Rscript ci_sim_051118.R 10 8.5 8.5 0.006 2
