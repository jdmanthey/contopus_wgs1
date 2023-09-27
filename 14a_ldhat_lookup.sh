#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=lookup
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

# eastern
../complete -n 24 -rhomax 100 -n_pts 101 -theta 0.008



#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=lookup2
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

# western
../complete -n 30 -rhomax 100 -n_pts 101 -theta 0.004
