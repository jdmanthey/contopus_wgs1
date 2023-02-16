#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=satsuma
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=128
#SBATCH --time=48:00:00
#SBATCH --mem=450G

module load gcc

cd /lustre/scratch/jmanthey/15_contopus/20_satsuma

export SATSUMA2_PATH=/home/jmanthey/satsuma2/bin/

/home/jmanthey/satsuma2/bin/Chromosemble -t /lustre/scratch/jmanthey/15_contopus/20_satsuma/GCF_009829145.1_bChiLan1.pri_genomic.fna -q /lustre/scratch/jmanthey/15_contopus/20_satsuma/GCF_003031625.1_ASM303162v1_genomic.fna -o flycatcher_chromosemble -n 120 -s 1
