#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter3
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-32

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )

# define main working directory
workdir=/lustre/scratch/jmanthey/15_contopus

# two species (only non-admixed individuals) for stairway and fastsimcoal (only biallelic and invariant)
# run twice to get the invariant counts
# variant
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep nonadmixed.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --mac 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_nonadmixed_variant
# variant plus invariant
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep nonadmixed.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_nonadmixed_all

#eastern subset
vcftools --vcf ${workdir}/09_demography/${input_array}_nonadmixed_variant.recode.vcf --keep eastern.txt --mac 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_eastern
#western subset
vcftools --vcf ${workdir}/09_demography/${input_array}_nonadmixed_variant.recode.vcf --keep western.txt --mac 1 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_western





