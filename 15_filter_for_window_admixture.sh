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


# filter for windowed admixture
# up to 20% missing
# mac = 2
# no outgroup
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep ingroup.txt --max-missing 0.8 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/14_50kbp_admixture/${input_array}

# bgzip and index the vcf files
bgzip -c ${workdir}/14_50kbp_admixture/${input_array}.recode.vcf > ${workdir}/14_50kbp_admixture/${input_array}.recode.vcf.gz

tabix -p vcf ${workdir}/14_50kbp_admixture/${input_array}.recode.vcf.gz


