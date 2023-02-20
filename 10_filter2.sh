#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter2
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-32

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )

# define main working directory
workdir=/lustre/scratch/jmanthey/15_contopus


# find fixed differences (only mac >= 15)
# (filter to start search for diagnostic SNPs)
# snps
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --max-missing 0.9 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 15 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/07_find_fixed_differences/${input_array}_snps
#indels
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --max-missing 0.9 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 15 --max-maf 0.49 --keep-only-indels --recode --recode-INFO-all --out ${workdir}/07_find_fixed_differences/${input_array}_indels

# bgzip and index the vcf files
# snps
bgzip -c ${workdir}/07_find_fixed_differences/${input_array}_snps.recode.vcf > ${workdir}/07_find_fixed_differences/${input_array}_snps.recode.vcf.gz

tabix -p vcf ${workdir}/07_find_fixed_differences/${input_array}_snps.recode.vcf.gz
#indels
bgzip -c ${workdir}/07_find_fixed_differences/${input_array}_indels.recode.vcf > ${workdir}/07_find_fixed_differences/${input_array}_indels.recode.vcf.gz

tabix -p vcf ${workdir}/07_find_fixed_differences/${input_array}_indels.recode.vcf.gz

# run bcftools to simplify the vcftools output for the 10kbp spacing for each dataset
# snps
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/07_find_fixed_differences/${input_array}_snps.recode.vcf.gz > ${workdir}/07_find_fixed_differences/${input_array}_snps.simple.vcf
#indels
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/07_find_fixed_differences/${input_array}_indels.recode.vcf.gz > ${workdir}/07_find_fixed_differences/${input_array}_indels.simple.vcf


# per species up to one missing individual for recombination rate (only biallelic)
#eastern
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep eastern.txt --max-missing 0.9 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/08_relernn/${input_array}_eastern
#western
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep western.txt --max-missing 0.9 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/08_relernn/${input_array}_western


# per species no missing for stairway plots (only biallelic and invariant)
# run twice to get the invariant counts
# variant
#eastern
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep eastern.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_eastern_variant
#western
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep western.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_western_variant
# variant plus invariant
#eastern
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep eastern.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_eastern_all
#western
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep western.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_western_all


# two species (only non-admixed individuals) for fastsimcoal (only biallelic and invariant)
# run twice to get the invariant counts
# variant
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep nonadmixed.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_nonadmixed_variant
# variant plus invariant
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --keep nonadmixed.txt --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_demography/${input_array}_nonadmixed_all






