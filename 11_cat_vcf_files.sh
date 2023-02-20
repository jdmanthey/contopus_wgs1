interactive -p quanah

cd /lustre/scratch/jmanthey/15_contopus/07_find_fixed_differences

grep "#C" CHR_4_snps.recode.vcf > header.txt

cat *_snps.simple* >> finding_fixed_snps.simple.vcf

cat *_indels.simple* >> finding_fixed_indels.simple.vcf
