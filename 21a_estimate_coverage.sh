#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=cont_depth
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

samtools depth -a KU31900_C_virens_final.bam KU31901_C_BXsordidulus_final.bam KU31902_C_sordidulus_final.bam \
KU31903_C_F2_final.bam KU31904_C_sordidulus_final.bam KU31905_C_virens_final.bam KU31906_C_sordidulus_final.bam \
KU31907_C_sordidulus_final.bam KU31908_C_virens_final.bam KU31912_C_virens_final.bam KU31913_C_virens_final.bam \
KU31914_C_sordidulus_final.bam KU31915_C_sordidulus_final.bam KU31916_C_sordidulus_final.bam KU31918_C_sordidulus_final.bam \
KU31919_C_sordidulus_final.bam KU31920_C_sordidulus_final.bam KU31921_C_sordidulus_final.bam KU31922_C_sordidulus_final.bam \
KU31923_C_sordidulus_final.bam KU31945_C_virens_final.bam KU31946_C_virens_final.bam KU31948_C_virens_final.bam \
KU37529_C_sordidulus_final.bam KU37797_C_virens_final.bam KU37798_C_virens_final.bam KU37801_C_virens_final.bam \
KU37802_C_virens_final.bam KU37803_C_put-hybrid_final.bam KU37805_C_put-hybrid_final.bam KU37806_C_put-hybrid_final.bam \
KU37809_C_put-hybrid_final.bam KU37810_C_put-hybrid_final.bam KU8948_C_pertinax_final.bam KU8990_C_pertinax_final.bam > \
contopus_coverage.txt

# after the above script finishes
# break up the depth files into single column files for each individual (locations dropped)

while read -r name1 number1; do
	number2=$((number1 + 2));
  cut contopus_coverage.txt -f $number2 > ${name1}_depth.txt;
done < popmap.txt
