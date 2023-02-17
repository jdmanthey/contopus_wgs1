# start interactive session
interactive -p quanah -c 4 -m 8G

# move to directory
workdir=/lustre/scratch/jmanthey/15_contopus
cd ${workdir}/05_structure

# combine all chromosome vcfs into one vcf for the two different datasets
# 5kbp thinning
grep "#" CHR_30_structure_5kbpthin.recode.vcf > 5kbpthin.vcf
for i in $( ls *structure_5kbpthin* ); do grep -v "#" $i >> 5kbpthin.vcf; done
# no thinning
grep "#" CHR_30_structure_nothin.recode.vcf > nothin.vcf
for i in $( ls *structure_nothin* ); do grep -v "#" $i >> nothin.vcf; done


# make chromosome map for the vcfs
grep -v "#" 5kbpthin.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt


# run vcftools for the combined vcfs
# 5kbp thinning
vcftools --vcf 5kbpthin.vcf  --plink --chrom-map chrom_map.txt --out 5kbpthin 
# no thinning
vcftools --vcf nothin.vcf  --plink --chrom-map chrom_map.txt --out nothin 


# convert  with plink
# 5kbp thinning
plink --file ${workdir}/05_structure/5kbpthin --recode12 --allow-extra-chr --out ${workdir}/05_structure/5kbpthin_plink
# no thinning
plink --file ${workdir}/05_structure/nothin --recode12 --allow-extra-chr --out ${workdir}/05_structure/nothin_plink


# run pca on each dataset
# 5kbp thinning
plink --file ${workdir}/05_structure/5kbpthin_plink --pca --allow-extra-chr --out ${workdir}/05_structure/5kbpthin_plink_pca
# no thinning
plink --file ${workdir}/05_structure/nothin_plink --pca --allow-extra-chr --out ${workdir}/05_structure/nothin_plink_pca


# run admixture on each dataset
# 5kbp thinning
for K in 3; do admixture --cv ${workdir}/05_structure/5kbpthin_plink.ped $K; done
# no thinning
for K in 3; do admixture --cv ${workdir}/05_structure/nothin_plink.ped $K; done



