# another PCA without the outgroup as requested by a reviewer

source activate bcftools

vcftools --vcf nothin.vcf --max-missing 0.9 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 2 \
--keep ingroup.txt --max-maf 0.49 --thin 5000 --remove-indels --recode --recode-INFO-all --out pca_ingroup

# make chromosome map for the vcfs
grep -v "#" pca_ingroup.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the vcf to plink format
vcftools --vcf pca_ingroup.recode.vcf  --plink --chrom-map chrom_map.txt --out pca 

# convert  with plink
plink --file pca --recode12 --allow-extra-chr --out pca_plink

# run pca 
plink --file pca_plink --pca --allow-extra-chr --out pca_plink_pca
