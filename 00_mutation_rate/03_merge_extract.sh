
# pull out the cds regions (in 3 jobs or 1, time depending for the next few steps)
bedtools intersect -a Contopus_KU31919.vcf.gz -b manakin_cds_subset.gff > contopus.vcf

bedtools intersect -a Neopipo_SRR9946662.vcf.gz -b manakin_cds_subset.gff > neopipo.vcf

bedtools intersect -a Onych_SRR9947257.vcf.gz -b manakin_cds_subset.gff > onych.vcf

# extract the header for each 
gunzip -cd Contopus_KU31919.vcf.gz | grep "#" > header_c.vcf

gunzip -cd Neopipo_SRR9946662.vcf.gz | grep "#" > header_n.vcf

gunzip -cd Onych_SRR9947257.vcf.gz | grep "#" > header_o.vcf


# combine header and vcf files
cat header_c.vcf contopus.vcf > contopus2.vcf

cat header_n.vcf neopipo.vcf > neopipo2.vcf

cat header_o.vcf onych.vcf > onych2.vcf

# bgzip and tabix
bgzip contopus2.vcf

tabix contopus2.vcf.gz

bgzip neopipo2.vcf

tabix neopipo2.vcf.gz

bgzip onych2.vcf

tabix onych2.vcf.gz

# merge the vcf files
bcftools merge -m id contopus2.vcf.gz neopipo2.vcf.gz onych2.vcf.gz > genes.vcf

# remove sites with missing data, require a minimum depth of 3
vcftools --vcf genes.vcf --max-missing 1.0 --minDP 3 --recode --recode-INFO-all --out genes_filtered

# simplify the vcf file for memory ease
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' genes_filtered.recode.vcf > genes_filtered.simple.vcf


