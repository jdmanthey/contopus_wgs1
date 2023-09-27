cd /lustre/scratch/jmanthey/15_contopus/08_relernn

grep "#" CHR_29_western.recode.vcf > western_ldhat.vcf

grep "#" CHR_29_eastern.recode.vcf > eastern_ldhat.vcf

for i in $( ls *eastern.recode.vcf ); do grep -v "#" $i >> eastern_ldhat.vcf; done

for i in $( ls *western.recode.vcf ); do grep -v "#" $i >> western_ldhat.vcf; done

# bgzip and tabix index vcfs
bgzip eastern_ldhat.vcf
#tabix
tabix -p vcf eastern_ldhat.vcf.gz

bgzip western_ldhat.vcf
#tabix
tabix -p vcf western_ldhat.vcf.gz

# make header files
gunzip -cd eastern_ldhat.vcf.gz | head -n100000 | grep "#" > eastern_header.txt
gunzip -cd western_ldhat.vcf.gz | head -n100000 | grep "#" > western_header.txt

