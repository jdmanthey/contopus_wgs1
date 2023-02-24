cd /lustre/scratch/jmanthey/15_contopus/08_relernn

grep "#" CHR_29_western.recode.vcf > western_relernn.vcf

grep "#" CHR_29_eastern.recode.vcf > eastern_relernn.vcf

for i in $( ls *eastern.recode.vcf ); do grep -v "#" $i >> eastern_relernn.vcf; done

for i in $( ls *western.recode.vcf ); do grep -v "#" $i >> western_relernn.vcf; done
