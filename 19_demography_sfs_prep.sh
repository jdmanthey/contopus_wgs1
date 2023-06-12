# don't use sex chromosomes for demography
rm CHR_W*
rm CHR_Z*

# cat all files for each grouping
grep "#" CHR_10_western.recode.vcf > western.vcf
for i in $( ls *western.recode* ); do cat $i >> western.vcf; done

grep "#" CHR_10_eastern.recode.vcf > eastern.vcf
for i in $( ls *eastern.recode* ); do cat $i >> eastern.vcf; done


grep "#" CHR_10_nonadmixed_all.recode.vcf > all.vcf
for i in $( ls *nonadmixed_all.recode* ); do cat $i >> all.vcf; done


grep "#" CHR_10_nonadmixed_variant.recode.vcf > variant_2pop.vcf
for i in $( ls *nonadmixed_variant.recode* ); do cat $i >> variant_2pop.vcf; done

# remove all the individual chromosome files
rm *recode*

# count the number of sites in each file
grep -v "#" eastern.vcf | wc -l			# 4528955
grep -v "#" western.vcf | wc -l			# 3311371
grep -v "#" all.vcf | wc -l				# 240591606
grep -v "#" variant_2pop.vcf | wc -l	# 6329302


# number of invariant sites (if needed to manually fill for FSC2 SFS):
# 2pop: 240591606 - 8702513 = 231889093
# eastern: 240591606 - 4528955 = 236062651
# western: 240591606 - 6329302 = 234262304

# use at least 60GB RAM

conda activate easySFS

# convert vcf to sfs files
~/easySFS/easySFS.py -i variant_2pop.vcf -p /lustre/scratch/jmanthey/15_contopus/09_demography/sfs_2pop.txt -a -f --proj 24,30

~/easySFS/easySFS.py -i eastern.vcf -p /lustre/scratch/jmanthey/15_contopus/09_demography/sfs_east.txt -a -f --proj 24

~/easySFS/easySFS.py -i western.vcf -p /lustre/scratch/jmanthey/15_contopus/09_demography/sfs_west.txt -a -f --proj 30
