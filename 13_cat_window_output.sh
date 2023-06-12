cd /lustre/scratch/jmanthey/15_contopus/06_window_stats/windows

# combine the output for different analyses into a single file each
# first add a header for each file
grep 'pop1' CHR_10__10000001__10050000__stats.txt > ../window_heterozygosity.txt
grep 'pop1' CHR_10__10000001__10050000__stats.txt > ../window_pi.txt
grep 'pop1' CHR_10__10000001__10050000__stats.txt > ../window_Dxy.txt
grep 'pop1' CHR_10__10000001__10050000__stats.txt > ../window_Fst.txt

# add the relevant stats to each file
for i in $( ls *txt ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done
for i in $( ls *txt ); do grep 'pi' $i >> ../window_pi.txt; done
for i in $( ls *txt ); do grep 'Dxy' $i >> ../window_Dxy.txt; done
for i in $( ls *txt ); do grep 'Fst' $i >> ../window_Fst.txt; done



