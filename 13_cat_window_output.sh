cd /lustre/scratch/jmanthey/15_contopus/06_window_stats/windows

# combine the output for different analyses into a single file each
# first add a header for each file
grep 'pop1' CHR_10__10000001__10100000__stats.txt > ../window_heterozygosity.txt
grep 'pop1' CHR_10__10000001__10100000__stats.txt > ../window_theta.txt
grep 'pop1' CHR_10__10000001__10100000__stats.txt > ../window_pi.txt
grep 'pop1' CHR_10__10000001__10100000__stats.txt > ../window_Tajima_D.txt
grep 'pop1' CHR_10__10000001__10100000__stats.txt > ../window_Dxy.txt
grep 'pop1' CHR_10__10000001__10100000__stats.txt > ../window_Fst.txt

# add the relevant stats to each file
for i in $( ls *txt ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done
for i in $( ls *txt ); do grep 'theta' $i >> ../window_theta.txt; done
for i in $( ls *txt ); do grep 'pi' $i >> ../window_pi.txt; done
for i in $( ls *txt ); do grep 'Tajima_D' $i >> ../window_Tajima_D.txt; done
for i in $( ls *txt ); do grep 'Dxy' $i >> ../window_Dxy.txt; done
for i in $( ls *txt ); do grep 'Fst' $i >> ../window_Fst.txt; done



