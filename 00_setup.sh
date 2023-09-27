# define main working directory
workdir=/lustre/scratch/jmanthey/15_contopus

# make output directories
cd ${workdir}

mkdir 00_fastq
mkdir 01_cleaned
mkdir 01_bam_files
mkdir 02_vcf
mkdir 03_vcf
mkdir 04_filtered_vcf
mkdir 05_structure
mkdir 06_window_stats
mkdir 06_window_stats/windows
mkdir 07_find_fixed_differences
mkdir 08_ldhat
mkdir 08_ldhat/windows
mkdir 09_demography
mkdir 10_filter
mkdir 12_filter_vcf
mkdir 14_50kbp_admixture
mkdir 14_50kbp_admixture/windows
mkdir 20_satsuma
