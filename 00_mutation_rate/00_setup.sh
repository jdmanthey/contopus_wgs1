# define main working directory
workdir=/lustre/scratch/jmanthey/15_contopus/30_mutation_rate

# make output directories
cd ${workdir}

mkdir 00_fastq
mkdir 01_cleaned
mkdir 01_bam_files
mkdir 02_vcf
mkdir 03_vcf
mkdir 04_filtered_vcf

module load intel java bwa samtools singularity

# index reference
bwa index GCF_009829145.1_bChiLan1.pri_genomic.fna

samtools faidx GCF_009829145.1_bChiLan1.pri_genomic.fna

cd 

java -jar picard.jar CreateSequenceDictionary R=/lustre/scratch/jmanthey/15_contopus/30_mutation_rate/GCF_009829145.1_bChiLan1.pri_genomic.fna O=/lustre/scratch/jmanthey/15_contopus/30_mutation_rate/GCF_009829145.1_bChiLan1.pri_genomic.dict

