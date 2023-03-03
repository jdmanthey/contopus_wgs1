
	# popmap = individual base names of fastq files, one line per individual
	# make sure reference is indexed with bwa and samtools before use, and use CreateSequenceDictionary in GATK 
	options(scipen=999)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################	
##############    START OPTIONS     #####################################

	# project specific options
	project_directory <- "/lustre/scratch/jmanthey/19_mexico/catharus"
	directory_name <- "11_genotype_catharus"
	reference_genome_location <- "/home/jmanthey/references/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna"
	cluster <- "quanah"
	output_name <- "catharus"
	popmap <- "popmap_catharus.txt"
	individuals <- read.table(popmap, sep="\t")
	individuals[,1] <- as.character(individuals[,1])
	faidx <- read.table("catharus.fai", stringsAsFactors=F)
	singularity_cache <- "/lustre/work/jmanthey/singularity-cachedir"
	name_of_gatk_singularity_image <- "gatk_4.2.3.0.sif"

	# define minimum and maximum genotyping job sizes
	min_scaffold_size <- 2000000
	max_genotype_job_size <- 5000000
	max_individual_genotype_job_size <- 100000000
	
	# define number of cores for each step
	ncores_step1 <- 8
	ncores_step2 <- 8
	ncores_step3 <- 12
	# define memory per core for each step
	mem_step1 <- "8G"
	mem_step2 <- "8G"
	mem_step3 <- "8G"

##############    END OPTIONS     #######################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

	# calculate total java memory per job (= memory available * 0.9)
	total_mem1 <- ceiling(ncores_step1 * as.numeric(strsplit(mem_step1, "G")[[1]]) * 0.9)
	total_mem2 <- ceiling(ncores_step2 * as.numeric(strsplit(mem_step2, "G")[[1]]) * 0.9)
	total_mem3 <- ceiling(ncores_step3 * as.numeric(strsplit(mem_step3, "G")[[1]]) * 0.9)

	# make directories
	dir.create(directory_name)
	dir.create(paste0(directory_name, "/01_gatk_split"))
	dir.create(paste0(directory_name, "/02b_gatk_database"))
	dir.create(paste0(directory_name, "/03b_group_genotype_database"))

	# subset the index
	faidx <- faidx[faidx[,2] >= min_scaffold_size, ]

	# finds scaffolds too big to genotype at once
	faidx_keep <- faidx[faidx[,2] < max_individual_genotype_job_size,1:2]
	faidx_change <- faidx[faidx[,2] >= max_individual_genotype_job_size,1:2]
	# paste the interval to use to each of the faidx objects
	faidx_keep <- cbind(faidx_keep, faidx_keep[,1], rep(1, nrow(faidx_keep)), faidx_keep[,2])
	faidx_change <- cbind(faidx_change, rep("x", nrow(faidx_change)))
	new_faidx_change <- c()
	for(a in 1:nrow(faidx_change)) {
		a_rep <- faidx_change[a,]
		a_breaks <- floor(as.numeric(a_rep[1,2]) / 2)
		a_break1 <- c(a_rep[1,1], a_breaks, paste0(a_rep[1,1], ":1-", a_breaks), 1, a_breaks)
		a_break2 <- c(paste0(a_rep[1,1], "b"), as.numeric(a_rep[1,2]) - a_breaks, paste0(a_rep[1,1], ":", a_breaks + 1, "-", as.numeric(a_rep[1,2])), a_breaks + 1, as.numeric(a_rep[1,2]))
		new_faidx_change <- rbind(new_faidx_change, a_break1, a_break2)
	}
	colnames(faidx_keep) <- c("id", "length", "interval", "start", "end")
	colnames(new_faidx_change) <- c("id", "length", "interval", "start", "end")
	faidx <- rbind(faidx_keep, new_faidx_change)
	faidx[,3] <- as.character(faidx[,3])
	faidx[,1] <- as.character(faidx[,1])
	faidx[,2] <- as.numeric(faidx[,2])
	faidx[,4] <- as.numeric(faidx[,4])
	faidx[,5] <- as.numeric(faidx[,5])
	faidx <- na.omit(faidx)
	
	# write the helper files for the 1st genotyping step
	for(a in 1:nrow(individuals)) {
		if(a == 1) {
			helper1 <- cbind(faidx[,1], as.character(faidx[,3]), rep(individuals[a,1], nrow(faidx)))
		} else {
			helper1 <- rbind(helper1, cbind(faidx[,1], as.character(faidx[,3]), rep(individuals[a,1], nrow(faidx))))
		}
	}
	# write the chromosome/scaffold to genotype for each job
	write(helper1[,2], file=paste0(directory_name, "/01_gatk_split", "/helper1.txt"), ncolumns=1)
	# write the individual to genotype for each job
	write(helper1[,3], file=paste0(directory_name, "/01_gatk_split", "/helper2.txt"), ncolumns=1)
	# write the chromosome name for each job
	write(helper1[,1], file=paste0(directory_name, "/01_gatk_split", "/helper1b.txt"), ncolumns=1)
	
	# step 1
	# genotype all individuals using GATK, one array job using the above two helper files

	a.script <- paste0(directory_name, "/01_gatk_split/step1_array.sh")
	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste0("#SBATCH --job-name=", "step1"), file=a.script, append=T)
	write("#SBATCH --nodes=1", file=a.script, append=T)
	write(paste0("#SBATCH --ntasks=", ncores_step1), file=a.script, append=T)
	write(paste0("#SBATCH --partition=", cluster), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write(paste0("#SBATCH --mem-per-cpu=", mem_step1), file=a.script, append=T)
	write(paste0("#SBATCH --array=1-", nrow(helper1)), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load singularity", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("chr_array=$( head -n${SLURM_ARRAY_TASK_ID} helper1.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("ind_array=$( head -n${SLURM_ARRAY_TASK_ID} helper2.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("name_array=$( head -n${SLURM_ARRAY_TASK_ID} helper1b.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste0('export SINGULARITY_CACHEDIR="', singularity_cache, '"'), file=a.script, append=T)
	write("", file=a.script, append=T)
	
	#gatk 4.0
	a_name <- paste0(project_directory, "/01_bam_files/", "${ind_array}", "_final.bam")
	gatk_command <- paste0('singularity exec $SINGULARITY_CACHEDIR/', name_of_gatk_singularity_image, ' gatk --java-options "-Xmx', total_mem1, 'g" HaplotypeCaller -R ', reference_genome_location, " -I ", a_name, " -ERC GVCF -O ", project_directory, "/02_vcf/", "${name_array}", "._${ind_array}_.g.vcf", " --QUIET --intervals ", "${chr_array}")
	write(gatk_command, file=a.script, append=T)
	
	
	
	
	# write the helper files for the 2nd genotyping step
	# write the chromosome/scaffold to database for each job
	write(faidx[,1], file=paste0(directory_name, "/02b_gatk_database", "/helper3.txt"), ncolumns=1)
	write(faidx[,3], file=paste0(directory_name, "/02b_gatk_database", "/helper3b.txt"), ncolumns=1)

	
	
	# step 2
	# create genotyping database for each of the chromosomes
	a.script <- paste0(directory_name, "/02b_gatk_database/step2_array.sh")

	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste0("#SBATCH --job-name=", "step2"), file=a.script, append=T)
	write("#SBATCH --nodes=1", file=a.script, append=T)
	write(paste0("#SBATCH --ntasks=", ncores_step2), file=a.script, append=T)
	write(paste0("#SBATCH --partition=", cluster), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write(paste0("#SBATCH --mem-per-cpu=", mem_step2), file=a.script, append=T)
	write(paste0("#SBATCH --array=1-", nrow(faidx)), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load singularity", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("name_array=$( head -n${SLURM_ARRAY_TASK_ID} helper3.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("interval_array=$( head -n${SLURM_ARRAY_TASK_ID} helper3b.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste0('export SINGULARITY_CACHEDIR="', singularity_cache, '"'), file=a.script, append=T)
	write("", file=a.script, append=T)
		
	#make list of all vcfs to database
	for(b in 1:nrow(individuals)) {
		if(b == 1) {
			vcf_total <- paste0("-V ", project_directory, "/02_vcf/", "${name_array}", "._", individuals[b,1], "_.g.vcf")
		} else {
			vcf_total <- paste0(vcf_total, " -V ", project_directory, "/02_vcf/", "${name_array}", "._", individuals[b,1], "_.g.vcf")
		}
	}
	
	#gatk 4.0
	gatk_command <- paste0('singularity exec $SINGULARITY_CACHEDIR/', name_of_gatk_singularity_image, ' gatk --java-options "-Xmx', total_mem2, 'g" GenomicsDBImport --genomicsdb-shared-posixfs-optimizations ', vcf_total, " --genomicsdb-workspace-path ", project_directory, "/02_vcf/", "${name_array}", " -L ", "${interval_array}")
				
	write(gatk_command, file=a.script, append=T)
		




	# write the helper files for the 3rd genotyping step
	helper <- c()
	job_suffixes <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
	for(a in 1:nrow(faidx)) {
		a_rep <- faidx[a,]
		job_number <- ceiling(a_rep[1,2] / max_genotype_job_size)
		job_number <- job_suffixes[1:job_number]
		
		if(length(job_number) > 1) {
			# define interval to group genotype if more than one job
			interval_start <- a_rep[1,4]
			interval_end <- a_rep[1,4] + max_genotype_job_size - 1
			for(b in 1:length(job_number)) {
				helper <- rbind(helper, c(faidx[a,1], 
					paste0(strsplit(faidx[a,3], ":")[[1]][1], ":", interval_start, "-", interval_end),
					paste0(faidx[a,1], "__", job_number[b])
					))
					interval_start <- interval_start + max_genotype_job_size
				if((b+1) != length(job_number)) {
					interval_end <- interval_end + max_genotype_job_size
				} else {
					interval_end <- faidx[a,5]
				}
			}
		} else {
			helper <- rbind(helper, c(faidx[a,1], 
				paste0(faidx[a,1], ":", 1, "-", faidx[a,2]),
				paste0(faidx[a,1])
				))
		}
	}
	# write the chromosome/scaffold to genotype for each job
	write(helper[,1], file=paste0(directory_name, "/03b_group_genotype_database", "/helper4.txt"), ncolumns=1)
	# write the interval range to genotype for each job
	write(helper[,2], file=paste0(directory_name, "/03b_group_genotype_database", "/helper5.txt"), ncolumns=1)
	# write the output base name for each job
	write(helper[,3], file=paste0(directory_name, "/03b_group_genotype_database", "/helper6.txt"), ncolumns=1)
	
	
	


	

	# step 3
	# group genotyping for each interval
	a.script <- paste0(directory_name, "/03b_group_genotype_database/step3_array.sh")
	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste0("#SBATCH --job-name=", "step3"), file=a.script, append=T)
	write("#SBATCH --nodes=1", file=a.script, append=T)
	write(paste0("#SBATCH --ntasks=", ncores_step3), file=a.script, append=T)
	write(paste0("#SBATCH --partition=", cluster), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write(paste0("#SBATCH --mem-per-cpu=", mem_step3), file=a.script, append=T)
	write(paste0("#SBATCH --array=1-", nrow(helper)), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load singularity", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("chr_array=$( head -n${SLURM_ARRAY_TASK_ID} helper4.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("interval_array=$( head -n${SLURM_ARRAY_TASK_ID} helper5.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("name_array=$( head -n${SLURM_ARRAY_TASK_ID} helper6.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste0('export SINGULARITY_CACHEDIR="', singularity_cache, '"'), file=a.script, append=T)
	write("", file=a.script, append=T)
	
	#gatk 4.0
	gatk_command <- paste0('singularity exec $SINGULARITY_CACHEDIR/', name_of_gatk_singularity_image, ' gatk --java-options "-Xmx', total_mem3, 'g" GenotypeGVCFs --genomicsdb-shared-posixfs-optimizations -R ', reference_genome_location, " -V gendb://", project_directory, "/02_vcf/", "${chr_array}", " --include-non-variant-sites -O ", project_directory, "/03_vcf/", "${name_array}", ".g.vcf", " -L ", "${interval_array}")
	write(gatk_command, file=a.script, append=T)