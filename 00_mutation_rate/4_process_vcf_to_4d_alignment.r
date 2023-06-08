options(scipen=999)

# read in gff
gff <- read.table("manakin_cds_subset.gff", sep="\t", stringsAsFactors=F)

# read in vcf
vcf <- read.table("genes_filtered.simple.vcf", sep="\t", stringsAsFactors=F)

# number of individuals 
n_inds <- 3

# out_directory
temp_dir <- "temp_fasta"
dir.create(temp_dir)

# subset the gff to those chromosomes that were genotyped
gff <- gff[gff[,1] %in% unique(vcf[,1]), ]

# header
vcf_header <- c("CHROM", "POS", "REF", "ALT", "Contopus_KU31919",	"Neopipo_SRR9946662",	 "Onych_SRR9947257")
colnames(vcf) <- vcf_header
reference_name <- "Chiroxiphia"

# keep only invariant and biallelic sites
vcf <- vcf[nchar(vcf[,4]) == 1, ]

# define unique cds in gff file
cds <- unique(gff[,9])

###############################
###############################
###############################
###############################
###############################
# loop for each unique cds to determine whether to keep + write fasta
###############################
###############################
###############################
###############################
###############################

for(a in 1:length(cds)) {
	if(a %% 100 == 0) {print(a)} # progress
	# subset gff
	a_gff <- gff[gff[,9] == cds[a],]
	# keep going if the CDS is divisible by 3 (codon length and full)
	if(sum(a_gff[,5] - a_gff[,4] + 1) %% 3 == 0) {
		# filter vcf for this cds
		a_vcf <- list()
		for(b in 1:nrow(a_gff)) {
			a_vcf[[b]] <- vcf[vcf[,1] == a_gff[b,1] & vcf[,2] >= a_gff[b,4] & vcf[,2] <= a_gff[b,5],]
		}
		a_vcf <- do.call("rbind", a_vcf)
		
		# check if the nrow of this subset is at least 95% the total length of the cds
		if(nrow(a_vcf) >= (0.95 * sum(a_gff[,5] - a_gff[,4] + 1))) {
			# determine which sites are missing and fill them in with N
			cds_positions <- c()
			for(b in 1:nrow(a_gff)) {
				cds_positions <- c(cds_positions, seq(from=a_gff[b,4], to=a_gff[b,5], by=1))
			}
			cds_missing <- cds_positions[cds_positions %in% a_vcf[,2] == FALSE]
			if(length(cds_missing) > 0) {
				for(b in 1:length(cds_missing)) {
					a_vcf <- rbind(a_vcf, c(a_vcf[1,1], cds_missing[b], "N", ".", rep("0/0", n_inds)))
				}
			}
			
			# order the matrix in the correct order based on orientation of the cds (and rename rownames in order)
			a_vcf[,2] <- as.numeric(a_vcf[,2])
			if(a_gff[1,7] == "+") {
				a_vcf <- a_vcf[order(a_vcf[,2], decreasing=F),]
			} else {
				a_vcf <- a_vcf[order(a_vcf[,2], decreasing=T),]
			}
			rownames(a_vcf) <- seq(from=1, to=nrow(a_vcf), by=1)
			
			# randomly select alleles for heterozygotes
			for(b in 5:ncol(a_vcf)) {
				b_hets <- sample(c("0/0", "1/1"), length(a_vcf[a_vcf[,b] == "0/1",b]), replace=T)
				a_vcf[a_vcf[,b] == "0/1",b] <- b_hets
			}
			
			# determine alleles for each site
			allele_1 <- a_vcf[,3]
			allele_2 <- a_vcf[,4]
			# complement the alleles if the cds is in reverse orientation
			if(a_gff[1,7] == "-") {
				allele_1 <- gsub("A", 1, allele_1)
				allele_1 <- gsub("C", 2, allele_1)
				allele_1 <- gsub("G", 3, allele_1)
				allele_1 <- gsub("T", 4, allele_1)
				allele_1 <- gsub(1, "T", allele_1)
				allele_1 <- gsub(2, "G", allele_1)
				allele_1 <- gsub(3, "C", allele_1)
				allele_1 <- gsub(4, "A", allele_1)
				allele_2 <- gsub("A", 1, allele_2)
				allele_2 <- gsub("C", 2, allele_2)
				allele_2 <- gsub("G", 3, allele_2)
				allele_2 <- gsub("T", 4, allele_2)
				allele_2 <- gsub(1, "T", allele_2)
				allele_2 <- gsub(2, "G", allele_2)
				allele_2 <- gsub(3, "C", allele_2)
				allele_2 <- gsub(4, "A", allele_2)
			}
			
			# get the sequences for each individual
			a_genotypes <- a_vcf[,5:ncol(a_vcf)]
			seq_names <- list()
			seqs <- list()
			for(b in 1:n_inds) {
				b_rep <- a_genotypes[,b]
				b_rep[b_rep == "0/0"] <- allele_1[b_rep == "0/0"]
				b_rep[b_rep == "1/1"] <- allele_2[b_rep == "1/1"]
				seqs[[b]] <- paste0(b_rep, collapse="")
				seq_names[[b]] <- paste0(">", colnames(a_genotypes)[b])
			}
			# add the reference genome info
			b <- b + 1
			seqs[[b]] <- paste0(allele_1, collapse="")
			seq_names[[b]] <- paste0(">", reference_name)
			
			# write the fasta file to the temporary directory
			out_name <- paste0(temp_dir, "/", a, ".fasta")
			for(b in 1:(n_inds + 1)) {
				if(b == 1) {
					write(seq_names[[b]], out_name)
				} else {
					write(seq_names[[b]], out_name, append=T)
				}
				write(seqs[[b]], out_name, append=T)
			}		
		}
	}  
}

###############################
###############################
###############################
###############################
###############################
# function to identify 4 fold degenerate sites
###############################
###############################
###############################
###############################
###############################

codon_table <- read.table("codon_table.txt", header=T, stringsAsFactors=F)
determine_4d <- function(xxxx) {
	# remove codons with N
	if(length(grep("N", xxxx)) == 0) {
		# check if matches have 4D sites and return them
		if(codon_table[match(xxxx[1], codon_table[,1]),3] == "no") {
			return("")
		} else if(codon_table[match(xxxx[1], codon_table[,1]),2] == codon_table[match(xxxx[2], codon_table[,1]),2] & codon_table[match(xxxx[1], codon_table[,1]),2] == codon_table[match(xxxx[3], codon_table[,1]),2] & codon_table[match(xxxx[1], codon_table[,1]),2] == codon_table[match(xxxx[4], codon_table[,1]),2]) {
			return_object <- substr(xxxx, 3, 3)
		} else {
			return("")
		}
	}
}

###############################
###############################
###############################
###############################
###############################
# get 4-fold degenerate sites and make an alignment
###############################
###############################
###############################
###############################
###############################
library(seqinr)

# get fasta list
fasta_list <- list.files(path=temp_dir, full.names=T)

# loop for each fasta file
locus_results <- list()
for(a in 1:length(fasta_list)) {
	a_rep <- read.fasta(fasta_list[a])
	
	# determine how many codons and put into a list for each codon
	codons <- length(a_rep$Chiroxiphia) / 3
	
	# make codon list
	codon_list <- list()
	for(b in 1:codons) {
		codon_list[[b]] <- c(paste(toupper(a_rep$Contopus_KU31919[(b*3-2):(b*3)]), collapse=""),
							paste(toupper(a_rep$Neopipo_SRR9946662[(b*3-2):(b*3)]), collapse=""),
							paste(toupper(a_rep$Onych_SRR9947257[(b*3-2):(b*3)]), collapse=""),
							paste(toupper(a_rep$Chiroxiphia[(b*3-2):(b*3)]), collapse=""))
	}
	
	# use the determine_4d function to return four fold degenerate sites for this locus
	fourd_sites <- list()
	for(b in 1:length(codon_list)) {
		fourd_sites[[b]] <- determine_4d(codon_list[[b]])
	}
	
	# remove null results
	fourd_sites2 <- list()
	for(b in 1:length(fourd_sites)) {
		if(b == 1) { b_count <- 1 }
		 if(length(fourd_sites[[b]]) > 1) {
			fourd_sites2[[b_count]] <- fourd_sites[[b]]
			b_count <- b_count + 1
		}
	}
	
	# add 4d sites to total list
	locus_results[[a]] <- c(paste0(sapply(fourd_sites2, "[[", 1), collapse=""),
							paste0(sapply(fourd_sites2, "[[", 2), collapse=""),
							paste0(sapply(fourd_sites2, "[[", 3), collapse=""),
							paste0(sapply(fourd_sites2, "[[", 4), collapse=""))
	
	# print progress
	if(a %% 100 == 0) {
		print(floor(a / length(fasta_list) * 100))
	}
}

# concatenate results across loci
seqs <- list()
seq_names <- seq_names
for(a in 1:length(seq_names)) {
	seqs[[a]] <- paste0(sapply(locus_results, "[[", a), collapse="")
}
# total of X sites

output_name <- "_total_4d_sites.fasta"
for(a in 1:length(seq_names)) {
	if(a == 1) {
		write(seq_names[[a]], output_name)
	} else {
		write(seq_names[[a]], output_name, append=T)
	}
	write(seqs[[a]], output_name, append=T)
}














