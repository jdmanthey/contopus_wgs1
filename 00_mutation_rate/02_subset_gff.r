
x <- read.table("GCF_009829145.1_bChiLan1.pri_genomic.gff", stringsAsFactors=F, sep="\t", fill=T)

# keep only cds
x2 <- x[x[,3] == "CDS",]

# change row names
rownames(x2) <- seq(from=1, to=nrow(x2), by=1)

# get list of all genes
gene_list <- sapply(strsplit(x2[,9], ";"), "[[", 6)
uniq_gene_list <- unique(gene_list)

# keep only basic attribute data
x2[,9] <- sapply(strsplit(x2[,9], ";"), "[[", 1)

# keep only largest isoform per gene
keep <- list()
for(a in 1:length(uniq_gene_list)) {
	a_rep <- x2[gene_list == uniq_gene_list[a],]
	a_cds_uniq <- unique(a_rep[,9])
	if(length(a_cds_uniq) == 1) {
		# if only one isoform, keep it
		keep[[a]] <- a_rep[1,9]
	} else {
		# find counts for each isoform
		b_rep <- c()
		for(b in 1:length(a_cds_uniq)) {
			b_rep <- c(b_rep, length(a_rep[a_rep[,9] == a_cds_uniq[b], 9]))
		}
		# keep the larger isoform
		keep[[a]] <- a_cds_uniq[b_rep == max(b_rep)][1]
	}
}
keep <- unlist(keep)
x3 <- x2[x2[,9] %in% keep, ]

# write the cds to file
write.table(x3, file="manakin_cds_subset.gff", sep="\t", row.names=F, col.names=F, quote=F)


