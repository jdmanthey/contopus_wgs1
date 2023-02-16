options(scipen=999)

library(Biostrings)
genome <- readDNAStringSet("pseudochromosomes.fasta")

# filter scaffolds >= 500,000 bp
keep <- genome@ranges@width >= 500000
genome <- genome[keep]

# keep first 35 chromosomes
genome <- genome[1:35]

# rename the chromosomes to simpler names
genome_names <- genome@ranges@NAMES
genome_names <- paste0("CHR_", sapply(strsplit(sapply(strsplit(genome_names, "chromosome_"), "[[", 2), ","), "[[", 1))
genome@ranges@NAMES <- genome_names

# write output
writeXStringSet(genome, file="flycatcher_rearranged.fa")
