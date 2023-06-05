options(scipen=999)

x_files <- list.files(pattern="*__b.2.Q")

output <- c("chr", "start", "end", "individual", "admix1", "admix2")

individuals <- scan("ingroup.txt", what="character")
ind_number <- length(individuals)

write(output, file="admixture_windows.txt", sep="\t", ncolumns=6)

for(a in 1:length(x_files)) {
	a_rep <- read.table(x_files[[a]])
	output_df <- data.frame(chr=as.character(rep(as.character(strsplit(x_files[a], "__")[[1]][1]), ind_number)),
							start=as.numeric(rep(as.character(strsplit(x_files[a], "__")[[1]][2]), ind_number)),
							end=as.numeric(rep(as.numeric(strsplit(x_files[a], "__")[[1]][3]), ind_number)),
							individuals=as.character(individuals),
							admix1=as.numeric(a_rep[,1]),
							admix2=as.numeric(a_rep[,2]))
	write.table(output_df, file="admixture_windows.txt", sep="\t", quote=F, row.names=F, col.names=F, append=T)

	
}

