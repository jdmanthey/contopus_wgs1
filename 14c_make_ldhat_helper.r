	options(scipen=999)
	
	output_name <- "ldhat_helper.txt"
	
	# read in reference index
	# filtered to only include genotyped chromosomes
	ref_index <- read.table("flycatcher_rearranged.fa.fai", stringsAsFactors=F)
	# don't include W, 31, 32, 33 (all small and not included in other stats)
	ref_index <- ref_index[ref_index[,1] != "CHR_31" & ref_index[,1] != "CHR_32" & ref_index[,1] != "CHR_33" & ref_index[,1] != "CHR_W",]
	
	# define window size
	window_size <- 5000000
	
	# overlap of windows size
	overlap <- 50000
	
	# define intervals and write to helper file
	helper1 <- list() # job number
	helper2 <- list() # chromosome 
	helper3 <- list() # start 
	helper4 <- list() # end
	counter <- 1
	for(a in 1:nrow(ref_index)) {
		a_start <- 1
		a_end <- window_size
		a_max <- ref_index[a,2]
		a_chromosome <- ref_index[a,1]
		
		# determine number of windows for this chromosome
		if(window_size >= a_max) {
			a_windows <- 1
		} else {
			end_test <- a_end
			a_windows <- 1
			while(end_test < ref_index[a,2]) {
				a_windows <- a_windows + 1
				end_test <- end_test + window_size - overlap
			}
		}
		
		# loop for defining helper info for each window
		for(b in 1:a_windows) {
			if(b == a_windows & b != 1) {
				a_start <- a_start + window_size - overlap
				a_end <- a_max
				
				helper1[[counter]] <- counter
				helper2[[counter]] <- a_chromosome
				helper3[[counter]] <- a_start
				helper4[[counter]] <- a_end
				
				counter <- counter + 1
			} else if(b == 1 & b != a_windows) {
				a_start <- 1
				a_end <- window_size
				
				helper1[[counter]] <- counter
				helper2[[counter]] <- a_chromosome
				helper3[[counter]] <- a_start
				helper4[[counter]] <- a_end
				
				counter <- counter + 1
			} else if(b == 1 & b == a_windows) {
				a_start <- 1
				a_end <- a_max
				
				helper1[[counter]] <- counter
				helper2[[counter]] <- a_chromosome
				helper3[[counter]] <- a_start
				helper4[[counter]] <- a_end
				
				counter <- counter + 1
			} else {
				a_start <- a_start + window_size - overlap
				a_end <- a_start + window_size - 1
				
				helper1[[counter]] <- counter
				helper2[[counter]] <- a_chromosome
				helper3[[counter]] <- a_start
				helper4[[counter]] <- a_end
				
				counter <- counter + 1
			}
		}
	}
		
		output <- data.frame(job=as.numeric(unlist(helper1)), chr=as.character(unlist(helper2)), start=as.numeric(unlist(helper3)), end=as.numeric(unlist(helper4)))
		
		write.table(output, file=output_name, row.names=F, col.names=F, quote=F, sep="\t")
		
		
