options(scipen=999)

# this summarizes the output as estimates of 4Nr per bp

# window size for summarization
window_size <- 50000

# read in the helper file
helper <- read.table("ldhat_helper.txt")

# how much overlap on files?
# find from the helper file
file_overlap <- helper[1,4] - helper[2,3] + 1

# directory of files to read from
output_directory <- "output"

# prefix on output files
prefix <- "eastern__"

# unique chromosomes to loop through
uni_chr <- unique(helper[,2])

# loop for each chromosome
output <- list()
counter <- 1
for(a in 1:length(uni_chr)) {
	chr_values <- list()
	print(a)
	# how many output files for this chromosome?
	a_help <- helper[helper[,2] == uni_chr[a],]
	a_num_windows <- nrow(a_help)
	# summarize if only one file otherwise loop for each file
	if(a_num_windows == 1) {
			b_results <- read.table(paste0(output_directory, "/", prefix, a_help[1,2], "__", a_help[1,3], "__", a_help[1,4], ".res.txt"), header=T)
			b_mapping <- read.table(paste0(output_directory, "/", prefix, a_help[1,2], "__", a_help[1,3], "__", a_help[1,4], ".mapping"), header=F)
			
			# combine results and mapping
			b_results[,1] <- b_mapping[,1]
			
			# remove the first line (always super inflated)
			b_results <- b_results[2:nrow(b_results),]
			
			# keep only the locus position and median value for each site
			b_results <- b_results[,c(1,3)]
			
			chr_values[[1]] <- b_results
			
	} else {
		for(b in 1:a_num_windows) {
			# read in output and mapping file
			b_results <- read.table(paste0(output_directory, "/", prefix, a_help[b,2], "__", a_help[b,3], "__", a_help[b,4], ".res.txt"), header=T)
			b_mapping <- read.table(paste0(output_directory, "/", prefix, a_help[b,2], "__", a_help[b,3], "__", a_help[b,4], ".mapping"), header=F)
			
			# combine results and mapping
			b_results[,1] <- b_mapping[,1]
			
			# remove the first line (always super inflated)
			b_results <- b_results[2:nrow(b_results),]
			
			# keep only the locus position and median value for each site
			b_results <- b_results[,c(1,3)]
			
			# different situations for first, last, and middle output files to account for overlap
			if(b == 1) {
				b_results <- b_results[b_results[,1] <= (a_help[b,4] - (file_overlap/2)),]
				chr_values[[b]] <- b_results
			} else if(b == a_num_windows) {
				b_results <- b_results[b_results[,1] > (a_help[b,3] + (file_overlap/2)),]
				chr_values[[b]] <- b_results
			} else {
				b_results <- b_results[b_results[,1] <= (a_help[b,4] - (file_overlap/2)),]
				b_results <- b_results[b_results[,1] > (a_help[b,3] + (file_overlap/2)),]
				chr_values[[b]] <- b_results
			}
		}
	}
	
	# combine all window values
	chr_values <- do.call(rbind, chr_values)
	
	# summarize all output for this chromosome in windows
	n_chr_windows <- floor(max(a_help[,4]) / window_size)
	window_start <- 1
	window_end <- window_size
	for(b in 1:n_chr_windows) {
		b_chr_values <- chr_values[chr_values[,1] >= window_start & chr_values[,1] <= window_end,]
		if(nrow(b_chr_values) > 0) {
			b_chr_values <- mean(b_chr_values[,2]) / 1000
		} else {
			b_chr_values <- NA
		}
		output[[counter]] <- c(a_help[1,2], window_start, window_end, b_chr_values)
		
		window_start <- window_start + window_size
		window_end <- window_end + window_size
		counter <- counter + 1
	}

}

total_output <- do.call(rbind, output)
total_output <- data.frame(chr=as.character(total_output[,1]), start=as.numeric(total_output[,2]), end=as.numeric(total_output[,3]), rho=as.numeric(total_output[,4]))

plot(total_output$rho, pch=19, cex=0.2)

write.table(total_output, file=paste0(prefix, "rho_windows.txt"), row.names=F, quote=F)
























