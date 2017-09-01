# Call Genomes to Identify Ideal Probe Sites
# Thomas J Moutinho Jr, July 2017

### QUESTIONS FOR PAUL ###
# Is there a better way to clean the RAM that was previously used after removing a variable? I've used gc(), 
#   but that doesn't seem to get everything. 

# Are there any noticable parts of this that could be done faster? 

# When I align the optimal probes to all of the genomes, I don't remove the regions of the genomes that have 
#   N's (instead of A,T,G, or C). I reason that it is better to keep those regions because when we allow mismatches
#   the probes will be able to match in those regions if there is sufficient alignment. What do you think? 

### Clear Global Environment ###
rm(list=ls())
gc()

### Load Libraries ###
library(Biostrings) #Need Bioconductor Package
library(dplyr)
library(ggplot2)
library(tictoc)
library(parallel)
library(primer3)

### Load Functions ###
source("generate_opti_probe_list.R") ## Edit the name of the column for the genome file names.
source("return_multi_vars.R")
source("match_probes.R")

### Run function to find the optimal probes ###
genome_fastas <- c("ASF356.fna", "ASF360.fna", "ASF361.fna", "ASF457.fna", 
                   "ASF492.fna", "ASF500.fna", "ASF502.fna", "ASF519.fna")
# tic()
c(opti_percent_gc, probe_counts, opti_probes) := generate_opti_probe_list(genome_fastas) # 21 seconds per genome
# toc()

### Get top x number of probes for each member ###
top_num = 11 ## This number sets how many probes for each species are tested for alignment (with mismatches) to the rest of the genomes.

columns = c()
n = length(top_num*length(genome_fastas))
top_of_each <- data.frame(matrix(vector(), n, 0,
                           dimnames=list(c(), columns)),
                           stringsAsFactors=F)
for(i in 1:length(genome_fastas)){
  df <- head(opti_probes,top_num)
  one_mem <- filter(opti_probes, opti_probes$ASF == genome_fastas[i])
  one_mem_order <- one_mem[order(one_mem$score, decreasing = FALSE),] %>% head(top_num)
  top_of_each <- rbind(top_of_each, one_mem_order)
}

### Parallelized code to check how many mismatches each probe can have without aligning to genomic DNA ###
mismatch_num = 10 ## 25% of the probe can be mismatched. 
probe_or_no_set <- vector(mode = "numeric", nrow(top_of_each));
unique_to_self_set <- vector(mode = "numeric", nrow(top_of_each));

for(k in 1:length(genome_fastas)){
  tic()
  ### Make a DNAStringSet of the fragments to feed into vcountPattern ###
  fragments <- DNAStringSet()
  range <- (1+(top_num*(k-1))):((top_num*(k-1))+top_num)
  for(n in range){
    i <- match(top_of_each$ASF[n], genome_fastas)
    j <- top_of_each$Contig[n]
    range <- IRanges(start = top_of_each$Index[n], end = top_of_each$Index[n]+39)
    genomes <- sapply(genome_fastas, readDNAStringSet)
    fragment <- extractAt(genomes[[i]][[j]], range)
    fragments <- append(fragments, fragment)
  }
  ### Check for uniqueness in own genome ###
  genome <- readDNAStringSet(genome_fastas[k])
  n = length(genome)
  genome_seq_set <- DNAStringSet()  
  for (i in 1:n){
    df_contig_member <- filter(top_of_each, top_of_each$ASF == genome_fastas[k], top_of_each$Contig == i)
    df_contig_member_ordered <- df_contig_member[order(df_contig_member$Index),]
    set_of_index_vals <- df_contig_member_ordered$Index
    if(length(set_of_index_vals) > 0){
      index_ranges <- vector(mode = "numeric")
      for(w in 1:length(set_of_index_vals)){
        index_set <- (set_of_index_vals[w]-mismatch_num):(set_of_index_vals[w]+39+mismatch_num)
        index_ranges <- append(index_ranges, index_set)
      }
      index_full <- 1:(length(genome[[i]])-39)
      index <- index_full[-c(index_ranges)]
      range  <- IRanges(start = index, end = index+39)
      genome_seq_set <- append(genome_seq_set, extractAt(genome[[i]], range))
    }
    else if(length(set_of_index_vals) == 0){
      index <- 1:(length(genome[[i]])-39)
      range  <- IRanges(start = index, end = index+39)
      genome_seq_set <- append(genome_seq_set, extractAt(genome[[i]], range))
    }
  }
  sequence_set_one <- genome_seq_set # This set of sequences is from one member of the set, but without the sequences that are where the optimal probes will match to. 
  
  ### Create a genome set to assay the probes against ###
  genomes <- sapply(genome_fastas[-k], readDNAStringSet)
  n = length(genomes)
  genomic_seq_set <- DNAStringSet()  
  for (i in 1:n){
    m = length(genomes[[i]])
    for (j in 1:m){
      index <- 1:(length(genomes[[i]][[j]])-39)
      range  <- IRanges(start = index, end = index+39)
      genomic_seq_set <- append(genomic_seq_set, extractAt(genomes[[i]][[j]], range))
    }
  }
  sequence_set <- genomic_seq_set # This set of sequences does not contain data from one member of the set. 
  
  ### Parallel Code ###
  no_cores <- 4 #detectCores()/4
  clust <- makeCluster(no_cores)
  eval <- clusterEvalQ(clust, {source("match_probes.R")
                               library(Biostrings)
                               library(tictoc)})
  clusterExport(clust, c("sequence_set", "sequence_set_one","mismatch_num"))
  range <- (1+(top_num*(k-1))):((top_num*(k-1))+top_num)
  probe_or_no_set[range] <- parSapply(clust, fragments, function(x) match_probes(x, sequence_set, mismatch_num))
  unique_to_self_set[range] <- parSapply(clust, fragments, function(x) match_probes(x, sequence_set_one, mismatch_num))
  stopCluster(clust)
  toc()
  gc()
}
top_of_each_mm_uts <- cbind(top_of_each, probe_or_no_set, unique_to_self_set)

probe_options_all_data <- filter(top_of_each_mm_uts, top_of_each_mm_uts$probe_or_no_set == 1, top_of_each_mm_uts$unique_to_self_set == 1)
PROBE_OPTIONS <- probe_options_all_data[ , c("C|G","Index","Contig","1st_half","2nd_half","ASF")]

### Print results to Excel file ###
library(xlsx)
write.xlsx2(PROBE_OPTIONS, "PROBE_OPTIONS.xlsx")

