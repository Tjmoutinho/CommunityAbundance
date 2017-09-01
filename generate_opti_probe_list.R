generate_opti_probe_list <- function(genomes){
  
  percent_vec <- numeric()
  count_vec <- numeric()
  columns = c("C|G", "N", "Index", "Contig", "lig", "1st_half", "2nd_half", "tf", "Freq_1st_half", "ASF")
  kmers_optimal_allasf <- data.frame(matrix(vector(), 0, length(columns),
                                     dimnames=list(c(), columns)),
                                     stringsAsFactors=F)
  
  for(j in 1:length(genomes)){
    genome <- genomes[j]
    
    ### Read in Genome ###
    ASF_member <- readDNAStringSet(genome)
    
    ### Chop Genome up into 40 bp kmers ###
    view.width <- 40
    letters <- c("CG","N")
    columns = c("C|G", "N", "Index", "Contig")
    kmers_40bp_ID_C_all <- data.frame(matrix(vector(), 0, 4,
                                      dimnames=list(c(), columns)),
                                      stringsAsFactors=F)
    n <- length(ASF_member)
    for (i in 1:n){
      kmers_40bp <- letterFrequencyInSlidingView(ASF_member[[i]], view.width, letters, as.prob=T)
      a <- seq_len(nrow(kmers_40bp))
      kmers_40bp_ID <- cbind(kmers_40bp, Index=a)
      kmers_40bp_ID_C <- cbind(kmers_40bp_ID, Contig=i)
      kmers_40bp_ID_C_all <- rbind(kmers_40bp_ID_C_all, kmers_40bp_ID_C)
      rm(kmers_40bp, a, kmers_40bp_ID)
    }
    
    ### Delete Kmers comtaining at least one 'N'###
    kmers_40bp_noN <- filter(kmers_40bp_ID_C_all, N == 0)
    
    ### Isolate probes with -AA- in middle for optimal ligation ###
    n = length(ASF_member)
    lig_site <- DNAStringSet()
    kmer20a   <- DNAStringSet()
    kmer20b   <- DNAStringSet()
    for (i in 1:n){
      kmers_40bp_noN_ci <- filter(kmers_40bp_noN, Contig == i)
      range <- IRanges(start = kmers_40bp_noN_ci$Index+19, end = kmers_40bp_noN_ci$Index+20)
      range2 <- IRanges(start = kmers_40bp_noN_ci$Index, end = kmers_40bp_noN_ci$Index+19)
      range3 <- IRanges(start = kmers_40bp_noN_ci$Index+20, end = kmers_40bp_noN_ci$Index+39)
      lig_site <- append(lig_site, extractAt(ASF_member[[i]], range))
      kmer20a <- append(kmer20a, extractAt(ASF_member[[i]], range2))
      kmer20b <- append(kmer20b, extractAt(ASF_member[[i]], range3))
    }
    kmers <- cbind(kmers_40bp_noN, lig = lig_site, `1st_half` = kmer20a, `2nd_half` = kmer20b)
    
    ### Filter kmers based on optimal ligation probe ###
    tf <- kmers$lig == 'AA'
    kmers_lig_tf <- cbind(kmers, tf)
    
    ### Filter kmers based on GC content being identical in 40bp kmer and 20bp kmers. ###
    letters <- c("CG")
    freq <- letterFrequency(kmer20a, letters, as.prob=TRUE)
    colnames(freq) <- c("Freq_1st_half")
    kmers_lig_tf_freq <- cbind(kmers_lig_tf, freq)
    
    # Apply halves GC frequency and ligation probe center filters
    kmers_opti_lig <- filter(kmers_lig_tf_freq, `tf`== TRUE)
    kmers_optimal <- filter(kmers_opti_lig, kmers_opti_lig$Freq_1st_half == kmers_opti_lig$`C|G`)
    asf_name_vec <- rep(genome, nrow(kmers_optimal))
    kmers_optimal_asf <- cbind(kmers_optimal, ASF = asf_name_vec)
    
    ### Return two values of interest: the percent with the max probes and how many probes there are in that group. ###
    kmers_gc_range <- filter(kmers_optimal_asf, `C|G` >= 0.35 & `C|G` <= 0.65)
    freq <- as.data.frame(table(unlist(kmers_gc_range$`C|G`)))
    gc_opti <- filter(freq, freq$Freq == max(freq$Freq))

    if (j == 1) {
      count_matrix <- matrix(ncol=length(genomes),nrow=nrow(freq))
      rownames(count_matrix) <- freq$Var1
      colnames(count_matrix) <- genomes
    }
    ### Create data frames of important info to keep with each revolution in loop ###
    count_matrix[,j] <- c(freq$Freq)
    kmers_optimal_allasf <- rbind(kmers_optimal_allasf, kmers_gc_range)
    sprintf("Genome %g: done", j)
  }
  
  all_probe_options <- count_matrix
  min_vec <- apply(count_matrix,1,min)
  Opti_GC <- as.numeric(levels(gc_opti$Var1))[which.max(min_vec)]
  c <- filter(kmers_optimal_allasf, `C|G` == Opti_GC)
  
  ### Add global alignment scores to dataframe to limit the chance of hairpin formation and homo-dimerization
  score <- pairwiseAlignment(DNAStringSet(c$`1st_half`), reverseComplement(DNAStringSet(c$`2nd_half`)), type = "global", scoreOnly = TRUE)
  kmers_optimal_allasf_opti <- cbind(c, score)
  
  return(list(Opti_GC, all_probe_options, kmers_optimal_allasf_opti))
}

###
# https://www.idtdna.com/pages/decoded/decoded-articles/ask-alex/ask-alex/2012/04/11/how-can-i-determine-if-a-hairpin-in-my-oligo-is-too-strong-to-allow-hybridization-
# Using free OligoAnalyzer® software, part of the IDT SciTools programs, enter your oligonucleotide sequence
# and choose “Hairpin.” The software will generate a series of possible hairpin structures. You can arrange 
# these structures in order of decreasing melting temperature (Tm). If the highest hairpin Tm is at or above 
# your annealing temperature, that hairpin is likely to impede hybridization.
###

### Create histogram of kmer GC content ###
# ggplot(data=kmers_opti_lig, aes(kmers_opti_lig$`C|G`)) + geom_histogram(breaks=seq(0, 1, by = 0.005))
# mean(kmers_opti_lig$`C|G`)
# sd(kmers_opti_lig$`C|G`)
# unique(kmers_opti_lig$`C|G`)

### Make a bunch of variables ###
# x <- as.list(rnorm(10))
# names(x) <- paste("a", 1:length(x), sep = "")
# list2env(x , envir = .GlobalEnv)
