## Function to report if a single fragment matches any seqences in the set with a set number of mismatches

match_probes <- function(fragment, sequence_set, mismatch_num) {
  counts <- c()
  # tic()
  counts <- vcountPattern(fragment, sequence_set, max.mismatch = mismatch_num) + 
    vcountPattern(reverseComplement(fragment), sequence_set, max.mismatch = mismatch_num)
  # toc()
  tf <- counts != 0
  if (sum(tf) == 0) {probe <- 1} else if (sum(tf) > 0) {probe <- 0} else probe <- 2
  return(probe)
}
