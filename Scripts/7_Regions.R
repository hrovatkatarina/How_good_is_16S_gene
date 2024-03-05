library(Biostrings)
library(stringr)
library(dplyr)

###############################################################################
# Get start and stop positions of variable regions based on primer pairs
# and make separate fasta files of all different V-regions for each genus
###############################################################################

# Import Multiple Sequence Alignment

file_name <- "full16S_MSAall.fasta"

fasta_file <- file_name
fasta_text <- readLines(fasta_file)
sub_fasta_text <- gsub("U", "T", fasta_text)
sub_fasta_file <- tempfile() 
writeLines(sub_fasta_text, sub_fasta_file)

dna_seqs <- readDNAStringSet(sub_fasta_file)

V1F <- DNAString("AGAGTTTGATCMTGGCTCAG")
V3R <- reverseComplement(DNAString("TTACCGCGGCKGCTGGCACG"))
V3F <- DNAString("CCTACGGGNGGCWGCAG")
V4R <- reverseComplement(DNAString("GGACTACNVGGGTWTCTAAT"))
V4F <- DNAString("GTGYCAGCMGCCGCGGTAA")
V5R <- reverseComplement(DNAString("CCGYCAATTYMTTTRAGTTT"))
V6F <- DNAString("AATTGACGGGGRCCCGC")
V8R <- reverseComplement(DNAString("ACGGGCRGTGWGTRCAA"))
V9R <- reverseComplement(DNAString("TACGGYTACCTTGTTAYGACTT"))


extract_and_write <- function(primer_f, primer_r, prefix) { 
  seq_f <- c()
  subseqs <- lapply(seq_along(dna_seqs), function(i) {
    position_f <- matchPattern(primer_f, dna_seqs[[i]], fixed = FALSE, max.mismatch = 3)
    position_r <- matchPattern(primer_r, dna_seqs[[i]], fixed = FALSE, max.mismatch = 3)
    
    if (length(position_f) > 0 && length(position_r) > 0) { 
      seq_f <<- c(seq_f, i)
      start <- start(position_f[1])
      end <- start(position_r[1]) + length(primer_r)
      # Extract the subsequence
      subseq(dna_seqs[[i]], start = start, end = end)
    } else {
      NULL
    }
  })
  
  subseqs <- Filter(Negate(is.null), subseqs)
  sequences <- DNAStringSet(subseqs)
  names(sequences) <- names(dna_seqs)[seq_f]
  
  # Split the sequences by genus
  sequences_by_genus <- split(sequences, sapply(strsplit(names(sequences), " "), `[`, 2))
  
  # Write each group of sequences to a separate file
  lapply(names(sequences_by_genus), function(genus) {
    writeXStringSet(sequences_by_genus[[genus]], format = "fasta", 
                    filepath = paste0("Data\\",genus,"\\",genus, "_", prefix, ".fasta"))
  })
}

extract_and_write(V1F, V3R, "V1V3")
extract_and_write(V3F, V4R, "V3V4")
extract_and_write(V4F, V4R, "V4")
extract_and_write(V4F, V5R, "V4V5")
extract_and_write(V6F, V8R, "V6V8")
extract_and_write(V6F, V9R, "V6V9")

# Write full 16S sequences by genus

seq_by_genus <- split(dna_seqs, sapply(strsplit(names(dna_seqs), " "), `[`, 2))

lapply(names(seq_by_genus), function(genus) {
    writeXStringSet(seq_by_genus[[genus]], format = "fasta", filepath = paste0("Data\\",genus,"\\",genus,"_full16S.fasta"))
})