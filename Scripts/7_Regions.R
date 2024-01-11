library(Biostrings)
library(stringr)
library(dplyr)

###############################################################################
# Get start and stop positions of variable regions based on primer pairs
###############################################################################

# Import Multiple Sequence Alignment

file_name <- "Data/full16S_MSAall.fasta"

fasta_file <- file_name
fasta_text <- readLines(fasta_file)
fixed_fasta_text <- gsub("U", "T", fasta_text)
fixed_fasta_file <- tempfile() 
writeLines(fixed_fasta_text, fixed_fasta_file)

dna_seqs <- readDNAStringSet(fixed_fasta_file)

extract_position <- function(primer, start_position, reverse_complement) {
  # If the primer needs to be reverse complement
  if (reverse_complement) {
    primer <- DNAString(primer)
    primer <- reverseComplement(primer)
  }
  
  position <- vmatchPattern(primer, dna_seqs, fixed = FALSE, max.mismatch = 3)
  start_position <- start(position[[1]])
  
  # If the primer is reverse_complement, add the length of the primer
  if (reverse_complement) {
    start_position <- start_position + length(primer)
  }
  
  return(start_position)
}

V1F <- extract_position("AGAGTTTGATCMTGGCTCAG", TRUE, FALSE)
V3R <- extract_position("TTACCGCGGCKGCTGGCACG", FALSE, TRUE)
V3F <- extract_position("CCTACGGGNGGCWGCAG", TRUE, FALSE)
V4R <- extract_position("GACTACHVGGGTATCTAATCC", FALSE, TRUE)
V4F <- extract_position("GTGYCAGCMGCCGCGGTAA", TRUE, FALSE)
V5R <- extract_position("CCGYCAATTYMTTTRAGTTT", FALSE, TRUE)
V6F <- extract_position("AATTGACGGGGRCCCGC", TRUE, FALSE)
V8R <- extract_position("ACGGGCRGTGWGTRCAA", FALSE, TRUE)
V9R <- extract_position("TACGGYTACCTTGTTAYGACTT", FALSE, TRUE)
START <- 1
END <- length(dna_seqs[[1]])

##########################################################################
# Make separate fasta files of all different V-regions for each genus
##########################################################################

extract_and_write <- function(start, end, prefix) {
  # Extract the subsequence
  subseqs <- lapply(seq_along(dna_seqs), function(i) {
    subseq(dna_seqs[[i]], start = start, end = end)
  })
  sequences <- DNAStringSet(subseqs, use.names = TRUE)
  names(sequences) <- names(dna_seqs)
  
  # Split the sequences by genus
  sequences_by_genus <- split(sequences, sapply(strsplit(names(sequences), " "), `[`, 2))
  
  # Write each group of sequences to a separate file
  lapply(names(sequences_by_genus), function(genus) {
    writeXStringSet(sequences_by_genus[[genus]], format = "fasta", filepath = paste0("Data/",genus,"/",genus, "_", prefix, ".fasta"))
  })
}

extract_and_write(START, END, "full16S")
extract_and_write(V1F, V3R, "V1V3")
extract_and_write(V3F, V4R, "V3V4")
extract_and_write(V4F, V4R, "V4")
extract_and_write(V4F, V5R, "V4V5")
extract_and_write(V6F, V8R, "V6V8")
extract_and_write(V6F, V9R, "V6V9")