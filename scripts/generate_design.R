library(Biostrings)

design_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]

design <- read.table(design_file, header = TRUE, sep = "\t")
design$names <- gsub(" ", "_", design$names) # Some names have spaces

# Convert to DNAStringSet
a <- DNAStringSet(design$seq)
names(a) <- design$names

# The snakemake pipeline uses the reverse complement of designed seqs
a_rev <- reverseComplement(a)
a_sub <- subseq(a_rev, start=26, end=175) # Remove primers

writeXStringSet(a_sub, file = out_file, width = 150)
