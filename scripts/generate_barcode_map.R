library(dplyr)
library(Biostrings)
library(data.table)

# Get barcodes for sequence
input_file <- snakemake@input

d <- fread(input_file[[1]], header = FALSE, data.table = FALSE)

n <- dim(d)[2]

d$count <- rowSums(d[3:n])
d <- d[, c(1, 2, n+1)]
colnames(d) <- c("variant", "barcode", "count")

count <- split(d$count, factor(d$barcode))
threshold <- 0.5

ambiguous <- NULL
n <- length(count)
for (i in 1:n) {
  if (length(count[[i]]) > 1) {
    x <- count[[i]]
    if (max(x) < threshold) {
      ambiguous <- c(ambiguous, names(count[i]))
    }
  }
}

df_filt <- d %>% 
  filter(count >= 3, !barcode %in% ambiguous, !is.na(variant))

#length(unique(df_filt$variant))

#var_summary <- df_filt %>% group_by(variant) %>%
#  summarize(n = n())
#hist(var_summary$n)

seqstr <- DNAStringSet(df_filt$barcode) 
seqstr <- reverseComplement(seqstr)
names(seqstr) <- df_filt$variant
writeXStringSet(seqstr, file=snakemake@output[[1]]) 


