library(dplyr)


# Make vector of file names
input_files <- snakemake@input
samples <- gsub(".reads.txt", "", input_files)
DNA_idx <- grep("DNA", input_files)
RNA_idx <- grep("RNA", input_files)


# Combine reads from all DNA and RNA files
combine_reads <- function(input_files) {
  reads <- read.table(input_files[[1]], header = F)
  colnames(reads) <- c(samples[1], "variant", "barcode")
  
  for (i in 2:length(input_files)) {
    r <- read.table(input_files[[i]], header = F)
    colnames(r) <- c(samples[i], "variant", "barcode")
    reads <- merge(reads, r, by = c("variant", "barcode"), all = TRUE)
  }
  reads[is.na(reads)] <- 0 # Set NA to 0
  return(reads)
}


# Filter reads
filter_reads <- function(reads) {
  reads <- reads %>%
    mutate(DNA_sum = rowSums(reads[DNA_idx + 2]),
           RNA_sum = rowSums(reads[RNA_idx + 2])) %>%
    filter(DNA_sum > length(DNA_idx), RNA_sum > 0) %>%
    group_by(variant) %>%
    mutate(n_bc = n()) %>%
    filter(n_bc >= length(DNA_idx)) %>%
    mutate(ratio = log2(RNA_sum / DNA_sum),
           ratio_med = median(ratio),
           ratio_diff = ratio - ratio_med,
           mad = mad(ratio))
    
  n_bin <- 20
  qs <- quantile(log10(reads$RNA_sum), 0:n_bin / n_bin)
  reads <- reads %>%
    mutate(bin = cut(log10(RNA_sum), qs, include.lowest = TRUE, labels = 1:n_bin)) %>%
    group_by(bin) %>%
    filter(ratio_diff < 5 * mad(ratio_diff)) %>%
    as.data.frame()
  return(reads)
}


# Summarize counts across barcodes
summarize_reads <- function(reads_filt) {
  reads_summary <- reads_filt %>%
    group_by(variant) %>%
    mutate(n_bc = 1) %>%
    summarize_at(vars(samples, n_bc), sum)
  return(as.data.frame(reads_summary))
}


# Parse output filename
output_name <- snakemake@output[[1]]
output_prefix <- gsub(".counts.by.variant.txt", "", output_name)


# Process read without filtering
reads_all <- combine_reads(input_files)
write.table(reads_all, file = paste0(output_prefix, ".counts.by.barcode.txt"), col.names=T, row.names=F, quote=F)
reads_all_summary <- summarize_reads(reads_all)
write.table(reads_all_summary, file = output_name, col.names=T, row.names=F, quote=F)


# Process reads with filtering
reads_all_filt <- filter_reads(reads_all)
write.table(reads_all_filt, file = paste0(output_prefix, ".counts.by.barcode.filtered.txt"), col.names=T, row.names=F, quote=F)
reads_all_summary <- summarize_reads(reads_all_filt)
write.table(reads_all_summary, file = paste0(output_prefix, ".counts.by.variant.filtered.txt"), col.names=T, row.names=F, quote=F)



