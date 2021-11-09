library(tibble)
library(dplyr)
library(DECIPHER)
library(Biostrings)

path <- "/mnt/scratch/f0008425/SSU/filtered"

set.seed(1390)
CORES = parallel::detectCores()

seqtab <- readRDS(file='seq_table.RDS')
seqtab <- t(seqtab)
asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = CORES)
d <- DECIPHER::DistanceMatrix(aln, processors = CORES)
clusters <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU 
  processors = CORES)

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>%
  add_column(sequence = asv_sequences)
merged_seqtab <- seqtab %>% 
  t %>%
  rowsum(clusters$cluster) %>%
  t
# Optional renaming of clusters to OTU<cluster #>
colnames(merged_seqtab) <- paste0("OTU", colnames(merged_seqtab))

write.csv(t(merged_seqtab), file="OTU_table.csv")
write.csv(clusters, file="cluster_results.csv")


