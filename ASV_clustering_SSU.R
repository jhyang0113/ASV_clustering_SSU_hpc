library(tibble)
library(dplyr)
library(DECIPHER)
library(Biostrings)

path <- "/mnt/scratch/f0008425/SSU/filtered"

set.seed(1390)
CORES = parallel::detectCores()

# Load dataset
seqtab <- readRDS(file='seq_table.RDS')
taxtab <- readRDS(file='tax_table.RDS')
seqtab <- t(seqtab)
asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

# ASV clustering (use `cutoff = 0.03` for a 97% OTU)
alnseq <- DECIPHER::AlignSeqs(dna, processors = nprocessor)
dm <- DECIPHER::DistanceMatrix(alnseq, processors = nprocessor)
clusters <- DECIPHER::IdClusters(dm, method = "complete", cutoff = 0.03, processors = nprocessor)

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

# Representative sequences extracting (new OTU, representative ASV, its sequences)
newdf<-as.data.frame(cbind.data.frame(clusters$cluster, rownames(clusters),clusters$sequence))
colnames(newdf)<-c("OTU","ASV","sequence") 
newdf.sort<-newdf %>% arrange(OTU) %>% distinct(OTU,.keep_all = TRUE)

# Create tax table
taxlist<-as.data.frame(paste("ASV",newdf.sort$ASV, sep=""))
colnames(taxlist)<-c("ASV")
entiretax <- as.data.frame(paste("ASV",seq(1:nrow(taxtab)), sep=""))
colnames(entiretax)<-c("ASV")
entiretax <- cbind(entiretax, as.data.frame(taxtab))
tax_table_otu <- merge(taxlist, entiretax[,c("ASV", "Kingdom", "Phylum","Class","Order", "Family", "Genus", "Species")], by="ASV")
                      
# Create the final dataset for further analysis
write.csv(t(merged_seqtab), file="otu_table_otu.csv")
write.csv(tax_table_otu, file="tax_table_otu.csv", row.names = FALSE)
write.csv(clusters, file="cluster_results.csv")
write.csv(newdf.sort,file="representative_seqs.csv", row.names=FALSE, col.names=FALSE,quote = FALSE)
