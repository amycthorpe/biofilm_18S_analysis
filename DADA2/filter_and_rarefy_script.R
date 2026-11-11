# Filter and rarefy
# Amy Thorpe
# Samples: biofilms


# set amplicon
amplicon <- "18S" # 16S, 18S, ITS, rbcL


# set rarefaction depth (check reads per sample etc first)
depths <- c("16S" = 10000,
            "18S" = 10000)

depth <- depths[amplicon]


# load libraries, set wd
library(plyr)
library(dplyr)
library(vegan)
library(seqinr)

set.seed(123)


# import data
tax_tab <- read.table(paste0("data/DADA2_output/taxonomy_merged_", amplicon, ".tsv"),
                      header=T, row.names=1, check.names=F, sep="\t")
counts_tab <- read.table(paste0("data/DADA2_output/counts_merged_", amplicon, ".tsv"),
                      header=T, row.names=1, check.names=F, sep="\t")
meta <- read.table("data/DADA2_output/seq_metadata_all_amplicons.csv", header=T, sep=",", fill=T, check.names=F)

id_col <- paste0(amplicon, "_ID")
colnames(counts_tab) <- sub("_S[0-9]+(_L[0-9]+)?_R1_001$", "", colnames(counts_tab)) # drop seq id suffix
meta <- meta[!is.na(meta[[id_col]]) & meta[[id_col]] != "", ] # remove non target amplicon samples
setdiff(colnames(counts_tab), meta[[id_col]])


# remove prefix (k__ etc) from taxonomy names if needed
tax_tab[] <- lapply(tax_tab, function(x) {
  gsub("^[a-z]__+", "", x)
})


# ASVs per sample
ASVs_per_sample <- as.data.frame(colSums(counts_tab))
ASVs_per_sample_meta <- merge(ASVs_per_sample, meta, by.x="row.names", by.y=id_col)


# filtering
# filter based on taxonomy
tax_tab[is.na(tax_tab)] <- "unassigned" # NA must be converted to unassigned first

if (amplicon == "16S") {
  tax_tab_filt <- subset(tax_tab, Kingdom == "Bacteria")
  tax_tab_filt <- subset(tax_tab_filt, Order != "Chloroplast")

} else if (amplicon == "18S") {
  tax_tab_filt <- subset(tax_tab, Domain == "Eukaryota")
  tax_tab_filt <- subset(tax_tab_filt, Division != "Rhodophyta_X")
  tax_tab_filt <- subset(tax_tab_filt, Phylum != "Metazoa")
  tax_tab_filt <- subset(tax_tab_filt, Class != "Embryophyceae")
  
} else if (amplicon == "ITS") {
  tax_tab_filt <- subset(tax_tab, Kingdom == "Fungi")
  
} else if (amplicon == "rbcL") {
  tax_tab_filt <- subset(tax_tab, Phylum=="Bacillariophyta")
}

# subset counts based on filtered tax (so both tax and counts tables have the same samples in after filtering)
counts_tab_filt <- counts_tab[rownames(counts_tab)%in%rownames(tax_tab_filt),]

# filter based on abundance
counts_tab_filt2 <- counts_tab_filt[rowSums(counts_tab_filt)>5,] # remove ASVs with <X reads
counts_tab_filt3 <- counts_tab_filt2[rowSums(counts_tab_filt2>0)>=5,] # remove ASVs that only appear in <X samples


# rarefaction
# check row sums (samples need to be rows so transpose first)
counts_tab_filt3 <- t(counts_tab_filt3)
rowsums <- as.data.frame(rowSums(counts_tab_filt3))
plot(sort(rowsums[,1]))

min(rowsums[,1])
mean(rowsums[,1])
hist(rowSums(counts_tab_filt3))

# rarefaction curve
plot <- rarecurve(counts_tab_filt3, step = 1000, col = "blue", cex = 0.7, subsample = 1000, label = FALSE, xlab = "Number of Sequences",
                ylab = "Observed Species Richness",
                main = "Rarefaction Curve") #, xlim = c(0, 20000))

# rarefy
counts_tab_rare <- rrarefy(counts_tab_filt3, depth)
counts_tab_rare <- as.data.frame(counts_tab_rare)

# remove samples below rarefied depth (if any)
counts_tab_rare <- counts_tab_rare[rowSums(counts_tab_rare[])==depth,]

# check rarefied rowsums - min and mean should now be the same
rowsums <- as.data.frame(rowSums(counts_tab_rare))
rowsums_meta <- merge(rowsums, meta, by.x="row.names", by.y=id_col)

plot(sort(rowsums[,1]))

min(rowsums[,1]) 
mean(rowsums[,1])

# check proportion samples dropped
((ncol(counts_tab) - nrow(counts_tab_rare)) - length(n_blanks)) / ncol(counts_tab) * 100

# remove 0 abundance ASVs from counts after rarefying (if samples below rarefaction depth were removed, these are ASVs only in those samples)
counts_tab_rare <- counts_tab_rare[,colSums(counts_tab_rare[])>0]

# subset tax to new rarefied ASV list (only if you removed samples)
counts_tab_rare <- as.data.frame(t(counts_tab_rare))
tax_tab_rare <- tax_tab[rownames(tax_tab)%in%rownames(counts_tab_rare),]

# transpose so samples are rows
counts_tab_rare <- as.data.frame(t(counts_tab_rare))


# save filtered and rarefied tables
write.csv(counts_tab_rare, paste0("data/processed_2021_23/ASV_counts_", amplicon, "_processed.csv"))
write.csv(tax_tab_rare, paste0("data/processed_2021_23/ASV_tax_", amplicon, "_processed.csv"))

