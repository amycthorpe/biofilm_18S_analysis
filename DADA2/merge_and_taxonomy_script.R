# Merge and assign taxonomy
# Amy Thorpe
# Samples: biofilms
# Amplicon: 18S


# set variables - only need to change this part of the script
amplicon <- "18S" # 16S, 18S, ITS, rbcL

input <- "/hdd0/amytho/EA_biofilm_nextseqruns" # path to data
run1_path <- file.path(input, "all_amplicons_repeated_for_consistency/run_2021_22", amplicon, "DADA2_output/")
run2_path <- file.path(input, "all_amplicons_repeated_for_consistency/run_2023", amplicon, "DADA2_output/")
run3_path <- file.path(input, "all_amplicons_repeated_for_consistency/run_optimisation", amplicon, "DADA2_output/")

run_paths <- list(run1_path, run2_path, run3_path)

output <- paste0("/hdd0/amytho/EA_biofilm_nextseqruns/all_amplicons_repeated_for_consistency/merged/", amplicon) # path to output

minfold <- 12 # set minfoldoverabundance depending on % chimera
minboot <- 80 # set bootstrap value for taxonomy assignment
db_paths <- list("16S" = "/hdd0/amytho/training_datasets/silva_nr99_v138.2_toGenus_trainset.fa.gz", # downloaded from https://zenodo.org/records/14169026
                 "18S" = "/hdd0/amytho/training_datasets/pr2_version_5.1.1_SSU_dada2.fasta.gz", # downloaded from https://zenodo.org/records/14169026
                 "ITS" = "/hdd0/amytho/training_datasets/sh_general_release_dynamic_all_19.02.2025.fasta") # downloaded from https://unite.ut.ee/repository.php
# rbcL handled differently


# load libraries and set seed
library(dada2)
library(data.table)
library(dplyr)
library(tibble)

set.seed(123)


# load seq tabs (samples need to be rows)
seqtabs <- lapply(run_paths, function(p) {
  fread(file.path(p, "seqtab.csv")) %>%
    column_to_rownames(var="V1") %>%
    as.matrix()
})


# for 16S
# seqtabs <- lapply(run_paths, function(p) {
#   f <- list.files(p, pattern = "^seqtab.*\\.csv$", full.names = TRUE)[1]
#   fread(f) %>%
#     column_to_rownames(var = "V1") %>%
#     as.matrix()
# })

seqtabs[1:3] <- lapply(seqtabs[1:3], t) # transpose first three seqtabs so samples are rows
seqtabs[[3]] <- seqtabs[[3]][grepl("^Z", rownames(seqtabs[[3]])), ] # drop non zymo samples from optimisation
seqtabs[[3]] <- seqtabs[[3]][ , colSums(seqtabs[[3]]) != 0] # drop ASVs only found in non zymo samples (now zero abundance)


# merge runs
seqtab.all <- do.call(mergeSequenceTables, seqtabs)


# remove chimeras after merge
seqtab.all.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance=minfold)

dim(seqtab.all.nochim) # check table dimensions
sum(seqtab.all.nochim)/sum(seqtab.all) # check proportion non-chimera


# transpose and save merged seqtab
seqtab.all.nochim_t <- t(seqtab.all.nochim)
write.csv(seqtab.all.nochim_t, file.path(output, paste0("seqtab_merged_", amplicon, ".csv")))


# assign taxonomy
if (amplicon == "rbcL") {
  # diatbarcode
  if (!requireNamespace("diatbarcode", quietly = TRUE)) {
    remotes::install_github("fkeck/diatbarcode")
  }

  library(diatbarcode)

  tax_fas <- diatbarcode::download_diatbarcode(flavor = "rbcl312_dada2_tax")

  taxa_full <- assignTaxonomy(seqtab.all.nochim, tax_fas$path, minBoot = minboot,
                              taxLevels = c("Empire","Kingdom","Subkingdom","Phylum","Class","Order","Family","Genus","Species"),
                              outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)

  taxa <- taxa_full$tax
  boot <- taxa_full$boot
} else {
  db_path <- db_paths[[amplicon]]
  taxa <- assignTaxonomy(seqtab.all.nochim, db_path, minBoot = minboot, multithread = TRUE)
}

asv_seqs <- colnames(seqtab.all.nochim)
asv_headers <- vector(dim(seqtab.all.nochim)[2], mode="character")

for (i in 1:dim(seqtab.all.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


# get fasta
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(output, paste0("ASVs_merged_", amplicon, ".fa")))


# get counts tab
asv_tab <- t(seqtab.all.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(output, paste0("counts_merged_", amplicon, ".tsv")),
            sep="\t", quote=F, col.names=NA)


# needed to run 16S taxonomy assignment in chunks due to dataset size
# chunk_size <- 50000
# n_asv <- ncol(seqtab.all.nochim)
# idx <- split(seq_len(n_asv), ceiling(seq_len(n_asv) / chunk_size)) # split column indices into chunks
# 
# tax_list <- vector("list", length(idx))
# 
# for (i in seq_along(idx)) { # loop over each chunk
#   cat("Processing chunk", i, "of", length(idx), "\n")
#   
#   seqtab_chunk <- seqtab.all.nochim[, idx[[i]], drop = FALSE]
#   
#   tax_list[[i]] <- assignTaxonomy(
#     seqtab_chunk,
#     db_path,
#     minBoot = minboot,
#     multithread = TRUE
#   )
# }

# taxa <- do.call(rbind, tax_list) # combine taxonomy tables
# 
# taxa <- taxa[colnames(seqtab.all.nochim), ] # reorder to match original ASV order
# 
# taxa.print <- taxa # removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print) # check taxonomy

# get tax tab
rownames(taxa) <- gsub(pattern=">", replacement="", x=asv_headers) 
write.table(taxa, file.path(output, paste0("taxonomy_merged_", amplicon, ".tsv")),
            sep="\t", quote=F, col.names=NA)

if (amplicon == "18S") {
  # set standard 18S column names
  colnames(taxa) <- c("Domain", "Supergroup", "Division", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  write.table(taxa, file.path(output, paste0("taxonomy_merged_", amplicon, ".tsv")),
              sep = "\t", quote = F, col.names = NA)
} else if (amplicon %in% c("16S", "ITS", "rbcL")) {
  write.table(taxa, file.path(output, paste0("taxonomy_merged_", amplicon, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
}



