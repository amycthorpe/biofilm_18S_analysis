# DADA2 for big data https://benjjneb.github.io/dada2/bigdata_paired.html
# Amy Thorpe
# Samples: biofilms 2021-22, 2023, optimisation samples
# Amplicon: 18S, ITS, rbcL (16S done separately for manuscript in /EA_16S_all_DADA2_repeated_for_manuscript, and 2024 samples done separately in /all_panel - both maintain consistent parameters)
# Platform: NextSeq 2000 (JHI)
# Flow cell: P1
# Cutadapt: yes


# set variables - only need to change this part of the script
run <- "run_2021_22" # run_2021_22, run_2023, optimisation
amplicon <- "18S" # 18S, ITS, rbcL

truncLens <- list("16S" = c(230, 220),
                  "18S" = c(230, 220),
                  "ITS" = c(230, 220),
                  "rbcL" = c(230, 220))

truncLen <- truncLens[[amplicon]]

input <- paste0("/hdd0/amytho/EA_biofilm_nextseqruns/", run, "/", amplicon, "/fastq") # path to data
output <- paste0("/hdd0/amytho/EA_biofilm_nextseqruns/", run, "/", amplicon, "/DADA2_output") # path to output


# set variables for taxonomy assignment if not merging tables first
minfold <- 4 # set minfoldoverabundance depending on % chimera
minboot <- 80 # set bootstrap value for taxonomy assignment
db_paths <- list("16S" = "/hdd0/amytho/training_datasets/silva_nr99_v138.2_toGenus_trainset.fa.gz", # downloaded from https://zenodo.org/records/14169026
                 "18S" = "/hdd0/amytho/training_datasets/pr2_version_5.1.1_SSU_dada2.fasta.gz", # downloaded from https://zenodo.org/records/14169026
                 "ITS" = "/hdd0/amytho/training_datasets/sh_general_release_dynamic_all_19.02.2025.fasta") # downloaded from https://unite.ut.ee/repository.php
# rbcL handled differently


# load libraries and set seed
library(dada2)
set.seed(123)


# file parsing
pathF <- file.path(input, "cutadapt/trimmed/forward") # check files
pathR <- file.path(input, "cutadapt/trimmed/reverse")
list.files(pathF)
list.files(pathR)

filtpathF <- file.path(pathF, "filtered") # set up subdirectory for filtered reads
filtpathR <- file.path(pathR, "filtered")

fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq")) # sort ensures forward/reverse reads are in same order
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq"))

sample.names <- sub("\\.fastq\\.gz$", "", basename(fastqFs)) # get sample names

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.") # check files


# filtering: optimise parameters by sequencing run
fastqFs_full <- sort(list.files(pathF, pattern="_R1_001.fastq", full.names=T)) # need full.names=T for plots to work
fastqRs_full <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names=T))

# plotQualityProfile(fastqFs_full[1:12]) # check where Q score drops below Q30 for truncating
# plotQualityProfile(fastqRs_full[1:12])

out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                     rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                     truncLen=truncLen,
                     maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE, # high stringency
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

write.csv(out, file.path(output, "filtered_out.csv"))


# file parsing
filtFs <- list.files(filtpathF, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_R2_001.fastq", full.names = TRUE)

sample.namesF <- sub("\\.fastq\\.gz$", "", basename(filtFs)) # get filtered sample names
sample.namesR <- sub("\\.fastq\\.gz$", "", basename(filtRs))

names(filtFs) <- sample.names
names(filtRs) <- sample.names


# learn the error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE) # a million 100nt reads (or 100M total bases) is more than adequate to learn the error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)


# sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose = TRUE)
  ddF <- dada(derepF, err=errF, multithread=TRUE, verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose = TRUE)
  ddR <- dada(derepR, err=errR, multithread=TRUE, verbose = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)


# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # check table dimensions
table(nchar(getSequences(seqtab))) # inspect distribution of sequence lengths 

write.csv(seqtab, file = file.path(output, "seqtab.csv"))


# proceed to merge seqtabs across runs by amplicon, chimera removal, then taxonomy assignment


# set variables - only need to change this part of the script
run <- "run_2023" # run_2021_22, run_2023, optimisation
input <- paste0("/hdd0/amytho/EA_biofilm_nextseqruns/", run, "/", amplicon, "/fastq") # path to data
output <- paste0("/hdd0/amytho/EA_biofilm_nextseqruns/", run, "/", amplicon, "/DADA2_output") # path to output


# set variables for taxonomy assignment if not merging tables first
minfold <- 4 # set minfoldoverabundance depending on % chimera
minboot <- 80 # set bootstrap value for taxonomy assignment
db_paths <- list("16S" = "/hdd0/amytho/training_datasets/silva_nr99_v138.2_toGenus_trainset.fa.gz", # downloaded from https://zenodo.org/records/14169026
                 "18S" = "/hdd0/amytho/training_datasets/pr2_version_5.1.1_SSU_dada2.fasta.gz", # downloaded from https://zenodo.org/records/14169026
                 "ITS" = "/hdd0/amytho/training_datasets/sh_general_release_dynamic_all_19.02.2025.fasta") # downloaded from https://unite.ut.ee/repository.php
# rbcL handled differently


# load libraries and set seed
library(dada2)
set.seed(123)


# file parsing
pathF <- file.path(input, "cutadapt/trimmed/forward") # check files
pathR <- file.path(input, "cutadapt/trimmed/reverse")
list.files(pathF)
list.files(pathR)

filtpathF <- file.path(pathF, "filtered") # set up subdirectory for filtered reads
filtpathR <- file.path(pathR, "filtered")

fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq")) # sort ensures forward/reverse reads are in same order
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq"))

sample.names <- sub("\\.fastq\\.gz$", "", basename(fastqFs)) # get sample names

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.") # check files


# filtering: optimise parameters by sequencing run
fastqFs_full <- sort(list.files(pathF, pattern="_R1_001.fastq", full.names=T)) # need full.names=T for plots to work
fastqRs_full <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names=T))

# plotQualityProfile(fastqFs_full[1:12]) # check where Q score drops below Q30 for truncating
# plotQualityProfile(fastqRs_full[1:12])

out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                     rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                     truncLen=truncLen,
                     maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE, # high stringency
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

write.csv(out, file.path(output, "filtered_out.csv"))


# file parsing
filtFs <- list.files(filtpathF, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_R2_001.fastq", full.names = TRUE)

sample.namesF <- sub("\\.fastq\\.gz$", "", basename(filtFs)) # get filtered sample names
sample.namesR <- sub("\\.fastq\\.gz$", "", basename(filtRs))

names(filtFs) <- sample.names
names(filtRs) <- sample.names


# learn the error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE) # a million 100nt reads (or 100M total bases) is more than adequate to learn the error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)


# sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose = TRUE)
  ddF <- dada(derepF, err=errF, multithread=TRUE, verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose = TRUE)
  ddR <- dada(derepR, err=errR, multithread=TRUE, verbose = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)


# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # check table dimensions
table(nchar(getSequences(seqtab))) # inspect distribution of sequence lengths 

write.csv(seqtab, file = file.path(output, "seqtab.csv"))


# proceed to merge seqtabs across runs by amplicon, chimera removal, then taxonomy assignment


# set variables - only need to change this part of the script
run <- "run_optimisation" # run_2021_22, run_2023, run_optimisation
input <- paste0("/hdd0/amytho/EA_biofilm_nextseqruns/", run, "/", amplicon, "/fastq") # path to data
output <- paste0("/hdd0/amytho/EA_biofilm_nextseqruns/", run, "/", amplicon, "/DADA2_output") # path to output


# set variables for taxonomy assignment if not merging tables first
minfold <- 4 # set minfoldoverabundance depending on % chimera
minboot <- 80 # set bootstrap value for taxonomy assignment
db_paths <- list("16S" = "/hdd0/amytho/training_datasets/silva_nr99_v138.2_toGenus_trainset.fa.gz", # downloaded from https://zenodo.org/records/14169026
                 "18S" = "/hdd0/amytho/training_datasets/pr2_version_5.1.1_SSU_dada2.fasta.gz", # downloaded from https://zenodo.org/records/14169026
                 "ITS" = "/hdd0/amytho/training_datasets/sh_general_release_dynamic_all_19.02.2025.fasta") # downloaded from https://unite.ut.ee/repository.php
# rbcL handled differently


# load libraries and set seed
library(dada2)
set.seed(123)


# file parsing
pathF <- file.path(input, "cutadapt/trimmed/forward") # check files
pathR <- file.path(input, "cutadapt/trimmed/reverse")
list.files(pathF)
list.files(pathR)

filtpathF <- file.path(pathF, "filtered") # set up subdirectory for filtered reads
filtpathR <- file.path(pathR, "filtered")

fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq")) # sort ensures forward/reverse reads are in same order
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq"))

sample.names <- sub("\\.fastq\\.gz$", "", basename(fastqFs)) # get sample names

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.") # check files


# filtering: optimise parameters by sequencing run
fastqFs_full <- sort(list.files(pathF, pattern="_R1_001.fastq", full.names=T)) # need full.names=T for plots to work
fastqRs_full <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names=T))

# plotQualityProfile(fastqFs_full[1:12]) # check where Q score drops below Q30 for truncating
# plotQualityProfile(fastqRs_full[1:12])

out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                     rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                     truncLen=truncLen,
                     maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE, # high stringency
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

write.csv(out, file.path(output, "filtered_out.csv"))


# file parsing
filtFs <- list.files(filtpathF, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_R2_001.fastq", full.names = TRUE)

sample.namesF <- sub("\\.fastq\\.gz$", "", basename(filtFs)) # get filtered sample names
sample.namesR <- sub("\\.fastq\\.gz$", "", basename(filtRs))

names(filtFs) <- sample.names
names(filtRs) <- sample.names


# learn the error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE) # a million 100nt reads (or 100M total bases) is more than adequate to learn the error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)


# sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose = TRUE)
  ddF <- dada(derepF, err=errF, multithread=TRUE, verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose = TRUE)
  ddR <- dada(derepR, err=errR, multithread=TRUE, verbose = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)


# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # check table dimensions
table(nchar(getSequences(seqtab))) # inspect distribution of sequence lengths 

write.csv(seqtab, file = file.path(output, "seqtab.csv"))


# proceed to merge seqtabs across runs by amplicon, chimera removal, then taxonomy assignment
