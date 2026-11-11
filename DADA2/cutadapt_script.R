# cutadapt


# make and copy fastq files to separate project/amplicon folders as required (command line)


# load libraries
library(dada2)
library(ShortRead)
library(Biostrings)


# set amplicon/run
amplicon <- "18S"


# set paths
path <- paste0("/prj/freshwaters/biofilms/", amplicon, "/fastq/")
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))


# primer seqs
primers <- list(
  # bacterial 16S (Walters 2015 modified)
  "16S" = list(FWD = "GTGYCAGCMGCCGCGGTAA",  # 19 bases
               REV = "GGACTACNVGGGTWTCTAAT"),  # 20 bases
  
  # eukaryotic 18S (Mangot 2012)
  "18S" = list(FWD = "CGCGGTAATTCCAGCTCCA",  # 19 bases
               REV = "TTGGYRAATGCTTTCGC"),  # 17 bases
  
  # fungal ITS (Ihrmark 2012)
  "ITS" = list(FWD = "GTGARTCATCGAATCTTTG",  # 19 bases
               REV = "TCCTCCGCTTATTGATATGC"),  # 20 bases
  
  # diatom rbcL (Kelly 2020)
  "rbcL" = list(FWD = "ATGCGTTGGAGAGARCGTTTC",  # 21 bases
                REV = "GATCACCTTCTAATTTACCWACAACTG"),  # 26 bases
  
  # fish 12S (Kelly and Riaz 2014)
  "12S" = list(FWD = "ACTGGGATTAGATACCCC",  # 18 bases
               REV = "TAGAACAGGCTCCTCTAG"))  # 18 bases

FWD <- primers[[amplicon]]$FWD
REV <- primers[[amplicon]]$REV


# find primer orientations
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString)) 
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


# pre filter N's
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


# count primer hits
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


# cutadapt
cutadapt <- "/home/amytho/miniconda3/envs/cutadapt/bin/cutadapt"
# cutadapt <- "/opt/miniconda3/envs/cutadapt/bin/cutadapt"
system2(cutadapt, args = "--version") # run shell commands from R


# remove primers
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) # trim FWD and the reverse-complement of REV off of R1 (forward reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) # trim REV and the reverse-complement of FWD off of R2 (reverse reads)


# run cutadapt in parallel
library(parallel)

num_cores <- detectCores() - 2

mclapply(seq_along(fnFs), function(i) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}, mc.cores = num_cores)

# check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# set minimum read length to prevent errors in DADA2 caused by 0 bp reads after cutadapt
filtpathF <- file.path(path.cut, "trimmed/forward")
filtpathR <- file.path(path.cut, "trimmed/reverse")
dir.create(filtpathF, recursive=TRUE)
dir.create(filtpathR, recursive=TRUE)

fnFs.trim <- file.path(filtpathF, basename(fnFs.cut))
fnRs.trim <- file.path(filtpathR, basename(fnRs.cut))

out <- filterAndTrim(fwd=fnFs.cut, filt=fnFs.trim,
                     rev=fnRs.cut, filt.rev=fnRs.trim,
                     minLen = 20,
                     compress=TRUE, verbose=TRUE, multithread = TRUE)

