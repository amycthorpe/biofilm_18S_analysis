# biofilm analysis for manuscripts - 18S
# titan on R bash neohuxley


# load libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(TITAN2)
library(lubridate)


set.seed(123)

dry_run <- FALSE # if TRUE then check samples and filtering only and skip TITAN


path <- "/prj/freshwaters/biofilms/metabarcoding/18S_2021_23"


# import data
ASV_tab <- data.frame(fread(file.path(path, "data/ASV_counts_18S_processed.csv"), sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.table(file.path(path, "data/ASV_tax_18S_processed.csv"), header=T, row.names=1, check.names=F, sep=",")
meta_all <- read.table(file.path(path, "data/seq_metadata_all_amplicons.csv"), header=T, check.names=F, sep=",")

meta <- meta_all[meta_all$`18S_ID` %in% rownames(ASV_tab), ] # subset metadata to seq ID (change depending on amplicon)
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="18S_ID") # make seq ID rownames

meta <- meta %>% # add season columm
  mutate(
    season = case_when(
      month(Sample_date) %in% c(03, 04, 05, 06, 07, 08) ~ "Spring",
      month(Sample_date) %in% c(09, 10, 11, 12, 01, 02) ~ "Autumn"
    )
  )

meta$Sample_date <- as.Date(meta$Sample_date, format = "%d/%m/%Y")  # change format to date

meta <- meta[!grepl("_2024-18S", rownames(meta), fixed = T), ] # drop 2024 samples


# convert to relative abundance
ASV_rel <- ASV_tab %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(across(everything(), ~ . / sum(.))) %>%
  t(.) # samples as rows


# environmental metadata
meta_merged <- read.csv(file.path(path, "data/meta_merged.csv"), header = T, check.names = F) %>%
  rename("18S_ID" = 1) %>% 
  select(`18S_ID`, Alkalinity_combined, Nitrate_N, DO, Water_temp) %>% 
  filter(!is.na(`18S_ID`) & `18S_ID` != "") %>% 
  column_to_rownames("18S_ID")

meta_merged <- as.data.frame(na.omit(meta_merged)) # remove rows with NA


# as loop
variables <- c("Alkalinity_combined", "Nitrate_N", "DO", "Water_temp")

titan_results <- list()

for (var in variables) {
  
  cat("\n====================================\n")
  cat("Running TITAN for:", var, "\n")
  cat("Started:", as.character(Sys.time()), "\n")
  
  # subset metadata to the current variable
  meta_sub <- meta_merged[var]
  
  # subset ASV table to matching samples
  ASV_sub <- as.data.frame(ASV_rel[rownames(ASV_rel) %in% rownames(meta_sub), ]) 
  meta_sub <- as.data.frame(meta_sub[intersect(rownames(meta_sub), rownames(ASV_rel)), , drop=F])
  ASV_sub <- ASV_sub[match(rownames(meta_sub), rownames(ASV_sub)), ]
  
  # increase minimum occurrence frequency to present in >10% sites
  min_occ <- ceiling(nrow(meta_sub) * 0.1)
  ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= min_occ]

  cat("Sample names and order match:", identical(rownames(meta_sub), rownames(ASV_sub4)), "\n")
  cat("Samples:", nrow(meta_sub), "\n")
  cat("Min occurrence (10%): present in >", min_occ, "samples\n")
  cat("ASVs after filtering:", ncol(ASV_sub4), "\n")
  
  if (dry_run) next # skip TITAN and saving
  
  # run TITAN
  titan_out <- titan(meta_sub, ASV_sub4, minSplt=5, numPerm=500, boot=T, nBoot=500, imax=F, ivTot=F, pur.cut=0.95, rel.cut=0.95, ncpus=52)
  
  # save in list
  titan_results[[var]] <- titan_out
  
  # save individual RDS
  saveRDS(titan_out, file = file.path(path, "titan", paste0("titan_", var, "_10perc.rds")))
  
  cat("Finished:", as.character(Sys.time()), "\n")
}


# run by season
seasons <- c("Spring", "Autumn")

titan_results_season <- list()

for (s in seasons) {
  
  # sample IDs for the current season
  season_ids <- rownames(meta)[meta$season == s]
  
  for (var in variables) {
    
    cat("\n====================================\n")
    cat("Running TITAN for:", var, "-", s, "\n")
    cat("Started:", as.character(Sys.time()), "\n")
    
    # subset metadata to the current variable and season
    meta_sub <- meta_merged[rownames(meta_merged) %in% season_ids, var, drop=F]
    
    # subset ASV table to matching samples
    ASV_sub <- as.data.frame(ASV_rel[rownames(ASV_rel) %in% rownames(meta_sub), ]) 
    meta_sub <- as.data.frame(meta_sub[intersect(rownames(meta_sub), rownames(ASV_rel)), , drop=F])
    ASV_sub <- ASV_sub[match(rownames(meta_sub), rownames(ASV_sub)), ]
    
    # minimum occurrence frequency, present in >10% of samples in this season
    min_occ <- ceiling(nrow(meta_sub) * 0.1)
    ASV_sub4 <- ASV_sub[, colSums(ASV_sub > 0) >= min_occ]
    
    cat("Sample names and order match:", identical(rownames(meta_sub), rownames(ASV_sub4)), "\n")
    cat("Samples:", nrow(meta_sub), "\n")
    cat("Min occurrence (10%): present in >", min_occ, "samples\n")
    cat("ASVs after filtering:", ncol(ASV_sub4), "\n")
    
    if (dry_run) next # skip TITAN and saving
    
    # run TITAN
    titan_out <- titan(meta_sub, ASV_sub4, minSplt=5, numPerm=500, boot=T, nBoot=500, imax=F, ivTot=F, pur.cut=0.95, rel.cut=0.95, ncpus=52)
    
    # save in list
    titan_results_season[[paste(var, tolower(s), sep = "_")]] <- titan_out
    
    # save individual RDS
    saveRDS(titan_out, file = file.path(path, "titan", paste0("titan_", var, "_10perc_", tolower(s), ".rds")))
    
    cat("Finished:", as.character(Sys.time()), "\n")
  }
}

