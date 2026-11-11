# biofilm analysis - 18S
# taxonomy clean up, ASV-genus labels for plots, and EukFunc assignment with corrections


# load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(data.table)
library(EukFunc)

set.seed(123)


# import data
tax_tab <- read.table("data/processed_2021_23/ASV_tax_18S_processed.csv", header=T, row.names=1, check.names=F, sep=",")
ranks <- colnames(tax_tab)


# remove _X's
x_pattern <- "_X+(?=_|$)" # _X, _XX, _XXX etc at the end of a name or before another _

x_issues <- tax_tab %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "rank", values_to = "before") %>%
  filter(grepl(x_pattern, before, perl = TRUE)) %>%
  mutate(rank = factor(rank, levels = ranks))

n_distinct(x_issues$ASV) # ASVs with _X's
table(x_issues$rank)

tax_x <- tax_tab
tax_x[] <- lapply(tax_x, function(x) gsub(x_pattern, "", x, perl = TRUE))

x_issues$after <- tax_x[cbind(x_issues$ASV, as.character(x_issues$rank))]
sum(grepl(x_pattern, unlist(tax_x), perl = TRUE)) # check, should be 0


# names repeated at lower ranks become unassigned (from Class down, Phylum may repeat Division)
repeat_rank <- matrix(FALSE, nrow = nrow(tax_x), ncol = length(ranks), dimnames = list(rownames(tax_x), ranks))

for (r in 5:9) { # Class to Species
  repeat_rank[, r] <- tax_x[[r]] != "unassigned" & rowSums(tax_x[, 1:(r - 1)] == tax_x[[r]]) > 0
}

repeat_issues <- as.data.frame(which(repeat_rank, arr.ind = TRUE)) %>%
  mutate(ASV = rownames(tax_x)[row], rank = factor(ranks[col], levels = ranks), before = tax_x[cbind(row, col)]) %>%
  select(ASV, rank, before) %>%
  arrange(rank)

n_distinct(repeat_issues$ASV) # ASVs with a repeated name
table(repeat_issues$rank)

tax_fix <- tax_x
tax_fix[repeat_rank] <- "unassigned"

repeat_issues$after <- tax_fix[cbind(repeat_issues$ASV, as.character(repeat_issues$rank))]

repeat_after <- matrix(FALSE, nrow = nrow(tax_fix), ncol = length(ranks))
for (r in 5:9) {
  repeat_after[, r] <- tax_fix[[r]] != "unassigned" & rowSums(tax_fix[, 1:(r - 1)] == tax_fix[[r]]) > 0
}
sum(repeat_after) # check, should be 0


# genus names that are still not genera (PR2 clades) and higher ranks left in other columns
clade_pattern <- "lineage|[Cc]lade|Group|Novel|[0-9_/-]" # genus names that are really PR2 clade names

genus_issues <- tax_fix %>%
  rownames_to_column("ASV") %>%
  filter(Genus != "unassigned", grepl(clade_pattern, Genus)) %>%
  select(ASV, Class, Order, Family, Genus)

n_distinct(genus_issues$ASV) # these keep their name but are not italicised


# labels for plotting
# unassigned genera are named by the last assigned rank above genus
last_rank <- apply(tax_fix[, 1:7] != "unassigned", 1, function(x) max(which(x)))
last_name <- tax_fix[cbind(seq_len(nrow(tax_fix)), last_rank)]

labels <- tax_fix %>%
  mutate(ASV = rownames(tax_fix),
         is_genus = Genus != "unassigned" & !grepl(clade_pattern, Genus), # a real genus, so italicised
         assigned_to = if_else(Genus != "unassigned", "Genus", ranks[last_rank]),
         genus_label = if_else(Genus != "unassigned", Genus, paste0("unassigned_", last_name)),
         ASV_label = paste0(ASV, ": ", genus_label),
         genus_label_italic = if_else(is_genus, glue::glue("italic('{genus_label}')"), glue::glue("'{genus_label}'")), 
         ASV_label_italic = if_else(is_genus, glue::glue("'{ASV}: '*italic('{genus_label}')"), glue::glue("'{ASV_label}'")))

table(labels$assigned_to)
table(labels$is_genus)
labels %>% filter(!is_genus) %>% count(genus_label, sort = TRUE) %>% head(20)


# EukFunc
# assign functions on the original PR2 names as EukFunc matches on these
tax_func <- tax_tab %>%
  rename(domain = Domain,
         supergroup = Supergroup,
         division = Division,
         subdivision = Phylum,
         class = Class,
         order = Order,
         family = Family,
         genus = Genus,
         species = Species)

func_tab <- assign_path(tax_func)
rownames(func_tab) <- rownames(tax_tab)

table(func_tab$main_functional_class, useNA = "always")

# correction rules
fill_phototroph_phyla <- c("Chlorophyta") # unassigned chlorophytes -> phototroph
fill_phototroph_classes <- c("Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae", "Phaeophyceae") # unassigned diatoms and brown algae -> phototroph
fill_heterotroph_phyla <- c("Fungi") #, "Choanoflagellata", "Discosea", "Apusomonada", "Tubulinea", "Apicomplexa", "Ichthyosporea") # unassigned ASVs stay unassigned at detailed level but heterotroph at broad level

trophic_levels <- c("Phototroph", "Mixotroph", "Consumer", "Parasite/symbiont", "Saprotroph", "Unassigned")
broad_levels <- c("Phototroph", "Heterotroph", "Mixotroph", "Unassigned")


# trophic_group: detailed group
# trophic_broad: phototroph / heterotroph / mixotroph / unassigned
trophic <- labels %>%
  mutate(main_functional_class = func_tab[ASV, "main_functional_class"],
         detailed_functional_class = func_tab[ASV, "detailed_functional_class"],
         eukfunc_class = if_else(is.na(main_functional_class) | main_functional_class == "", "unassigned", main_functional_class),
         trophic_group = case_when(
           grepl("mixotroph", detailed_functional_class) ~ "Mixotroph",
           eukfunc_class == "phototroph" ~ "Phototroph",
           eukfunc_class == "predator" ~ "Consumer",
           eukfunc_class == "symbiotroph" ~ "Parasite/symbiont",
           eukfunc_class == "saprotroph" ~ "Saprotroph",
           Class %in% fill_phototroph_classes ~ "Phototroph", # filled
           Phylum %in% fill_phototroph_phyla ~ "Phototroph", # filled
           TRUE ~ "Unassigned"),
         trophic_broad = case_when(
           trophic_group == "Phototroph" ~ "Phototroph",
           trophic_group %in% c("Consumer", "Parasite/symbiont", "Saprotroph") ~ "Heterotroph",
           trophic_group == "Mixotroph" ~ "Mixotroph",
           Phylum %in% fill_heterotroph_phyla ~ "Heterotroph", # filled (broad only)
           TRUE ~ "Unassigned"),
         assignment = case_when(
           grepl("mixotroph", detailed_functional_class) | eukfunc_class %in% c("phototroph", "predator", "symbiotroph", "saprotroph") ~ "EukFunc",
           trophic_group != "Unassigned" ~ "filled from taxonomy",
           trophic_broad == "Heterotroph" ~ "filled from taxonomy (broad only)",
           TRUE ~ "none"),
         trophic_group = factor(trophic_group, levels = trophic_levels),
         trophic_broad = factor(trophic_broad, levels = broad_levels))

table(trophic$trophic_group, trophic$assignment) # check
table(trophic$trophic_group, trophic$trophic_broad) # check


# export
write.csv(func_tab, "data/processed_2021_23/eukfunc.csv")
write.csv(tax_fix, "data/processed_2021_23/taxonomy_clean.csv") # ranks only, same layout as ASV_tax_18S_processed.csv
write.csv(trophic %>% select(-all_of(ranks)), "data/processed_2021_23/taxonomy_labels.csv", row.names = FALSE) # labels and trophic groups (join to taxonomy_clean.csv by ASV)

