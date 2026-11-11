# biofilm analysis - 18S
# community structure: composition, diversity, core community, abundance vs prevalence


# load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(tidyverse)
library(glue)
library(phyloseq)
library(vegan)
library(microbiome)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(lme4)
library(lmerTest) 


set.seed(123)
options(scipen = 999)


# import data
ASV_tab <- data.frame(fread("data/processed_2021_23/ASV_counts_18S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.csv("data/processed_2021_23/taxonomy_clean.csv", header=T, row.names=1, check.names=F)
meta_all <- read.table("data/metadata/seq_metadata_all_amplicons.csv", header=T, check.names=F, sep=",")
trophic <- read.csv("data/processed_2021_23/taxonomy_labels.csv", header=T, check.names=F)

dir.create("output/community_structure/figures", recursive = TRUE, showWarnings = FALSE)

meta <- meta_all[meta_all$`18S_ID` %in% rownames(ASV_tab), ] # subset metadata to seq ID (change depending on amplicon)
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="18S_ID") # make seq ID rownames

meta$Sample_date <- as.Date(meta$Sample_date, format = "%d/%m/%Y")  # change format to date (moved up so month() works)

meta <- meta %>% # add season column
  mutate(
    season = case_when(
      month(Sample_date) %in% c(03, 04, 05, 06, 07, 08) ~ "Spring",
      month(Sample_date) %in% c(09, 10, 11, 12, 01, 02) ~ "Autumn"
    )
  )

meta <- meta %>%
  mutate(season_year = paste(season, Year)) # get season x year col


# add labels to taxonomy
tax_tab2 <- tax_tab
tax_tab2$taxon <- rownames(tax_tab2)
tax_tab2 <- cbind(tax_tab2, trophic[match(rownames(tax_tab2), trophic$ASV), c("is_genus", "genus_label", "ASV_label_italic")])


# phyloseq object
TAX <- tax_table(as.matrix(tax_tab2))
ASV <- otu_table(ASV_tab, taxa_are_rows = F)
samples <- sample_data(meta)
ps <- phyloseq(ASV, TAX, samples)
ps

# convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# group at phylum level
ps_phylum<-tax_glom(ps_rel, "Phylum")

# melt
melt_rel<-psmelt(ps_phylum)


# top phyla
top_phyla <- melt_rel %>%
  filter(Phylum != "unassigned") %>% # unassigned goes into Others
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Phylum) %>%
  summarise(mean = mean(Abundance)) %>%
  arrange(desc(mean)) %>%
  filter(mean >= 0.0009) # %>% # filter by abundance threshold
# slice_head(n = 10)

melt_top <- melt_rel %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others")) %>% 
  filter(Abundance > 0) # remove zeros

phylum_abundance <- melt_top %>% 
  select(Sample, Abundance, Phylum)

write.csv(phylum_abundance, "output/community_structure/phylum_abundance.csv")

mean_abundance <- phylum_abundance %>%
  group_by(Phylum) %>%
  summarise(mean = mean(Abundance),
            SD = sd(Abundance))

phylum_abundance_seasonal <- melt_rel %>% 
  select(season, Sample, Abundance, Phylum) %>% 
  group_by(season, Phylum) %>% 
  summarise(mean = mean(Abundance),
            SD = sd(Abundance))

phylum_diff <- phylum_abundance_seasonal %>%
  select(season, Phylum, mean) %>%
  pivot_wider(names_from = season,
              values_from = mean,
              names_prefix = "mean_") %>%
  mutate(diff_autumn_minus_spring = mean_Autumn - mean_Spring,
         abs_diff = abs(mean_Autumn - mean_Spring),
         higher_in = if_else(mean_Autumn > mean_Spring, "Autumn", "Spring")) %>%
  arrange(desc(abs_diff))

top_phylum_counts <- melt_top %>%
  group_by(Sample, season) %>%
  slice_max(Abundance, n = 1) %>% 
  ungroup() %>%
  count(Phylum, season, name = "num_top")

# diatom_sample <- melt_rel %>%
#   mutate(diatom = if_else(Class %in% c("Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae"), Abundance, 0)) %>%
#   group_by(Sample) %>%
#   summarise(diatom = sum(diatom), .groups = "drop")
# 
# mean(diatom_sample$diatom); sd(diatom_sample$diatom)

# set colours
phylum_colour_map <- data.frame(
  Phylum = c("Bigyra",  "Cercozoa", "Chlorophyta", "Choanoflagellata", "Ciliophora", "Cryptophyta", "Fungi",   "Gyrista", "Streptophyta", "x__Others"),
  nord          = c("#B48EAD", "#D08770", "#A3BE8C", "#8FBCBB", "#EBCB8B", "#81A1C1", "#BF616A", "#5E81AC", "#4C566A", "#D9D9D9"),
  retro_coastal = c("#264653", "#F4A261", "#E9C46A", "#A8DADC", "#E76F51", "#B5838D", "#457B9D", "#2A9D8F", "#8AB17D", "#D9D9D9"),
  scandi        = c("#4F5D75", "#E8B4B8", "#9CAF88", "#D9C5A0", "#E1B858", "#9DC3C2", "#C8745E", "#7393B3", "#8E6C8A", "#D9D9D9"),
  meadow        = c("#4E7D6B", "#A889C2", "#8DB86B", "#F2A7B8", "#F2C75C", "#9ED0E6", "#E26D5C", "#6F8FCF", "#C98C6B", "#D9D9D9"),
  citrus_slate  = c("#3E4A61", "#C9A0DC", "#7FB77E", "#8CCFC4", "#5B7DB1", "#E98FA8", "#F6D36B", "#F28C6B", "#A7C4A0", "#D9D9D9"),
  stringsAsFactors = FALSE)

palette_choice <- "nord"
phylum_colour_map$colour <- phylum_colour_map[[palette_choice]]
rownames(phylum_colour_map) <- phylum_colour_map$Phylum
phylum_colours <- setNames(phylum_colour_map$colour, phylum_colour_map$Phylum)

write.csv(phylum_colour_map, "output/community_structure/phylum_colour_map.csv")


# community composition
season_labels <- melt_top %>%
  group_by(season_year) %>%
  arrange(Sample_date) %>%
  slice(1) %>%
  ungroup() %>%
  select(Sample, season_year) # get first sample of each season x year for labelling

melt_top <- melt_top[order(melt_top$Sample_date), ] # order chronologically
melt_top$Sample <- factor(melt_top$Sample, levels=unique(melt_top$Sample[order(melt_top$Sample_date)])) # reorder samples to match

composition_plot <- ggplot(melt_top, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  geom_col(width = 1, colour = NA) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "bottom") +
  labs(y = "Relative abundance") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks = season_labels$Sample,
                   labels = season_labels$season_year) +
  scale_fill_manual(values = phylum_colours,
                    guide = guide_legend(ncol = 2),
                    labels = function(x) ifelse(x == "x__Others", "Others", x)) +
  geom_vline(xintercept = which(levels(melt_top$Sample) %in% season_labels$Sample), 
             linetype = "dashed", color = "white", size = 0.2)

composition_plot


# seasonal composition (mean relative abundance across all samples in each season)
season_n <- table(meta$season) # samples per season

phylum_season <- melt_rel %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others")) %>%
  group_by(Sample, season, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>% # one value per phylum per sample (zeros kept)
  group_by(season, Phylum) %>%
  summarise(mean = mean(Abundance), .groups = "drop") %>%
  mutate(season = factor(season, levels = c("Spring", "Autumn")),
         label = if_else(mean >= 0.03, paste0(round(mean * 100), "%"), "")) # label phyla making up at least 3%

write.csv(phylum_season, "output/community_structure/phylum_season_mean.csv")

phylum_season_plot <- ggplot(phylum_season, aes(x = season, y = mean, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 1.6, colour = "grey15") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none") +
  labs(y = "Mean relative abundance") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 1) +
  scale_fill_manual(values = phylum_colours)

phylum_season_plot


# seasonal trophic composition
trophic_levels <- c("Phototroph", "Mixotroph", "Consumer", "Parasite/symbiont", "Saprotroph", "Unassigned") 

trophic_colours <- c("Phototroph" = "#A3BE8C",
                     "Mixotroph" = "#5E81AC",
                     "Consumer" = "#B48EAD",
                     "Parasite/symbiont" = "#EBCB8B",
                     "Saprotroph" = "#BF616A",
                     "Unassigned" = "grey85")

ASV_rel_tab <- ASV_tab / rowSums(ASV_tab) # relative abundance, samples as rows
trophic_abund <- as.data.frame(t(rowsum(t(ASV_rel_tab), group = trophic$trophic_group[match(colnames(ASV_rel_tab), trophic$ASV)]))) # sum ASVs per trophic group

trophic_season <- trophic_abund %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "trophic_group", values_to = "Abundance") %>%
  mutate(season = meta[Sample, "season"]) %>%
  group_by(season, trophic_group) %>%
  summarise(mean = mean(Abundance),
            sd = sd(Abundance), .groups = "drop") %>%
  mutate(season = factor(season, levels = c("Spring", "Autumn")),
         trophic_group = factor(trophic_group, levels = trophic_levels),
         label = if_else(mean >= 0.03, paste0(round(mean * 100), "%"), "")) # label groups making up at least 3%

write.csv(trophic_season, "output/community_structure/trophic_season_mean.csv")

trophic_season_plot <- ggplot(trophic_season, aes(x = season, y = mean, fill = trophic_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 1.6, colour = "grey15") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"),
        legend.position = "right") +
  labs(y = "Mean relative abundance", fill = "Trophic group") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 1) +
  scale_fill_manual(values = trophic_colours,
                    guide = guide_legend(ncol = 1))

trophic_season_plot

ggsave("output/community_structure/figures/trophic_season_bar.png", trophic_season_plot, width = 90, height = 90, units = "mm", dpi = 600)
ggsave("output/community_structure/figures/trophic_season_bar.pdf", trophic_season_plot, width = 90, height = 90, units = "mm", dpi = 600)


# ASVs driving the seasonal shift
n_top_ASVs <- 20 # number of ASVs shown

ASV_season <- as.data.frame(ASV_rel_tab) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(season = meta[Sample, "season"]) %>%
  group_by(ASV, season) %>%
  summarise(mean = mean(Abundance), .groups = "drop") %>% # mean across all samples in each season (zeros kept)
  pivot_wider(names_from = season, values_from = mean, names_prefix = "mean_") %>%
  mutate(diff_autumn_minus_spring = mean_Autumn - mean_Spring,
         log2_fold_change = log2((mean_Autumn + 1e-5) / (mean_Spring + 1e-5))) %>%
  left_join(tax_tab2 %>% select(taxon, Phylum, Genus, genus_label, ASV_label_italic), by = c("ASV" = "taxon")) %>%
  left_join(trophic %>% select(ASV, trophic_group), by = "ASV")

write.csv(ASV_season, "output/community_structure/ASV_season_difference.csv")

top_ASVs <- ASV_season %>%
  slice_max(abs(diff_autumn_minus_spring), n = n_top_ASVs) %>%
  mutate(ASV = fct_reorder(ASV, diff_autumn_minus_spring),
         Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others"))

ASV_bar_plot <- ggplot(top_ASVs, aes(x = diff_autumn_minus_spring * 100, y = ASV, fill = Phylum)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, size = 0.3, colour = "black") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.key.size = unit(0.6, "lines")) +
  labs(x = "Change in mean relative abundance (%)\n(spring | autumn)", fill = "Phylum") +
  scale_y_discrete(labels = function(x) parse(text = top_ASVs$ASV_label_italic[match(x, top_ASVs$ASV)])) + # ASV: genus (italic) or last assigned rank
  scale_fill_manual(values = phylum_colours,
                    labels = function(x) ifelse(x == "x__Others", "Others", x))

ASV_bar_plot


# genera driving the seasonal shift
n_top_genera <- 20 # number of genera shown

genus_names <- tax_tab2[colnames(ASV_rel_tab), "Genus"]
keep <- tax_tab2[colnames(ASV_rel_tab), "is_genus"] # true genera only (unassigned genera and class/order/clade names left out, as the ordination wascores)

genus_abund <- as.data.frame(t(rowsum(t(ASV_rel_tab[, keep]), group = genus_names[keep]))) # sum ASVs per genus

genus_season <- genus_abund %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "Genus", values_to = "Abundance") %>%
  mutate(season = meta[Sample, "season"]) %>%
  group_by(Genus, season) %>%
  summarise(mean = mean(Abundance), .groups = "drop") %>% # mean across all samples in each season (zeros kept)
  pivot_wider(names_from = season, values_from = mean, names_prefix = "mean_") %>%
  mutate(diff_autumn_minus_spring = mean_Autumn - mean_Spring,
         log2_fold_change = log2((mean_Autumn + 1e-5) / (mean_Spring + 1e-5)),
         Phylum = sapply(Genus, function(x) names(sort(table(tax_tab2$Phylum[tax_tab2$is_genus & tax_tab2$Genus %in% x]), decreasing = T))[1])) # most common phylum of the genus' ASVs

write.csv(genus_season, "output/community_structure/genus_season_difference.csv")

top_genera <- genus_season %>%
  slice_max(abs(diff_autumn_minus_spring), n = n_top_genera) %>%
  mutate(Genus = fct_reorder(Genus, diff_autumn_minus_spring),
         Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others"))

genus_bar_plot <- ggplot(top_genera, aes(x = diff_autumn_minus_spring * 100, y = Genus, fill = Phylum)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, size = 0.3, colour = "black") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black", face = "italic"), # genus names in italics
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.key.size = unit(0.6, "lines")) +
  labs(x = "Change in mean relative abundance (%)\n(spring | autumn)", fill = "Phylum") +
  scale_fill_manual(values = phylum_colours,
                    labels = function(x) ifelse(x == "x__Others", "Others", x))

genus_bar_plot

ggsave("output/community_structure/figures/genus_season_bar.png", genus_bar_plot, width = 100, height = 90, units = "mm", dpi = 600)
ggsave("output/community_structure/figures/genus_season_bar.pdf", genus_bar_plot, width = 100, height = 90, units = "mm", dpi = 600)


# diversity
diversity <- estimate_richness(ps) # calculate diversity
rownames(diversity) <- gsub("\\.","-",rownames(diversity)) # fix rownames

write.csv(diversity, "output/community_structure/diversity.csv")

meta <- rownames_to_column(meta, var="18S_ID")
diversity <- rownames_to_column(diversity, var="18S_ID")
diversity_meta <- merge(meta, diversity, by = "18S_ID") # add metadata

diversity_seasonal <- diversity_meta %>% 
  select(season, Sample_ID, Shannon, Observed) %>% 
  group_by(season) %>%
  summarise(mean_diversity = mean(Shannon),
            SD_diversity = sd(Shannon),
            mean_richness = mean(Observed),
            SD_richness = sd(Observed))

diversity_meta$season <- factor(diversity_meta$season, levels = c("Spring", "Autumn")) # set order of seasons

diversity_plot <- ggplot(diversity_meta, aes(x=season, y=Shannon, fill=season))+
  geom_boxplot(outlier.size = 0.5, alpha=0.8, size = 0.2)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none") +
  # facet_wrap(~ Year) +
  labs(y="Shannon diversity", title="Diversity")+
  scale_y_continuous(expand=c(0,0),limits=c(0,6),breaks=seq(0,6,by=2)) +
  scale_fill_manual(values=c("forestgreen", "steelblue"))

diversity_plot

richness_plot <- ggplot(diversity_meta, aes(x=season, y=Observed, fill=season))+
  geom_boxplot(outlier.size = 0.5, alpha=0.8, size = 0.2)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none") +
  # facet_wrap(~ Year) +
  labs(y="ASV richness", title="Richness")+
  scale_y_continuous(expand=c(0,0),limits=c(0,750),breaks=seq(0,750,by=200)) +
  scale_fill_manual(values=c("forestgreen", "steelblue"))

richness_plot

diversity_richness_plot <- diversity_plot + richness_plot +
  plot_layout(axes = "collect", ncol=2) &
  theme(plot.title = element_blank())

diversity_richness_plot

# diversity and richness anova
diversity_meta$Year <- factor(diversity_meta$Year)
diversity_meta$season <- factor(diversity_meta$season)

anova_model <- aov(Shannon ~ Year * season, data = diversity_meta)
summary(anova_model)

anova_model <- aov(Observed ~ Year * season, data = diversity_meta)
summary(anova_model)

# linear mixed models
options(contrasts = c("contr.sum", "contr.poly")) # sum-to-zero contrasts so the main effects of season/year are tested correctly when the interaction is also in the model (type III tests)

diversity_meta$Year <- factor(diversity_meta$Year)

shannon_model <- lmer(Shannon ~ season * Year + (1 | RSN_ID), data = diversity_meta) # shannon: continuous and roughly normal -> linear mixed model
anova(shannon_model) # F tests for season, Year and season:Year

richness_model <- glmer.nb(Observed ~ season * Year + (1 | RSN_ID), data = diversity_meta) # richness: counts, variance larger than the mean -> negative binomial GLMM
car::Anova(richness_model, type = 3) # Wald chi-square tests
emmeans::emmeans(richness_model, pairwise ~ Year)


# core community
ps_core <- core(ps_rel, detection = 0.0001, prevalence = 0.5) # core if relative abundance >detection in >prevalence % of samples
core_tax <- as.data.frame(tax_table(ps_core)) # extract core taxa
core_rel <- t(as.data.frame(otu_table(ps_core))) # extract core relative abundance

phylum_counts <- core_tax %>%
  count(Phylum, name = "count")

detections <- c(0.0001, 0.001, 0.01)
prevalences <- seq(0.05, 1, by = 0.05) # 0.5 = 50% prevalence

core_df <- lapply(detections, function(thresh) {
  prevalence <- rowMeans(core_rel >= thresh)
  data.frame(ASV = rownames(core_rel),
             detection = thresh,
             prevalence = prevalence)
}) %>%
  bind_rows() %>% 
  merge(., core_tax, by.x = "ASV", by.y = "row.names")

write.csv(core_df, "output/community_structure/core.csv")

ASV_order <- core_df %>%
  group_by(ASV) %>%
  summarise(mean_prev = mean(prevalence)) %>%
  arrange(mean_prev) %>%
  pull(ASV)

core_df$ASV <- factor(core_df$ASV, levels = ASV_order)

core_plot <- ggplot(core_df, aes(x = detection, y = ASV, fill = prevalence)) +
  geom_tile() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.8, "lines"),
        legend.position = "right",
        axis.text = element_text(size = 4.5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black")) +
  labs(x = "Detection threshold (relative abundance)") +
  scale_x_log10(breaks = detections) +
  scale_fill_gradientn(colours = rev(brewer.pal(9, "RdYlBu")), 
                       name = "Prevalence (%)",
                       limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                       labels = function(x) formatC(x * 100, format = "f", digits = 0)) +
  scale_y_discrete(labels = function(x) parse(text = core_df$ASV_label_italic[match(x, core_df$taxon)]))

core_plot

ggsave("output/community_structure/figures/core_heatmap.png", core_plot, width = 90, height = 90, units = "mm", dpi = 600)
ggsave("output/community_structure/figures/core_heatmap.pdf", core_plot, width = 90, height = 90, units = "mm", dpi = 600)


# abundance vs prevalence
ASV_rel <- ASV_tab %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(across(everything(), ~ . / sum(.))) %>% # convert to relative abundance
  t(.) # samples as rows

mean_abundance <- as.data.frame(colMeans(ASV_rel)) # calculate mean across all samples
mean_abundance <- mean_abundance %>% rownames_to_column(., var = "ASV")
colnames(mean_abundance)[2] <- "mean_rel_abund"

presence_absence <- (ASV_rel > 0) * 1 # convert coverage to binary presence/absence
occupancy <- as.data.frame(colSums(presence_absence)) # sum number of samples each ASV had abundance > 0 in
occupancy <- occupancy %>% rownames_to_column(., var = "ASV")
colnames(occupancy)[2] <- "occupancy"

total_samples <- nrow(ASV_rel) # calculate number of samples

occupancy$percentage_occupancy <- (occupancy$occupancy / total_samples) * 100 # convert to percentage of total samples

occupancy_abundance <- merge(occupancy, mean_abundance, by="ASV")
occupancy_tax <- merge(occupancy_abundance, tax_tab2, by.x="ASV", by.y="row.names")

scientific_10_labels <- function(x) {
  parse(text = paste0("10^", floor(log10(abs(x))))) # convert labels to scientific notation
}

occupancy_tax2 <- occupancy_tax %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla$Phylum, Phylum, "x__Others"))

occupancy_tax2$core <- ifelse(occupancy_tax2$ASV %in% core_tax$taxon, "Core", "Not core")

abundance_occupancy_log <- ggplot(occupancy_tax2, aes(x = percentage_occupancy, y = mean_rel_abund, colour = Phylum)) +
  geom_point(size = 0.8, alpha = 0.8, aes(shape = core)) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_blank(), # element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 4)) +
  labs(x = "Prevalence (%)", y = "Mean relative abundance") +
  scale_y_log10(labels = scientific_10_labels) +
  scale_colour_manual(values = phylum_colours,
                      labels = function(x) ifelse(x == "x__Others", "Others", x)) +
  scale_shape_manual(name = "Core ASV",
                     guide = guide_legend(ncol = 2),
                     values = c("Not core" = 16, "Core" = 17)) +
  geom_segment(aes(x = 50, xend = 50, y = 0.0001, yend = Inf),
               linetype = "dashed", color = "black", size = 0.2) +
  geom_segment(aes(x = 50, xend = Inf, y = 0.0001, yend = 0.0001),
               linetype = "dashed", color = "black", size = 0.2)

abundance_occupancy_log

mean(occupancy$percentage_occupancy < 10) * 100
mean(occupancy$percentage_occupancy < 1) * 100

no_y_axis <- theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank()) # y axis shared with the composition plot


# combined figure
legend_composition_h <- ggpubr::get_legend(composition_plot +
                                             guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
                                             theme(legend.direction = "horizontal",
                                                   legend.title.position = "top")) 
legend_composition_h <- patchwork::wrap_elements(full = legend_composition_h)

abundance_occupancy_inset <- abundance_occupancy_log +
  guides(colour = "none") + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.04), 
        legend.justification = c(1, 0),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank())

layout <- "
AAAAAAAAAABB
CCCCCCCEEEEE
DDDDDDDEEEEE
FFFFFFFEEEEE
"

combined_plot <- (composition_plot + theme(legend.position = "none")) +
  (phylum_season_plot + no_y_axis + theme(legend.position = "none")) +
  abundance_occupancy_inset +
  (diversity_richness_plot & theme(legend.position = "none")) +
  (genus_bar_plot + theme(legend.position = "none")) +
  legend_composition_h +
  plot_layout(design = layout, heights = c(1.5, 0.9, 0.9, 0.3))

combined_plot

ggsave("output/community_structure/figures/community_combined.png", combined_plot, width = 180, height = 190, units = "mm", dpi = 800)
ggsave("output/community_structure/figures/community_combined.pdf", combined_plot, width = 180, height = 190, units = "mm", dpi = 400)

