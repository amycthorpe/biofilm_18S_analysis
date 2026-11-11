# biofilm analysis - 18S
# NMDS and distance decay


# load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(geosphere)
library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)


set.seed(123)
options(scipen = 999)


# parameters
nmds_trymax <- 20 # random starts for metaMDS
n_perm <- 999 # permutations for envfit, adonis and Mantel tests

env_vars <- c("Alkalinity_combined", "Ammoniacal_N", "DOC", "Chloride", "Conductivity", "Nitrate_N", "Nitrite_N",
              "Orthophosphate", "DO", "SiO2", "Water_temp", "pH")
vector_p_max <- 0.05 # only draw water chemistry arrows with p <= threshold
vector_r2_min <- 0 # and R2 >= threshold
arrow_fill <- 0.9 # longest arrow as a proportion of the furthest sample from the origin

log_vars <- c("Ammoniacal_N", "DOC", "Chloride", "Conductivity", "Nitrate_N", "Nitrite_N", "Orthophosphate", "SiO2") # log(x + 1) (skewed)

top_n <- 20 # most abundant
wa_expand <- TRUE # expand weighted averages to the spread of the samples (otherwise widespread genera all sit near the origin)

distance_bin_km <- 25 # bin width for geographic distance decay
env_bin_width <- 0.5 # bin width for water chemistry distance decay
min_pairs_bin <- 50 # bins with fewer pairs than this are dropped from the plots

phototroph_colour <- "forestgreen" # phototroph-dominated (log ratio > 0)
heterotroph_colour <- "firebrick" # heterotroph-dominated (log ratio < 0)
ratio_pseudo <- 0.001 # added to both sides before log ratio (avoids log(0))
ratio_quantile <- 0.02 # colour scale squished to this and 1 - this quantile

season_colours <- c("Spring" = "forestgreen", "Autumn" = "steelblue")
geo_colour <- "#D08770"
env_colour <- "#D08770"
zero_line_colour <- "grey70"

new_names <- c("Alkalinity_combined" = "Alkalinity",
               "Ammoniacal_N" = "Ammoniacal-N",
               "Nitrate_N" = "Nitrate-N",
               "Nitrite_N" = "Nitrite-N",
               "Water_temp" = "Temperature") # tidy variable names for plots


# import data
ASV_tab <- data.frame(fread("data/processed_2021_23/ASV_counts_18S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.csv("data/processed_2021_23/taxonomy_clean.csv", header=T, row.names=1, check.names=F)
meta_all <- read.table("data/metadata/seq_metadata_all_amplicons.csv", header=T, check.names=F, sep=",")
trophic <- read.csv("data/processed_2021_23/taxonomy_labels.csv", header=T, row.names=1, check.names=F)
phylum_colours <- read.csv("output/community_structure/phylum_colour_map.csv", row.names=1)

dir.create("output/ordination/figures", recursive = TRUE, showWarnings = FALSE)

tax_tab2 <- tax_tab # add true genus flag (genus names, not unassigned or clade names)
tax_tab2$is_genus <- trophic[rownames(tax_tab2), "is_genus"]

meta <- meta_all[meta_all$`18S_ID` %in% rownames(ASV_tab), ] # subset metadata to seq ID (change depending on amplicon)
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="18S_ID") # make seq ID rownames

meta$Sample_date <- as.Date(meta$Sample_date, format = "%d/%m/%Y")  # change format to date

meta <- meta %>% # add season column
  mutate(
    season = case_when(
      month(Sample_date) %in% c(03, 04, 05, 06, 07, 08) ~ "Spring",
      month(Sample_date) %in% c(09, 10, 11, 12, 01, 02) ~ "Autumn"
    )
  )

meta <- select(meta, -c(1:6)) # remove unwanted ID cols


# convert to relative abundance
ASV_rel <- ASV_tab / rowSums(ASV_tab) # samples as rows


# import environmental metadata
water_chem <- read.csv("data/metadata/water_chemistry/wims_summaries_3months.csv", check.names=F) %>%
  select(2, alkalinity_combined, ends_with("_mean"))

chem_names <- c(
  "Extended_ID",
  "Alkalinity_combined",
  "Alkalinity_grans",
  "Alkalinity",
  "Ammonia_N_unionised",
  "Ammoniacal_N",
  "DOC",
  "Chloride",
  "Conductivity",
  "Nitrate_N",
  "Nitrite_N",
  "TN_oxidised",
  "TN",
  "Orthophosphate",
  "DO",
  "DO_sat",
  "TP",
  "SiO2",
  "Water_temp",
  "pH")

colnames(water_chem) <- chem_names # tidy colnames

water_chem <- water_chem %>% select(-TN) # TN has lots of NAs so remove

seq_ids <- meta_all %>%
  select(1:4)

water_chem <- water_chem %>%
  left_join(seq_ids, by=c("Extended_ID"="16S_ID")) %>%
  filter(`18S_ID` %in% rownames(ASV_tab)) %>% # subset to same samples
  column_to_rownames(., var="18S_ID") %>%
  select(-c("Extended_ID", "ITS_ID", "rbcL_ID"))

meta_merged <- merge(meta, water_chem, by="row.names") %>%
  column_to_rownames(., var="Row.names")

colnames(meta_merged) <- gsub("/", "_", colnames(meta_merged))

meta_merged <- meta_merged[match(rownames(ASV_tab), rownames(meta_merged)),] # same sample order as ASV table
meta_merged$season <- factor(meta_merged$season, levels = c("Spring", "Autumn"))
meta_merged$season_year <- paste(meta_merged$season, meta_merged$Year)

# write.csv(meta_merged, "data/metadata/meta_merged.csv")


# NMDS
ASV_hel <- decostand(ASV_tab, method = "hellinger")
nmds <- metaMDS(ASV_hel, distance = "bray", autotransform = F, trymax = nmds_trymax)
nmds

nmds_scores <- as.data.frame(scores(nmds, display = "sites")) # extract scores
write.csv(nmds_scores, "output/ordination/nmds_scores.csv")

nmds_scores_meta <- merge(nmds_scores, meta_merged, by="row.names") %>% # combine with meta
  column_to_rownames(.,"Row.names")


# envfit (water chemistry)
formula <- as.formula(paste("nmds ~", paste(env_vars, collapse="+"))) # set envfit formula

nmds.fit <- envfit(formula,
                   data=meta_merged,
                   perm=n_perm,
                   na.rm=T)

nmds.fit

df <- nmds.fit$vectors # get dataframe of vector r and p values
nmds.fit_df <- cbind(as.data.frame(df[["arrows"]]),
                     rvalue = df[["r"]],
                     pvalue = df[["pvals"]])

nmds.fit_df <- nmds.fit_df %>% # add significance levels
  mutate(significance = case_when(
    pvalue > 0.05 ~ "ns",
    pvalue <= 0.05 & pvalue > 0.01 ~ "*",
    pvalue <= 0.01 & pvalue > 0.001 ~ "**",
    pvalue <= 0.001 ~ "***"
  )) %>%
  arrange(desc(rvalue))

write.csv(nmds.fit_df, "output/ordination/envfit.csv")

vectors <- as.data.frame(scores(nmds.fit, "vectors")) # set up vector arrows
vectors <- vectors * arrow_fill * max(abs(nmds_scores)) / max(abs(vectors)) # scale arrows to the plot
vectors <- vectors[rownames(nmds.fit_df)[nmds.fit_df$pvalue <= vector_p_max & nmds.fit_df$rvalue >= vector_r2_min], ] # filter arrows
arrowhead = arrow(length = unit(0.02, "npc"))

vectors2 <- vectors
matched <- intersect(rownames(vectors2), names(new_names))
rownames(vectors2)[rownames(vectors2) %in% matched] <- new_names[matched]


# NMDS coloured by season with water chemistry vectors
plot_all <- ggplot(nmds_scores_meta, aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = season)) +
  geom_segment(data = vectors2, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               size = 0.3, arrow = arrowhead, colour = "black") +
  geom_text_repel(data = vectors2, aes(x = NMDS1, y = NMDS2),
                  label = row.names(vectors2), size = 2, colour = "black",
                  max.overlaps = Inf, segment.colour = NA,
                  bg.colour = "white", bg.r = 0.15, force = 2, box.padding = 0.3) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = "bottom",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_colour_manual(name = "Season", values = season_colours) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste("Stress =", round(nmds$stress, 2)),
           size = 2)

plot_all


# adonis (season and year, permuted within site so repeat visits aren't treated as independent)
bray <- vegdist(ASV_hel, method = "bray")

meta_adonis <- meta_merged %>%
  mutate(Year = factor(Year), RSN_ID = factor(RSN_ID))

adonis_result <- adonis2(bray ~ Year * season,
                         data = meta_adonis,
                         permutations = how(nperm = n_perm, blocks = meta_adonis$RSN_ID),
                         by = "terms") # terms in order: Year, season, Year:season

adonis_result

write.csv(as.data.frame(adonis_result), "output/ordination/adonis_season_year.csv")


# top genera weighted-average positions
genus_names <- tax_tab2[colnames(ASV_rel), "Genus"]
keep <- tax_tab2[colnames(ASV_rel), "is_genus"] # true genera only (unassigned genera and class/order/clade names left out)

genus_tab <- t(rowsum(t(ASV_rel[, keep]), genus_names[keep])) 
top_genera <- head(names(sort(colMeans(genus_tab), decreasing = T)), top_n)

genus_wa <- as.data.frame(wascores(as.matrix(nmds_scores[, c("NMDS1", "NMDS2")]), genus_tab[rownames(nmds_scores), top_genera], expand = wa_expand)) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(mean_abund = colMeans(genus_tab[, Genus]), 
         Phylum = sapply(Genus, function(x) names(sort(table(tax_tab2$Phylum[tax_tab2$is_genus & tax_tab2$Genus %in% x]), decreasing = T))[1])) 

colours <- setNames(phylum_colours$colour, rownames(phylum_colours)) 
genus_wa <- genus_wa %>%
  mutate(Phylum_colour = ifelse(Phylum %in% names(colours), Phylum, "x__Others"))

write.csv(genus_wa, "output/ordination/genus_weighted_average_scores.csv")

genus_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1, alpha = 0.4, shape = 16, colour = "grey85") +
  geom_point(data = genus_wa, aes(x = NMDS1, y = NMDS2, size = mean_abund, fill = Phylum_colour),
             shape = 21, colour = "black", stroke = 0.3, alpha = 0.9) +
  geom_text_repel(data = genus_wa, aes(x = NMDS1, y = NMDS2, label = Genus),
                  size = 1.8, fontface = "italic", colour = "black", max.overlaps = Inf,
                  segment.colour = "grey40", segment.size = 0.2, min.segment.length = 0.1,
                  bg.colour = "white", bg.r = 0.15, box.padding = 0.4, point.padding = 0.2, force = 2) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        legend.key.size = unit(0.6, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_fill_manual(name = "Phylum", values = colours, labels = function(x) gsub("x__", "", x)) + 
  scale_size_continuous(name = "Mean relative abundance", range = c(1.5, 5)) +
  guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(size = 2.5)),
         size = guide_legend(nrow = 1, order = 2))

genus_plot


# top ASVs weighted-average positions
asv_labels <- trophic[colnames(ASV_rel), "ASV_label_italic"] 
names(asv_labels) <- colnames(ASV_rel)

top_asvs <- head(names(sort(colMeans(ASV_rel), decreasing = T)), top_n)

asv_wa <- as.data.frame(wascores(as.matrix(nmds_scores[, c("NMDS1", "NMDS2")]), as.matrix(ASV_rel[rownames(nmds_scores), top_asvs]), expand = wa_expand)) %>%
  rownames_to_column(var = "ASV") %>%
  mutate(mean_abund = colMeans(ASV_rel[, ASV]),
         label = asv_labels[ASV],
         Phylum = tax_tab2[ASV, "Phylum"],
         Phylum_colour = ifelse(Phylum %in% names(colours), Phylum, "x__Others"))

write.csv(asv_wa, "output/ordination/ASV_weighted_average_scores.csv")

asv_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1, alpha = 0.4, shape = 16, colour = "grey85") +
  geom_point(data = asv_wa, aes(x = NMDS1, y = NMDS2, size = mean_abund, fill = Phylum_colour),
             shape = 21, colour = "black", stroke = 0.3, alpha = 0.9) +
  geom_text_repel(data = asv_wa, aes(x = NMDS1, y = NMDS2, label = label), parse = TRUE,
                  size = 1.8, colour = "black", max.overlaps = Inf,
                  segment.colour = "grey40", segment.size = 0.2, min.segment.length = 0.1,
                  bg.colour = "white", bg.r = 0.15, box.padding = 0.7, point.padding = 0.3, force = 3, max.time = 3, max.iter = 50000, seed = 5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        legend.key.size = unit(0.6, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_fill_manual(name = "Phylum", values = colours, labels = function(x) gsub("x__", "", x)) +
  scale_size_continuous(name = "Mean relative abundance", range = c(1.5, 5)) +
  guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(size = 2.5)),
         size = guide_legend(nrow = 1, order = 2))

asv_plot


# distance decay
env <- meta_merged %>%
  select(RSN_ID, season_year, lat, long, all_of(env_vars)) %>%
  mutate(across(all_of(log_vars), ~ log(. + 1))) %>%
  na.omit()

nrow(env) # samples with complete water chemistry
table(env$season_year)

geo_km <- distm(env[, c("long", "lat")], fun = distHaversine) / 1000 # great-circle distance between sites (km)
env_dist <- as.matrix(dist(scale(env[, env_vars]))) # water chemistry distance: Euclidean on scaled variables
bc <- as.matrix(bray)[rownames(env), rownames(env)] # Bray-Curtis on Hellinger
rownames(geo_km) <- colnames(geo_km) <- rownames(env)
rownames(env_dist) <- colnames(env_dist) <- rownames(env)

surveys <- sort(unique(env$season_year))

pairs_df <- data.frame()
mantel_results <- data.frame()

for (sy in surveys) {
  ids <- rownames(env)[env$season_year == sy]
  lower <- lower.tri(geo_km[ids, ids])

  pairs_df <- bind_rows(pairs_df,
                        data.frame(season_year = sy,
                                   geo_km = geo_km[ids, ids][lower],
                                   env_dist = env_dist[ids, ids][lower],
                                   similarity = 1 - bc[ids, ids][lower]))

  bc_d <- as.dist(bc[ids, ids])
  geo_d <- as.dist(geo_km[ids, ids])
  env_d <- as.dist(env_dist[ids, ids])

  m_geo <- mantel(bc_d, geo_d, method = "spearman", permutations = n_perm)
  m_env <- mantel(bc_d, env_d, method = "spearman", permutations = n_perm)
  m_geo_partial <- mantel.partial(bc_d, geo_d, env_d, method = "spearman", permutations = n_perm) # geography, controlling for water chemistry
  m_env_partial <- mantel.partial(bc_d, env_d, geo_d, method = "spearman", permutations = n_perm) # water chemistry, controlling for geography

  mantel_results <- bind_rows(mantel_results,
                              data.frame(season_year = sy,
                                         n_samples = length(ids),
                                         test = c("geographic", "water chemistry", "geographic | water chemistry", "water chemistry | geographic"),
                                         r = c(m_geo$statistic, m_env$statistic, m_geo_partial$statistic, m_env_partial$statistic),
                                         p = c(m_geo$signif, m_env$signif, m_geo_partial$signif, m_env_partial$signif)))
}

nrow(pairs_df) # number of pairs

write.csv(mantel_results, "output/ordination/mantel_by_survey.csv")

mantel_summary <- mantel_results %>% # mean across year x season
  group_by(test) %>%
  summarise(mean_r = mean(r),
            min_r = min(r),
            max_r = max(r),
            n_significant = sum(p <= 0.05),
            n_surveys = n(),
            .groups = "drop")

mantel_season <- mantel_results %>% # mean by season
  mutate(season = word(season_year, 1)) %>%
  group_by(season, test) %>%
  summarise(mean_r = mean(r),
            n_significant = sum(p <= 0.05),
            .groups = "drop")

write.csv(mantel_summary, "output/ordination/mantel_summary.csv")
write.csv(mantel_season, "output/ordination/mantel_season_summary.csv")

mantel_summary
mantel_season

mantel_table <- mantel_results %>%
  mutate(season = factor(word(season_year, 1), levels = c("Spring", "Autumn")),
         Year = word(season_year, 2),
         value = paste0(sprintf("%.2f", r), ifelse(p <= 0.001, "***", ifelse(p <= 0.01, "**", ifelse(p <= 0.05, "*", " ns"))))) %>%
  select(season, Year, n_samples, test, value) %>%
  pivot_wider(names_from = test, values_from = value) %>%
  arrange(season, Year)

write.csv(mantel_table, "output/ordination/table_mantel_by_survey.csv", row.names = F)

mean_r <- setNames(mantel_summary$mean_r, mantel_summary$test) # for plot labels


# geographic distance decay
decay_binned <- pairs_df %>%
  mutate(bin = floor(geo_km / distance_bin_km) * distance_bin_km + distance_bin_km / 2) %>%
  group_by(bin) %>%
  summarise(mean_similarity = mean(similarity),
            sd_similarity = sd(similarity),
            n = n(),
            .groups = "drop") %>%
  filter(n >= min_pairs_bin) # drop sparse bins

write.csv(decay_binned, "output/ordination/distance_decay_geo_binned.csv")

geo_decay_plot <- ggplot(decay_binned, aes(x = bin, y = mean_similarity)) +
  geom_ribbon(aes(ymin = mean_similarity - sd_similarity, ymax = mean_similarity + sd_similarity), alpha = 0.2, fill = geo_colour) +
  geom_line(linewidth = 0.4, colour = geo_colour) +
  geom_point(size = 0.6, colour = geo_colour) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2)) +
  labs(x = "Geographic distance (km)", y = "Community similarity") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = 2,
           label = paste0("Mantel r = ", sprintf("%.2f", mean_r["geographic"]),
                          "\npartial r, controlling for chemistry = ", sprintf("%.2f", mean_r["geographic | water chemistry"])))

geo_decay_plot


# water chemistry distance decay
env_binned <- pairs_df %>%
  mutate(bin = floor(env_dist / env_bin_width) * env_bin_width + env_bin_width / 2) %>%
  group_by(bin) %>%
  summarise(mean_similarity = mean(similarity),
            sd_similarity = sd(similarity),
            n = n(),
            .groups = "drop") %>%
  filter(n >= min_pairs_bin) # drop sparse bins

write.csv(env_binned, "output/ordination/distance_decay_env_binned.csv")

env_decay_plot <- ggplot(env_binned, aes(x = bin, y = mean_similarity)) +
  geom_ribbon(aes(ymin = mean_similarity - sd_similarity, ymax = mean_similarity + sd_similarity), alpha = 0.2, fill = env_colour) +
  geom_line(linewidth = 0.4, colour = env_colour) +
  geom_point(size = 0.6, colour = env_colour) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2)) +
  labs(x = "Water chemistry distance", y = "Community similarity") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = 2,
           label = paste0("Mantel r = ", sprintf("%.2f", mean_r["water chemistry"]),
                          "\npartial r, controlling for geography = ", sprintf("%.2f", mean_r["water chemistry | geographic"])))

env_decay_plot


# distance decay plots on the same scale
decay_y_limits <- range(c(decay_binned$mean_similarity - decay_binned$sd_similarity,
                          decay_binned$mean_similarity + decay_binned$sd_similarity,
                          env_binned$mean_similarity - env_binned$sd_similarity,
                          env_binned$mean_similarity + env_binned$sd_similarity))

geo_decay_plot <- geo_decay_plot + scale_y_continuous(limits = decay_y_limits)
env_decay_plot <- env_decay_plot + scale_y_continuous(limits = decay_y_limits)


# NMDS on the same axes
nmds_x_limits <- range(c(nmds_scores$NMDS1, vectors2$NMDS1, genus_wa$NMDS1))
nmds_y_limits <- range(c(nmds_scores$NMDS2, vectors2$NMDS2, genus_wa$NMDS2))

plot_all_shared <- plot_all + coord_cartesian(xlim = nmds_x_limits, ylim = nmds_y_limits)
genus_plot_shared <- genus_plot + coord_cartesian(xlim = nmds_x_limits, ylim = nmds_y_limits)


# main figures: NMDS [A][B] and distance decay [C][D] as separate files
fig_AB <- ((plot_all_shared + theme(aspect.ratio = 1)) | (genus_plot_shared +theme(aspect.ratio = 1))) +
  plot_layout(axes = "collect", axis_titles = "collect")

ggsave("output/ordination/figures/NMDS.png", fig_AB, width = 180, height = 110, units = "mm", dpi = 600)
ggsave("output/ordination/figures/NMDS.pdf", fig_AB, width = 180, height = 110, units = "mm", dpi = 600)

fig_CD <- (geo_decay_plot | env_decay_plot) +
  plot_layout(axes = "collect_y", axis_titles = "collect_y")

ggsave("output/ordination/figures/distance_decay.png", fig_CD, width = 180, height = 60, units = "mm", dpi = 600)
ggsave("output/ordination/figures/distance_decay.pdf", fig_CD, width = 180, height = 60, units = "mm", dpi = 600)

fig_ABCD <- ((plot_all_shared + theme(aspect.ratio = 1)) |
               (genus_plot_shared + theme(aspect.ratio = 1))) /
  (geo_decay_plot | env_decay_plot) +
  plot_layout(heights = c(1, 0.5),
              axes = "collect", axis_titles = "collect")

ggsave("output/ordination/figures/NMDS_distance_decay.png", fig_ABCD, width = 180, height = 160, units = "mm", dpi = 600)
ggsave("output/ordination/figures/NMDS_distance_decay.pdf", fig_ABCD, width = 180, height = 160, units = "mm", dpi = 600)


# NMDS coloured by water chemistry (log scale for skewed variables)
nmds_alk <- ggplot(nmds_scores_meta %>% arrange(!is.na(Alkalinity_combined)), aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = log(Alkalinity_combined + 1))) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.5, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2", title = "Alkalinity") +
  scale_colour_viridis_c(name = "log", option = "viridis", na.value = "grey80") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("R^2 == ", round(nmds.fit_df["Alkalinity_combined", "rvalue"], 2), " * ', ' ~ '", nmds.fit_df["Alkalinity_combined", "significance"], "'"), parse = TRUE,
           size = 2)

nmds_alk

nmds_DO <- ggplot(nmds_scores_meta %>% arrange(!is.na(DO)), aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = DO)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.5, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2", title = "DO") +
  scale_colour_viridis_c(name = "mg L⁻¹", option = "viridis", na.value = "grey80") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("R^2 == ", round(nmds.fit_df["DO", "rvalue"], 2), " * ', ' ~ '", nmds.fit_df["DO", "significance"], "'"), parse = TRUE,
           size = 2)

nmds_DO

nmds_nitrate <- ggplot(nmds_scores_meta %>% arrange(!is.na(Nitrate_N)), aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = log(Nitrate_N + 1))) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.5, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2", title = "Nitrate-N") +
  scale_colour_viridis_c(name = "log", option = "viridis", na.value = "grey80") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("R^2 == ", round(nmds.fit_df["Nitrate_N", "rvalue"], 2), " * ', ' ~ '", nmds.fit_df["Nitrate_N", "significance"], "'"), parse = TRUE,
           size = 2)

nmds_nitrate

nmds_temp <- ggplot(nmds_scores_meta %>% arrange(!is.na(Water_temp)), aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1, alpha = 0.6, shape = 16, aes(colour = Water_temp)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.5, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2", title = "Temperature") +
  scale_colour_viridis_c(name = "°C", option = "viridis", na.value = "grey80") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("R^2 == ", round(nmds.fit_df["Water_temp", "rvalue"], 2), " * ', ' ~ '", nmds.fit_df["Water_temp", "significance"], "'"), parse = TRUE,
           size = 2)

nmds_temp

supp_chem <- (nmds_alk | nmds_DO) / (nmds_nitrate | nmds_temp) +
  plot_layout(axes = "collect", axis_titles = "collect")

supp_chem

ggsave("output/ordination/figures/supp_NMDS_chemistry.png", supp_chem, width = 180, height = 180, units = "mm", dpi = 600)
ggsave("output/ordination/figures/supp_NMDS_chemistry.pdf", supp_chem, width = 180, height = 180, units = "mm", dpi = 600)


# phototroph:heterotroph ratio
write.csv(as.data.frame(table(trophic_group = trophic$trophic_group, trophic_broad = trophic$trophic_broad)) %>% filter(Freq > 0),
          "output/ordination/phototroph_heterotroph_categories.csv", row.names = F) # check which groups went where (ASV counts)

PH_tab <- t(rowsum(t(ASV_rel), as.character(trophic[colnames(ASV_rel), "trophic_broad"]))) %>% 
  as.data.frame()

ratio_scores <- nmds_scores_meta %>%
  select(NMDS1, NMDS2, season) %>%
  mutate(phototroph = PH_tab[rownames(.), "Phototroph"],
         heterotroph = PH_tab[rownames(.), "Heterotroph"],
         log_PH = log10((phototroph + ratio_pseudo) / (heterotroph + ratio_pseudo))) # > 0 phototroph-dominated, < 0 heterotroph-dominated

write.csv(ratio_scores, "output/ordination/phototroph_heterotroph_ratio.csv")

ratio_fit <- envfit(nmds_scores[rownames(ratio_scores), c("NMDS1", "NMDS2")], ratio_scores[, "log_PH", drop = F], permutations = n_perm) 
ratio_fit

ratio_summary <- ratio_scores %>%
  group_by(season) %>%
  summarise(mean_pct_phototroph = mean(phototroph) * 100,
            mean_pct_heterotroph = mean(heterotroph) * 100,
            median_log_PH = median(log_PH),
            prop_phototroph_dominant = mean(log_PH > 0),
            .groups = "drop") %>%
  mutate(wilcox_p = wilcox.test(ratio_scores$log_PH[ratio_scores$season == "Spring"], ratio_scores$log_PH[ratio_scores$season == "Autumn"])$p.value,
         envfit_r2 = ratio_fit$vectors$r,
         envfit_p = ratio_fit$vectors$pvals)

write.csv(ratio_summary, "output/ordination/phototroph_heterotroph_summary.csv")
ratio_summary

ratio_limit <- max(abs(quantile(ratio_scores$log_PH, c(ratio_quantile, 1 - ratio_quantile))))

ratio_plot <- ggplot(ratio_scores %>% arrange(abs(log_PH)), aes(x = NMDS1, y = NMDS2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = zero_line_colour, size = 0.3) +
  geom_point(size = 1.1, shape = 21, colour = "grey60", stroke = 0.1, aes(fill = log_PH)) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.position = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.5, "lines"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6, colour = "black")) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_fill_gradient2(name = "log10\n(phototroph / heterotroph)",
                       low = heterotroph_colour, mid = "white", high = phototroph_colour, midpoint = 0,
                       limits = c(-ratio_limit, ratio_limit), oob = scales::squish) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
         label = paste0("R^2 == ", round(ratio_fit$vectors$r, 2),
                        " * ', ' ~ '", case_when(ratio_fit$vectors$pvals > 0.05 ~ "ns",
                                                 ratio_fit$vectors$pvals > 0.01 ~ "*",
                                                 ratio_fit$vectors$pvals > 0.001 ~ "**",
                                                 TRUE ~ "***"), "'"), parse = TRUE, size = 2)

ratio_plot

ggsave("output/ordination/figures/supp_NMDS_phototroph_heterotroph.png", ratio_plot, height = 90, width = 90, units = "mm", dpi = 600)
ggsave("output/ordination/figures/supp_NMDS_phototroph_heterotroph.pdf", ratio_plot, height = 90, width = 90, units = "mm", dpi = 600)

