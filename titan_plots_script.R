# biofilm analysis for manuscripts - 18S
# TITAN analysis


# load libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(ggridges)
library(patchwork)
library(ggvenn)
library(glue)


set.seed(123)
options(scipen = 999)


# set thresholds
purity_cut <- 0.95 # TITAN purity and reliability
reliability_cut <- 0.95
z_cut <- 10 # minimum median z score (magnitude of response)
abund_cut <- 0.001 # minimum mean relative abundance

titan_files <- c(alkalinity = "Alkalinity_combined", nitrate = "Nitrate_N", DO = "DO", temp = "Water_temp") # variable name = file name in data/titan

out_dir <- "output/titan"
dir.create(file.path(out_dir, "figures"), recursive = T, showWarnings = F)


# import data
ASV_tab <- data.frame(fread("data/processed_2021_23/ASV_counts_18S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.table("data/processed_2021_23/taxonomy_clean.csv", header=T, row.names=1, check.names=F, sep=",")
trophic <- read.csv("data/processed_2021_23/taxonomy_labels.csv", header=T, row.names=1, check.names=F)

ASV_rel <- ASV_tab / rowSums(ASV_tab) # samples as rows

titan_res <- map(titan_files, ~ readRDS(glue("data/titan/titan_{.x}_10perc.rds")))
map(titan_res, ~ .x$sumz.cp) # check change points


# titan results: species summary and bootstrapped change points for each variable
sppmax_all <- map_dfr(names(titan_res), function(v) {
  as_tibble(titan_res[[v]]$sppmax, rownames = "taxon") %>% # extract summary table and clean colnames
    rename(Q5 = `5%`, Q10 = `10%`, Q50 = `50%`, Q90 = `90%`, Q95 = `95%`) %>%
    mutate(variable = v)
})

results_df <- map_dfr(names(titan_res), function(v) {
  sppmax_v <- sppmax_all %>% filter(variable == v)
  n_taxa <- dim(titan_res[[v]]$metricArray)[1] # get taxa
  n_boot <- dim(titan_res[[v]]$metricArray)[3] # get boot
  map_dfr(1:n_taxa, function(i) {
    tibble(taxon = sppmax_v$taxon[i],
           chk_pts = titan_res[[v]]$metricArray[i, 2, ], # metric 2 = change points
           bootstrap = 1:n_boot)}) %>%
    left_join(sppmax_v, by = "taxon") # merge
})

sppmax_all <- sppmax_all %>%
  filter(purity >= purity_cut, reliability >= reliability_cut) # only pure and reliable taxa

results_filtered <- results_df %>%
  filter(chk_pts > Q5, chk_pts < Q95) %>% # filter bootstrap points to be within Q5 and Q95
  filter(purity >= purity_cut, reliability >= reliability_cut)

write.csv(sppmax_all, file.path(out_dir, "titan_all.csv"))


# venn diagram
taxon_list <- list(Alkalinity = unique(results_filtered$taxon[results_filtered$variable == "alkalinity"]),
                   'Nitrate-N' = unique(results_filtered$taxon[results_filtered$variable == "nitrate"]),
                   Temperature = unique(results_filtered$taxon[results_filtered$variable == "temp"]),
                   DO = unique(results_filtered$taxon[results_filtered$variable == "DO"]))

titan_venn <- ggvenn(taxon_list, fill_color = c("skyblue", "lightpink", "lightgreen", "orange"),
                     stroke_size = 0.4, set_name_size = 4, text_size = 2.8)

ggsave(file.path(out_dir, "figures/titan_venn.png"), titan_venn, width = 120, height = 90, units = "mm", dpi = 600)
ggsave(file.path(out_dir, "figures/titan_venn.pdf"), titan_venn, width = 120, height = 90, units = "mm", dpi = 600)


# abundance and occupancy
abund_occ <- data.frame(taxon = colnames(ASV_rel),
                        mean_rel_abund = colMeans(ASV_rel), # mean relative abundance
                        occupancy = colSums(ASV_rel != 0)) %>% # count of samples where abundance > 0
  mutate(percent_occupancy = occupancy / nrow(ASV_rel) * 100)


# with core community
core_df <- read.table("output/community_structure/core.csv", header=T, row.names=1, check.names=F, sep=",")


# with taxonomy, labels, abundance and core
final_df <- results_filtered %>%
  left_join(abund_occ, by = "taxon") %>%
  left_join(tax_tab %>% rownames_to_column("taxon"), by = "taxon") %>%
  left_join(trophic %>% rownames_to_column("taxon") %>% select(taxon, taxon_label = ASV_label_italic, trophic_group, trophic_broad), by = "taxon") %>%
  mutate(core = ifelse(taxon %in% core_df$taxon, "Core", "Not core"))

write.csv(final_df, file.path(out_dir, "titan_abund.csv"))


# filter ASVs
df_filtered <- final_df %>%
  filter(z.median > z_cut) %>% # filter based on z median (magnitude of response)
  filter(mean_rel_abund > abund_cut) %>% # filter based on mean relative abundance
  mutate(taxon = fct_reorder(taxon, mean_rel_abund)) # sort


# titan prep
labels <- c("alkalinity" = "Alkalinity", # set tidy labels for facets
            "DO" = "DO",
            "nitrate" = "Nitrate-N",
            "temp" = "Temperature")

change_points <- data.frame(variable = names(titan_res), # extract change points
                            fsumz_inc = map_dbl(titan_res, ~ .x$sumz.cp["fsumz+", "cp"]),
                            fsumz_dec = map_dbl(titan_res, ~ .x$sumz.cp["fsumz-", "cp"]))

grad_dec <- scales::col_numeric(palette = c("#d1eeea", "#2a5674"), domain = range(df_filtered$zscore)) # assign colour scales
grad_inc <- scales::col_numeric(palette = c("#fde0c5", "#eb4a40"), domain = range(df_filtered$zscore)) 

df_filtered <- df_filtered %>%
  mutate(fill_color = ifelse(maxgrp == 1,
                             grad_dec(zscore), 
                             grad_inc(zscore)))

write.csv(df_filtered, file.path(out_dir, "df_filtered.csv"))


# titan ridge plots
titan_alk_plot <- ggplot(data = df_filtered %>% filter(variable == "alkalinity"),
                         aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 3) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"), 
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("Alkalinity CP (mg "*L^{-1}*" of CaCO"[3]*")"), fill = "Z score", title = "Alkalinity") +
  geom_vline(data = change_points %>% filter(variable == "alkalinity"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "alkalinity"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) +
  scale_x_continuous(breaks = seq(0, 250, by = 50)) # range set with coord_cartesian in the combined plot (zooms without dropping bootstrap points)


titan_DO_plot <- ggplot(data = df_filtered %>% filter(variable == "DO"),
                        aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("DO CP (mg "*L^{-1}*")"), fill = "Z score", title = "DO") +
  geom_vline(data = change_points %>% filter(variable == "DO"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "DO"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) +
  scale_x_continuous(breaks = seq(8, 13, by = 1)) # range set with coord_cartesian in the combined plot (zooms without dropping bootstrap points)


titan_nitrate_plot <- ggplot(data = df_filtered %>% filter(variable == "nitrate"),
                             aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("Nitrate-N CP (mg "*L^{-1}*")"), fill = "Z score", title = "Nitrate-N") +
  geom_vline(data = change_points %>% filter(variable == "nitrate"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "nitrate"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) #+
  # scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2))


titan_temp_plot <- ggplot(data = df_filtered %>% filter(variable == "temp"),
                          aes(x = chk_pts, y = taxon, fill = fill_color)) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        plot.title = element_blank()) + # element_text(size = 6, colour = "black", hjust=0.5)) +
  scale_fill_identity(guide = "legend", 
                      breaks = c("#d1eeea", "#2a5674", "#fde0c5", "#eb4a40"),
                      labels = c("Low 1", "High 1", "Low 2", "High 2")) + 
  labs(x = expression("Temperature CP (°C)"), fill = "Z score", title = "Temperature") +
  geom_vline(data = change_points %>% filter(variable == "temp"), aes(xintercept = fsumz_inc), linetype = "dashed", color = "#eb4a40", size = 0.4) +
  geom_vline(data = change_points %>% filter(variable == "temp"), aes(xintercept = fsumz_dec), linetype = "dashed", color = "#2a5674", size = 0.4) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)])) #+
  # scale_x_continuous(limits = c(6, 16.5), breaks = seq(6, 16, by = 2))


legend_inc <- ggplot(data = subset(df_filtered, maxgrp == "2"), aes(x = chk_pts, y = taxon, fill = zscore)) +
  facet_wrap(~ variable) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme(legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_fill_gradient(low = "#fde0c5", high = "#eb4a40") +
  labs(fill = "Z score (tolerant)")
legend_inc <- patchwork::wrap_elements(full = cowplot::get_legend(legend_inc)) # extract legend

legend_dec <- ggplot(data = subset(df_filtered, maxgrp == "1"), aes(x = chk_pts, y = taxon, fill = zscore)) +
  facet_wrap(~ variable) +
  geom_density_ridges(alpha = 0.8, colour = NA, rel_min_height = 0.01, scale = 1) +
  theme(legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines")) +
  scale_fill_gradient(low = "#d1eeea", high = "#2a5674") +
  labs(fill = "Z score (sensitive)")
legend_dec <- patchwork::wrap_elements(full = cowplot::get_legend(legend_dec)) # extract legend


# abundance
abundance_plot <- ggplot(df_filtered, aes(x = "Abundance", y = taxon, fill = mean_rel_abund)) +
  geom_tile(color = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(fill = "Relative abundance", title = "Abundance") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

legend_abundance <- patchwork::wrap_elements(full = cowplot::get_legend(abundance_plot)) # extract legend


# occupancy
occupancy_plot <- ggplot(df_filtered, aes(x = "Prevalence", y = taxon, fill = percent_occupancy)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(0.8, "lines"),
        legend.key.height = unit(0.6, "lines"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  labs(fill = "Prevalence (%)", title = "Prevalence") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

legend_occupancy <- patchwork::wrap_elements(full = cowplot::get_legend(occupancy_plot)) # extract legend


# modules
role_shapes <- c("Module hubs" = 22, "Connectors" = 24) # peripheral nodes shown as the tile only
pal_vec <- c("#9C9EDEFF", "#DE9ED6FF", "#CEDB9CFF", "#E7969CFF", "#E7CB94FF", # module colours (same order as the network figure)
             "#6B6ECFFF", "#CE6DBDFF", "#B5CF6BFF", "#D6616BFF", "#E7BA52FF",
             "#5254A3FF", "#A55194FF", "#8CA252FF", "#AD494AFF", "#BD9E39FF",
             "#393B79FF", "#7B4173FF", "#637939FF", "#843C39FF", "#8C6D31FF")

seasonal_modules <- map_dfr(c("Spring", "Autumn"), function(s) {
  read.csv(file.path("output/network", paste0("modules_", s, ".csv")), header = T, check.names = F, row.names = 1) %>%
    select(taxon = name, module, taxa_roles) %>%
    mutate(season = s,
           module = paste0(substr(s, 1, 1), module)) # SM1 = spring module 1, AM1 = autumn module 1
})

module_colours <- setNames(pal_vec[1:10], paste0("M", 1:10)) # module colours follow the module number in both seasons (SM1 and AM1 share a colour)

df_modules <- df_filtered %>% # one row per plotted ASV and season
  distinct(taxon) %>%
  tidyr::crossing(season = c("Spring", "Autumn")) %>%
  left_join(seasonal_modules, by = c("taxon", "season")) %>%
  mutate(column = factor(if_else(season == "Spring", "Spring module", "Autumn module"), levels = c("Spring module", "Autumn module")),
         module = factor(sub("^[SA]", "", module), levels = names(module_colours)), # SM1/AM1 -> M1
         role = factor(if_else(taxa_roles %in% names(role_shapes), taxa_roles, NA_character_), levels = names(role_shapes)),
         taxon = factor(taxon, levels = levels(df_filtered$taxon)))

module_plot <- ggplot(df_modules, aes(x = column, y = taxon)) +
  geom_tile(aes(fill = module), colour = "white", linewidth = 0.4) +
  geom_point(data = df_modules %>% filter(!is.na(role)), aes(shape = role, fill = module, size = role), colour = "black", stroke = 0.4) + 
  scale_size_manual(values = c("Module hubs" = 1.2, "Connectors" = 0.9), guide = "none") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_manual(values = module_colours, limits = levels(droplevels(df_modules$module[!is.na(df_modules$module)])), na.value = "white", na.translate = FALSE,
                    guide = guide_legend(nrow = 2, order = 1, override.aes = list(shape = NA, alpha = 1))) +
  scale_shape_manual(values = role_shapes, labels = c("Hubs", "Connectors"), na.translate = FALSE, drop = FALSE,
                     guide = guide_legend(order = 2, override.aes = list(fill = "grey70", size = 2))) +
  labs(fill = "Module", shape = "Network role", title = "Module") +
  scale_x_discrete(position = "top", limits = c("Spring module", "Autumn module")) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

legend_module <- patchwork::wrap_elements(full = cowplot::get_legend(module_plot + guides(shape = "none"))) # extract legends
legend_role <- patchwork::wrap_elements(full = cowplot::get_legend(module_plot + guides(fill = "none")))


# taxonomy
phylum_colour_map <- read.table("output/community_structure/phylum_colour_map.csv", check.names=F, row.names=1, header=T, sep = ",")
phylum_colours <- setNames(phylum_colour_map$colour, phylum_colour_map$Phylum)

tax_plot <- ggplot(df_filtered, aes(x = "Phylum", y = taxon, fill = Phylum)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.direction = "vertical",
        legend.title.position = "top",
        legend.key.size = unit(0.6, "lines"), 
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  labs(fill = "Phylum", title = "Phylum") +
  scale_x_discrete(position = "top") +
  scale_fill_manual(guide = guide_legend(ncol = 3),
                    values = phylum_colours) +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))

legend_tax <- patchwork::wrap_elements(full = cowplot::get_legend(tax_plot)) # extract legend


# core
core_plot <- ggplot(df_filtered, aes(x = "Core", y = taxon)) +
  geom_point(data = df_filtered %>% filter(core == "Core"), shape = 42, size = 4) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust=0),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()) +
  labs(fill = "Core", title = "Core") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = FALSE, labels = function(x) parse(text = df_filtered$taxon_label[match(x, df_filtered$taxon)]))


# combined plot
# ridges taller than one row (scale > 1) stretch the ridge panels' y axis above the top ASV, which shifts their rows
# down relative to the labels and tiles, so all panels share the same y range
y_top <- max(map_dbl(list(titan_alk_plot, titan_DO_plot, titan_nitrate_plot, titan_temp_plot),
                     ~ ggplot_build(.x)$layout$panel_params[[1]]$y.range[2])) - 0.6 # top of the tallest ridge (minus the discrete axis padding)

# column titles for the tile panels sit just above the top ASV (not at the top of the taller ridge panels)
n_rows <- nlevels(df_filtered$taxon)
title_angle <- 90 # angle of the column titles

core_plot <- core_plot + theme(axis.text.x = element_blank()) + coord_cartesian(ylim = c(1, y_top), clip = "off") +
  annotate("text", x = "Core", y = n_rows + 0.8, label = "Core", angle = title_angle, hjust = 0, size = 6 / .pt)
tax_plot <- tax_plot + theme(axis.text.x = element_blank()) + coord_cartesian(ylim = c(1, y_top), clip = "off") +
  annotate("text", x = "Phylum", y = n_rows + 0.8, label = "Phylum", angle = title_angle, hjust = 0, size = 6 / .pt)
module_plot <- module_plot + theme(axis.text.x = element_blank()) + coord_cartesian(ylim = c(1, y_top), clip = "off") +
  annotate("text", x = c("Spring module", "Autumn module"), y = n_rows + 0.8, label = c("Spring module", "Autumn module"), angle = title_angle, hjust = 0, size = 6 / .pt)
abundance_plot <- abundance_plot + theme(axis.text.x = element_blank()) + coord_cartesian(ylim = c(1, y_top), clip = "off") +
  annotate("text", x = "Abundance", y = n_rows + 0.8, label = "Abundance", angle = title_angle, hjust = 0, size = 6 / .pt)
occupancy_plot <- occupancy_plot + theme(axis.text.x = element_blank()) + coord_cartesian(ylim = c(1, y_top), clip = "off") +
  annotate("text", x = "Prevalence", y = n_rows + 0.8, label = "Prevalence", angle = title_angle, hjust = 0, size = 6 / .pt)

plots <- core_plot + tax_plot + module_plot + abundance_plot + occupancy_plot +
  (titan_alk_plot + coord_cartesian(xlim = c(0, 250), ylim = c(1, y_top))) + (titan_DO_plot + coord_cartesian(xlim = c(7.5, 13), ylim = c(1, y_top))) + # same rows in every panel
  (titan_nitrate_plot + coord_cartesian(ylim = c(1, y_top))) + (titan_temp_plot + coord_cartesian(ylim = c(1, y_top))) +
  plot_layout(ncol = 9, widths = c(1, 1, 2, 1, 1, 10, 10, 10, 10), guides = "collect", axes = "collect")


legends <- (legend_tax + legend_module + legend_role + plot_layout(widths=c(1,1,1))) / 
  (legend_abundance + legend_occupancy + legend_inc + legend_dec + plot_layout(widths=c(1,1,1,1)))

combined <- plot_spacer() / plots / legends + # spacer gives room above the panels for the column titles
  plot_layout(guides = "collect", heights = c(0.5, 6, 1)) &
  theme(legend.position = "none",
        plot.margin = margin(t = -1, r = 1, b = -1, l = 1))

combined

ggsave(file.path(out_dir, "figures/titan.pdf"), combined, width = 180, height = 210, units = "mm", dpi = 600) 
ggsave(file.path(out_dir, "figures/titan.png"), combined, width = 180, height = 210, units = "mm", dpi = 600)


# indicators whose 90% bootstrap interval lies entirely below or above the community change point for their direction
cp_position <- sppmax_all %>%
  left_join(change_points, by = "variable") %>%
  mutate(direction = if_else(maxgrp == 1, "Sensitive", "Tolerant"),
         community_cp = if_else(maxgrp == 1, fsumz_dec, fsumz_inc), # sensitive vs fsumz-, tolerant vs fsumz+
         position = case_when(Q95 < community_cp ~ "below",
                              Q5 > community_cp ~ "above",
                              TRUE ~ "overlaps")) %>%
  left_join(tax_tab %>% rownames_to_column("taxon") %>% select(taxon, Phylum, Genus), by = "taxon")

cp_position %>% count(variable, direction, position) %>%
  pivot_wider(names_from = position, values_from = n, values_fill = 0) # counts

below_cp <- cp_position %>% filter(position == "below") %>% arrange(variable, direction, zenv.cp) # ASVs

write.csv(below_cp, file.path(out_dir, "below_cp.csv"))

