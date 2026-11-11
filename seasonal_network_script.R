# biofilm analysis - 18S
# seasonal co-occurrence networks


# load libraries
library(ggplot2)
library(dplyr)
library(microeco)
library(phyloseq)
library(tidyverse)
library(data.table)
library(patchwork)
library(file2meco)
library(igraph)
library(ggraph)
library(xml2)
library(Hmisc)
library(RColorBrewer)
library(glue)

set.seed(123)
options(scipen = 999)


# import data
ASV_tab <- data.frame(fread("data/processed_2021_23/ASV_counts_18S_processed.csv", sep = ",", header = TRUE, data.table = FALSE), check.names=F, row.names=1)
tax_tab <- read.table("data/processed_2021_23/taxonomy_clean.csv", header=T, row.names=1, check.names=F, sep=",")

meta_merged <- read.csv("data/metadata/meta_merged.csv", header = T, row.names = 1, check.names = F)

ASV_rel <- ASV_tab / rowSums(ASV_tab)
meta_merged <- meta_merged[rownames(ASV_rel), ]

out_dir <- "output/network"
dir.create(file.path(out_dir, "figures"), recursive = T, showWarnings = F)


# set thresholds
seasons <- c("Spring", "Autumn")
prevalence <- 0.1 
abund_thres <- 0.0001
cor_cut <- 0.3
p_cut <- 0.01
min_module_size <- 3
n_zoom <- 4 

env_vars <- c("Alkalinity_combined", "Conductivity", "pH", "Nitrate_N", "Orthophosphate", "Ammoniacal_N",
              "DOC", "SiO2", "DO", "Water_temp") 
new_names <- c("Alkalinity_combined" = "Alkalinity", "Ammoniacal_N" = "Ammoniacal-N",
               "Nitrate_N" = "Nitrate-N", "Water_temp" = "Temperature")


# network and modules
modules_list <- list()
stats_list <- list()
edges_list <- list()

for (s in seasons) {

  # subset samples to season and apply the prevalence filter within the season
  ASV_season <- ASV_tab[rownames(ASV_tab) %in% rownames(meta_merged)[meta_merged$season == s], ]
  ASV_sub <- ASV_season[, colSums(ASV_season > 0) >= prevalence * nrow(ASV_season)]

  # phyloseq to microeco object
  TAX <- tax_table(as.matrix(tax_tab %>% mutate(taxon = rownames(tax_tab))))
  ASV <- otu_table(ASV_sub, taxa_are_rows = F)
  samples <- sample_data(meta_merged[rownames(ASV_sub), ])
  ps <- phyloseq(ASV, TAX, samples)

  # convert to relative abundance
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

  # microeco object
  meco <- phyloseq2meco(ps_rel)

  # compute network modules
  t1 <- trans_network$new(dataset = meco,
                          cor_method = "spearman",
                          filter_thres = abund_thres,
                          taxa_level = "OTU",
                          nThreads = 10)

  t1$cal_network(COR_p_thres = p_cut, COR_cut = cor_cut) 
  
  t1$res_network <- delete_edges(t1$res_network, E(t1$res_network)[E(t1$res_network)$label == "-"]) # positive associations only - remove negative edges 
  t1$res_network <- delete_vertices(t1$res_network, V(t1$res_network)[degree(t1$res_network) == 0]) # then remove ASVs left with no edges
  
  t1$cal_module(method = "cluster_fast_greedy")
  t1$get_node_table(node_roles = TRUE)
  modules <- t1$res_node_table

  modules <- modules %>%
    group_by(module) %>%
    filter(n() > min_module_size) %>%  
    ungroup()

  # export network modules
  write.csv(modules, file.path(out_dir, glue("modules_{s}.csv")))
  t1$save_network(filepath = file.path(out_dir, glue("module_network_{s}.gexf")))

  # network properties (modules of more than X taxa)
  filtered_net <- t1[["subset_network"]](node = modules$name, return_igraph = TRUE)
  t1$res_network <- filtered_net
  t1$cal_network_attr()
  
  # ASV-ASV edge table (positive edges between ASVs in the retained modules)
  tax_info <- modules %>% select(name, module, Phylum, Class, Order, Family, Genus)
  edges <- igraph::as_data_frame(filtered_net, what = "edges") %>%
    select(ASV_1 = from, ASV_2 = to, rho = weight) %>%
    left_join(tax_info %>% rename_with(~ paste0(.x, "_1")), by = c("ASV_1" = "name_1")) %>%
    left_join(tax_info %>% rename_with(~ paste0(.x, "_2")), by = c("ASV_2" = "name_2")) %>%
    mutate(across(matches("^(Phylum|Class|Order|Family|Genus)_"), ~ sub("^[a-z]__", "", .x)),
           season = s) %>%
    arrange(desc(rho))
  edges_list[[s]] <- edges

  stats_list[[s]] <- data.frame(season = s,
                                samples = nrow(ASV_sub),
                                ASVs_tested = ncol(t1$data_abund),
                                t1$res_network_attr %>% t() %>% as.data.frame())
  modules_list[[s]] <- modules
}

network_stats <- bind_rows(stats_list)
write.csv(network_stats, file.path(out_dir, "network_stats.csv"), row.names = F)
network_stats

edges_all <- bind_rows(edges_list)
write.csv(edges_all, file.path(out_dir, "edges.csv"), row.names = F)


# set colours and shapes
phylum_colour_map <- read.table("output/community_structure/phylum_colour_map.csv", check.names=F, row.names=1, header=T, sep = ",")
phylum_colours <- c(setNames(phylum_colour_map$colour, phylum_colour_map$Phylum), Others = "grey60")
season_colours <- c(Spring = "forestgreen", Autumn = "steelblue") 
role_shapes <- c("Peripheral nodes" = 21, "Connectors" = 24, "Module hubs" = 22)

# module colours
pal_vec <- c("#9C9EDEFF", "#DE9ED6FF", "#CEDB9CFF", "#E7969CFF", "#E7CB94FF",
             "#6B6ECFFF", "#CE6DBDFF", "#B5CF6BFF", "#D6616BFF", "#E7BA52FF",
             "#5254A3FF", "#A55194FF", "#8CA252FF", "#AD494AFF", "#BD9E39FF",
             "#393B79FF", "#7B4173FF", "#637939FF", "#843C39FF", "#8C6D31FF")

# shared text sizes
theme_fig <- theme(plot.title = element_text(size = 7, colour = "black"),
                   strip.text = element_text(size = 6, colour = "black"),
                   legend.text = element_text(size = 5, colour = "black"),
                   legend.title = element_text(size = 6, colour = "black"),
                   legend.key.size = unit(0.5, "lines"))


# modules, nodes and edges per season
modules <- map_dfr(seasons, function(s) {
  modules_list[[s]] %>%
    select(ASV = name, module, degree, taxa_roles) %>%
    mutate(season = s)
}) %>%
  left_join(tax_tab %>% rownames_to_column("ASV") %>% select(ASV, Phylum), by = "ASV") %>%
  mutate(Phylum_top = if_else(Phylum %in% names(phylum_colours) & Phylum != "unassigned", Phylum, "Others"), # phyla without a colour (and unassigned) as Others
         module = paste0(substr(season, 1, 1), module), # SM1 = spring module 1, AM1 = autumn module 1
         role = factor(replace_na(taxa_roles, "Peripheral nodes"), levels = names(role_shapes)),
         key_node = role != "Peripheral nodes",
         season = factor(season, levels = seasons))

module_order <- map(setNames(seasons, seasons), function(s) {
  m <- unique(modules$module[modules$season == s])
  m[order(as.numeric(sub("^[SA]M", "", m)))]
})
n_mod <- lengths(module_order)

read_gexf <- function(file) {
  x <- read_xml(file) %>% xml_ns_strip()
  nodes <- xml_find_all(x, "//nodes/node")
  node_tab <- data.frame(id = xml_attr(nodes, "id"),
                         ASV = xml_attr(nodes, "label"),
                         Abundance = as.numeric(xml_attr(xml_find_first(nodes, ".//attvalue[@for='att3']"), "value"))) 
  edges <- xml_find_all(x, "//edges/edge")
  edge_tab <- data.frame(source = xml_attr(edges, "source"),
                         target = xml_attr(edges, "target"),
                         sign = xml_attr(edges, "label"),
                         rho = as.numeric(xml_attr(edges, "weight"))) %>%
    left_join(node_tab %>% select(id, ASV1 = ASV), by = c("source" = "id")) %>%
    left_join(node_tab %>% select(id, ASV2 = ASV), by = c("target" = "id")) %>%
    select(ASV1, ASV2, sign, rho)
  list(nodes = node_tab %>% select(-id), edges = edge_tab)
}

net <- map(setNames(seasons, seasons), ~ read_gexf(file.path(out_dir, paste0("module_network_", .x, ".gexf"))))

# module networks
graphs <- map(setNames(seasons, seasons), function(s) {
  nodes_s <- modules %>% filter(season == s) %>%
    select(name = ASV, module, degree, Phylum_top, role, key_node) %>%
    mutate(module = factor(module, levels = module_order[[s]]))
  edges_s <- net[[s]]$edges %>%
    filter(ASV1 %in% nodes_s$name, ASV2 %in% nodes_s$name) %>%
    select(from = ASV1, to = ASV2, sign)
  graph_from_data_frame(edges_s, vertices = nodes_s, directed = FALSE)
})

# one node size scale (degree) for the networks and the module close-ups
size_scale <- scale_size_continuous(range = c(0.6, 4), limits = c(1, max(modules$degree)), breaks = c(10, 20, 40), name = "Degree")
shape_scale <- scale_shape_manual(values = role_shapes, limits = c("Module hubs", "Connectors", "Peripheral nodes"),
                                  labels = c("Hubs", "Connectors", "Peripheral"), name = "Role",
                                  guide = guide_legend(override.aes = list(size = 2, fill = "grey70")))


# networks coloured by module
net_plot <- function(s) {
  set.seed(42)
  ggraph(graphs[[s]], layout = "fr") +
    geom_edge_link(aes(colour = sign), edge_width = 0.1, alpha = 0.25, show.legend = FALSE) +
    geom_node_point(aes(filter = !key_node, fill = module, size = degree, shape = role), colour = "grey20", stroke = 0.15) +
    geom_node_point(aes(filter = key_node, fill = module, size = degree, shape = role), colour = "black", stroke = 0.4,
                    show.legend = c(fill = FALSE, size = FALSE, shape = TRUE)) +
    scale_edge_colour_manual(values = c("+" = "grey50", "-" = "#D6616B")) +
    scale_fill_manual(values = setNames(pal_vec[seq_along(module_order[[s]])], module_order[[s]]), name = "Module",
                      guide = guide_legend(ncol = 1, override.aes = list(size = 2, shape = 21))) + 
    size_scale +
    shape_scale +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) + 
    theme_void() +
    theme_fig +
    theme(plot.title = element_text(size = 7, colour = "black", hjust = 0.5),
          legend.position = "right",
          legend.justification = "top", 
          legend.box = "vertical",
          legend.box.just = "left",
          legend.spacing.y = unit(3.5, "mm"), 
          legend.key.spacing.y = unit(0.5, "mm")) + 
    guides(size = "none", shape = "none") + 
    labs(title = s)
}

net_spring <- net_plot("Spring")
net_autumn <- net_plot("Autumn") +
  guides(fill = guide_legend(ncol = 1, order = 1, override.aes = list(size = 2, shape = 21)),
         size = guide_legend(ncol = 1, order = 2, override.aes = list(shape = 21, fill = "grey70")),
         shape = guide_legend(ncol = 1, order = 3, override.aes = list(size = 2, fill = "grey70")))


# 4 largest modules per season, coloured by phylum
largest <- modules %>% count(season, module, name = "n_ASVs") %>% group_by(season) %>%
  slice_max(n_ASVs, n = n_zoom, with_ties = F) %>% ungroup()

phyla_all <- c(setdiff(names(phylum_colours), "Others")[setdiff(names(phylum_colours), "Others") %in% modules$Phylum_top], "Others")

sub_plot <- function(s, mods, nrow = 2) { 
  g <- induced_subgraph(graphs[[s]], V(graphs[[s]])$name[V(graphs[[s]])$module %in% mods])
  ends <- ends(g, E(g))
  g <- delete_edges(g, which(V(g)$module[match(ends[, 1], V(g)$name)] != V(g)$module[match(ends[, 2], V(g)$name)]))
  V(g)$module_lab <- factor(paste0(V(g)$module, " (", table(V(g)$module)[as.character(V(g)$module)], " ASVs)"),
                            levels = paste0(mods, " (", table(V(g)$module)[mods], " ASVs)"))
  set.seed(42)
  ggraph(g, layout = "fr") +
    geom_edge_link(edge_width = 0.1, colour = "grey50", alpha = 0.3) +
    geom_node_point(aes(filter = !key_node, fill = Phylum_top, size = degree, shape = role), colour = "grey20", stroke = 0.15) +
    geom_node_point(aes(filter = key_node, fill = Phylum_top, size = degree, shape = role), colour = "black", stroke = 0.4,
                    show.legend = c(fill = FALSE, size = FALSE, shape = TRUE)) +
    facet_nodes(~ module_lab, scales = "free", nrow = nrow) +
    scale_fill_manual(values = phylum_colours, limits = phyla_all, name = "Phylum",
                      guide = guide_legend(nrow = 2, override.aes = list(size = 2, shape = 21))) +
    size_scale +
    shape_scale +
    theme_void() +
    theme_fig +
    theme(panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3),
          panel.spacing = unit(1.5, "mm"),
          strip.text = element_text(size = 6, colour = "black", margin = margin(0.5, 0, 0.5, 0)),
          legend.position = "none")
}

sub_spring <- sub_plot("Spring", largest$module[largest$season == "Spring"])
sub_autumn <- sub_plot("Autumn", largest$module[largest$season == "Autumn"])


# number of ASVs of each phylum per module
node_counts <- modules %>%
  count(season, module, Phylum_top, name = "n_ASVs") %>%
  mutate(module = factor(module, levels = c(module_order$Spring, module_order$Autumn)),
         Phylum_top = factor(Phylum_top, levels = phyla_all))

count_shared <- ggplot(node_counts %>% mutate(module = factor(module, levels = c(module_order$Spring, module_order$Autumn))),
                       aes(x = module, y = n_ASVs, fill = Phylum_top)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.05) +
  facet_grid(~ season, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(linewidth = 0.2),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme_fig +
  theme(legend.position = "none",
        strip.text = element_blank()) + 
  labs(y = "Number of ASVs") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = phylum_colours, limits = phyla_all)


# module relative abundance vs water chemistry 
module_cor <- map_dfr(seasons, function(s) {
  m <- modules %>% filter(season == s)
  samples_s <- rownames(meta_merged)[meta_merged$season == s]
  module_abund <- as.data.frame(t(rowsum(t(ASV_rel[samples_s, m$ASV]), group = m$module)))[, module_order[[s]]]
  rcors <- rcorr(as.matrix(cbind(module_abund, meta_merged[samples_s, env_vars])), type = "spearman")
  data.frame(season = s,
             module = rep(module_order[[s]], times = length(env_vars)),
             variable = rep(env_vars, each = n_mod[[s]]),
             r = as.vector(rcors$r[module_order[[s]], env_vars]),
             p = as.vector(rcors$P[module_order[[s]], env_vars])) %>%
    mutate(p_adj = p.adjust(p, method = "BH")) 
}) %>%
  mutate(significance = case_when(p_adj <= 0.001 ~ "***", p_adj <= 0.01 ~ "**", p_adj <= 0.05 ~ "*", TRUE ~ ""))

heat_plot <- function(s) {
  ggplot(module_cor %>% filter(season == s) %>%
           mutate(variable = factor(recode(variable, !!!new_names), levels = rev(recode(env_vars, !!!new_names))),
                  module = factor(module, levels = module_order[[s]])),
         aes(x = module, y = variable, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = significance), size = 1.4) +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 5, colour = "black"),
          axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1),
          axis.ticks = element_line(linewidth = 0.2)) +
    theme_fig +
    theme(legend.position = "right") +
    labs(fill = expression("Spearman" ~ rho)) + 
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu")), limits = c(-1, 1))
}

heat_spring <- heat_plot("Spring") + theme(legend.position = "none")
heat_autumn <- heat_plot("Autumn") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")


# legends
legend_of <- function(p) {
  g <- ggplotGrob(p)
  keep <- grepl("guide-box", g$layout$name) & !vapply(g$grobs, inherits, logical(1), "zeroGrob")
  g$grobs[[which(keep)[1]]]
}

legend_top <- function(p) {
  g <- legend_of(p)
  g$vp <- grid::viewport(x = 0, y = unit(1, "npc") - unit(3.5, "mm"), just = c("left", "top")) 
  g
}

legend_phylum <- legend_of(count_shared + theme(legend.position = "bottom", legend.title.position = "left") + labs(fill = "Phylum") +
                             guides(fill = guide_legend(nrow = 2)))
legend_heat <- legend_of(heat_plot("Autumn") + theme(legend.position = "bottom", legend.direction = "horizontal", legend.title.position = "left",
                                                      legend.key.width = unit(1.2, "lines")))


# combined figure
row_nets <- (net_spring + theme(legend.position = "none")) + wrap_elements(full = legend_top(net_spring)) +
  (net_autumn + theme(legend.position = "none")) + wrap_elements(full = legend_top(net_autumn)) +
  plot_layout(widths = c(1, 0.22, 1, 0.22))
row_subs <- sub_spring + plot_spacer() + sub_autumn + plot_layout(widths = c(1, 0.03, 1)) 
row_bottom <- count_shared + heat_spring + heat_autumn +
  plot_layout(widths = c(n_mod[["Spring"]] + n_mod[["Autumn"]], n_mod[["Spring"]], n_mod[["Autumn"]])) 
row_bottom_legends <- wrap_elements(full = legend_phylum) + wrap_elements(full = legend_heat) + plot_layout(widths = c(1, 1))

module_figure <- wrap_elements(full = row_nets) / wrap_elements(full = row_subs) / wrap_elements(full = row_bottom) / row_bottom_legends +
  plot_layout(heights = c(0.98, 1.08, 0.72, 0.12)) &
  theme(plot.margin = margin(0, 1, 0, 1))
module_figure <- module_figure + plot_annotation(theme = theme(plot.margin = margin(0, 0, 3, 0, "mm"))) 

ggsave(file.path(out_dir, "figures/network_module_figure.png"), module_figure, width = 180, height = 235, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "figures/network_module_figure.pdf"), module_figure, width = 180, height = 235, units = "mm", dpi = 600, bg = "white")


# tables
write.csv(modules %>% select(-Phylum_top), file.path(out_dir, "module_membership.csv"), row.names = F)
write.csv(node_counts, file.path(out_dir, "module_phylum_counts.csv"), row.names = F)
write.csv(module_cor, file.path(out_dir, "module_chemistry_correlations.csv"), row.names = F)


df <- edges_all %>%
  filter(Genus_1 != Genus_2, Genus_1 != "", Genus_2 != "") %>%
  arrange(desc(rho))

df <- edges_all %>%
  filter(Phylum_1 != Phylum_2, Phylum_1 != "", Phylum_2 != "") %>%
  arrange(desc(rho))




# genus composition of modules 1-4 in each season (share of module abundance)
n_genera <- 10 # top genera per module, the rest (and ASVs unassigned at genus level) go in Others
genus_modules <- c(rbind(paste0("SM", 1:4), paste0("AM", 1:4))) # SM1, AM1, SM2, AM2... so spring is the left column and autumn the right
tax_labels <- read.csv("data/processed_2021_23/taxonomy_labels.csv", header = T, check.names = F)

genus_comp <- map_dfr(seasons, function(s) {
  m <- modules %>% filter(season == s, module %in% genus_modules)
  samples_s <- rownames(meta_merged)[meta_merged$season == s]
  data.frame(ASV = m$ASV, module = m$module, Phylum_top = m$Phylum_top,
             abund = colSums(ASV_rel[samples_s, m$ASV]))
}) %>%
  left_join(tax_labels %>% select(ASV, is_genus, genus_label, genus_label_italic), by = "ASV") %>%
  mutate(genus = if_else(is_genus, genus_label, "Others")) %>% # unassigned at genus level go in Others
  group_by(module, genus) %>%
  summarise(abund = sum(abund), n_ASVs = n(), label = first(genus_label_italic), Phylum_top = names(which.max(table(Phylum_top))), .groups = "drop") %>%
  group_by(module) %>%
  mutate(share = 100 * abund / sum(abund),
         rank = rank(-if_else(genus == "Others", -1, abund), ties.method = "first"),
         genus = if_else(rank <= n_genera, genus, "Others"),
         Phylum_top = if_else(genus == "Others", "Others", Phylum_top)) %>%
  mutate(label = if_else(genus == "Others", "'Others'", label)) %>%
  group_by(module, genus, label, Phylum_top) %>%
  summarise(share = sum(share), n_ASVs = sum(n_ASVs), .groups = "drop") %>%
  group_by(module) %>%
  mutate(order = if_else(genus == "Others", -1, share),
         module_lab = factor(paste0(module, " (", sapply(module, function(x) sum(modules$module == x)), " ASVs)"),
                             levels = paste0(genus_modules, " (", sapply(genus_modules, function(x) sum(modules$module == x)), " ASVs)")),
         Phylum_top = factor(Phylum_top, levels = phyla_all)) %>%
  ungroup() %>%
  mutate(y = reorder(paste(module, label, sep = "___"), order)) # order bars within each module, Others at the bottom

genus_plot <- ggplot(genus_comp, aes(x = share, y = y, fill = Phylum_top)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.05, width = 0.75) +
  geom_text(aes(label = n_ASVs), hjust = -0.3, size = 1.6) +
  facet_wrap(~ module_lab, scales = "free_y", ncol = 2) + # same x scale in all panels
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        axis.ticks = element_line(linewidth = 0.2),
        strip.background = element_blank()) +
  theme_fig +
  theme(legend.position = "bottom",
        legend.title.position = "left") +
  labs(x = "% of module abundance", fill = "Phylum") +
  scale_y_discrete(labels = function(x) parse(text = sub("^.*___", "", x))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = phylum_colours, limits = phyla_all[phyla_all %in% genus_comp$Phylum_top]) + # only phyla that are plotted
  guides(fill = guide_legend(nrow = 2))
genus_plot

ggsave(file.path(out_dir, "figures/module_genus_composition.png"), genus_plot, width = 180, height = 180, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "figures/module_genus_composition.pdf"), genus_plot, width = 180, height = 180, units = "mm", dpi = 600, bg = "white")
write.csv(genus_comp %>% select(module, genus, Phylum = Phylum_top, share, n_ASVs), file.path(out_dir, "module_genus_composition.csv"), row.names = F)

