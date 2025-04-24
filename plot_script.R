# A script to make the necesary plots for the article
# Packages ----
library(tidyverse)
library(igraph)
library(brainGraph)
library(patchwork)
library(bayestestR)
library(netrankr)
library(ggrepel)
library(corrplot)
library(ggraph)
library(tidygraph)
library(patchwork)
library(gtable)
library(gt)
library(gtExtras)
library(cowplot)
theme_set(jtools::theme_nice())

# Color & shape assignment ----

safe_colorblind_palette <-
  c(
    "#88CCEE",
             "#CC6677",
             "#DDCC77",
             "#117733",
             "#332288",
             "#AA4499",
             "#44AA99",
             "#999933",
             "#882255",
             "#661100",
             "#6699CC",
             "#888888"
  )

dolphin_color <- safe_colorblind_palette[1]
elephant_color <- safe_colorblind_palette[3]
baboon_color <- safe_colorblind_palette[4]

species_colors <- c("Elephant" = elephant_color, 
                    "Dolphin" = dolphin_color, 
                    "Baboon" = baboon_color)

dolphin_shape <- 21
elephant_shape <- 22
baboon_shape <- 23

species_shapes <- c("Elephant" = 22, 
                    "Dolphin" = 21, 
                    "Baboon" = 24)
# add function to change legend 
shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

# Network metrics plot ----
# Baboon Data
bab_w1_age <- read_rds("Data/networks/bab_w1_age.rds")
bab_w2_age <- read_rds("Data/networks/bab_w2_age.rds")
bab_w3_a_age <- read_rds("Data/networks/bab_w3_a_age.rds")
bab_w3_b_age <- read_rds("Data/networks/bab_w3_b_age.rds")
bab_w4_a_age <- read_rds("Data/networks/bab_w4_a_age.rds")
bab_w4_b_age <- read_rds("Data/networks/bab_w4_b_age.rds")
bab_new_age <- read_rds("Data/networks/bab_new_age.rds")

# Dolphin Networks 
dolphin_w1 <- read_rds("Data/networks/dolphin_w1.rds")
dolphin_w2 <- read_rds("Data/networks/dolphin_w2.rds")
dolphin_w3 <- read_rds("Data/networks/dolphin_w3.rds")
dolphin_w4 <- read_rds("Data/networks/dolphin_w4.rds")
dolphin_w5 <- read_rds("Data/networks/dolphin_w5.rds")
dolphin_w6 <- read_rds("Data/networks/dolphin_w6.rds")

# Elephant Networks 
elph_w1_age <- read_rds("Data/networks/elph_w1_age.rds")
elph_w2_age <- read_rds("Data/networks/elph_w2_age.rds")
elph_w3_age <- read_rds("Data/networks/elph_w3_age.rds")

all_networks <- list(dolphin_w1, 
                     dolphin_w2, 
                     dolphin_w3, 
                     dolphin_w4, 
                     dolphin_w5, 
                     dolphin_w6, 
                     bab_w1_age, 
                     bab_w2_age, 
                     bab_new_age,
                     bab_w3_a_age, 
                     bab_w3_b_age, 
                     bab_w4_a_age, 
                     bab_w4_b_age, 
                     elph_w1_age, 
                     elph_w2_age, 
                     elph_w3_age )

# Network Features Function
multiple_network_features <- function(net) {
  # Edge density 
  dens <- edge_density(net)
  # Now community detection and modularity
  louvain_communities <- cluster_louvain(net)
  mod <- modularity(net, membership(louvain_communities))
  # Average transitivity 
  avg_trans <- transitivity(net, type = 'average')
  # Average weighted distance 
  avg_wt_dist <- mean_distance_wt(net, 
                                  weights = E(net)$weight)
  avg_dist <- mean_distance(net)
  mean_deg <- mean(degree(net))
  
  # All results the function will produce
  net_features <- c(num_nodes = length(V(net)), 
                    num_edges = length(E(net)), 
                    avg_trans = avg_trans, 
                    mod = mod, 
                    dens = dens, 
                    avg_wt_dist = avg_wt_dist, 
                    avg_dist = avg_dist)
  
  return(net_features)
}

# Get the features
set.seed(76)
network_features_df <- map_df(all_networks, 
                              multiple_network_features)

# Create guiding columns
network_features_df <- network_features_df %>% 
  mutate(wave = c(
                  1:6, 
                  1:7, 
                  1:3), 
         species = c(
           rep("Dolphin", 6), 
           rep("Baboon", 7), 
           rep("Elephant", 3))) %>% 
  select(species, wave, everything())

# Create table of all metrics across waves
net_fts_table <- network_features_df %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  gt() %>% 
  tab_header(
    title = md("Network Features across Species and Waves")
  )  %>%  
  gt_highlight_rows(
    rows = c(1,7, 14),  
    fill = "lightgrey",
  ) %>% 
  cols_align(
    align = c("center"),
    columns = everything()
  ) %>% 
  cols_move_to_start(
     columns = c(species, 
                 wave, 
                 avg_dist, 
                 avg_trans, 
                 avg_wt_dist, 
                 dens, 
                 mod, 
                 num_nodes, 
                 num_edges)
  ) %>% 
  cols_label(
    species = "Species", 
    wave = "Wave", 
    avg_dist = "Average Distance",
    avg_trans = "Average Transitivity", 
    avg_wt_dist = "Average Weighted Distance", 
    dens = "Density",
    mod = "Modularity",
    num_nodes = "Number of Nodes", 
    num_edges = "Number of Edges"
  ) 

# save table 
net_fts_table %>% 
  gtsave("Figures/network_feature_table.png", 
         expand = 10)

# Table collapsed across waves 
avg_fts_table <- network_features_df %>% 
  group_by(species) %>% 
  select(-wave) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  gt() %>% 
  tab_header(
    title = md("Average Network Features across Species")
  ) %>% 
  cols_move_to_start(
    columns = c(species, 
                avg_dist, 
                avg_trans, 
                avg_wt_dist, 
                dens, 
                mod, 
                num_nodes, 
                num_edges)
  ) %>% 
  cols_align(
    align = c("center"),
    columns = everything()
  ) %>% 
  cols_label(
    species = "Species", 
    avg_dist = "Average Distance",
    avg_wt_dist = "Average Weighted Distance", 
    avg_trans = "Average Transitivity", 
    mod = "Modularity",
    dens = "Density", 
    num_nodes = "Number of Nodes", 
    num_edges = "Number of Edges"
  ) 

# save averages table 
avg_fts_table %>% 
  gtsave("Figures/avg_feature_table.png", 
         expand = 10)

network_fts_plt <- network_features_df %>% 
  rename(
    "Average Distance" = avg_dist, 
    "Average Transitivity" = avg_trans, 
    "Average Weighted Distance" = avg_wt_dist, 
    "Density" = dens, 
    "Modularity" = mod, 
    "Number of Edges" = num_edges, 
    "Number of nodes" = num_nodes
  ) %>% 
  pivot_longer(cols = 3:ncol(network_features_df), 
               values_to = "value", 
               names_to = "feature") %>% 
  ggplot(aes(x = wave, y = value, pch = species)) +
  geom_jitter(size = 3, 
              color = "black", 
              fill = "gray80", 
              width = 0.2, 
              alpha = 0.7) + 
  scale_shape_manual(
    values = species_shapes
  ) +
  scale_x_continuous(breaks=1:7,
                     labels=c(1:7)) +
  facet_wrap(~feature, scales = "free_y", ncol = 3) + 
  theme(
    legend.position = "top"
  ) + 
  labs(x = "Wave", 
       y = "", 
       title = "Network Features", 
       pch = "Species") 

p <- network_fts_plt + 
  guides(fill = guide_legend(title.position = "top",
                             label.position = "bottom",
                             nrow = 1)) +
  theme(legend.direction = "horizontal")

grid::grid.draw(shift_legend(p))
# Save the plot 
ggsave("Figures/network_fts_plt.png", 
       height = 10, 
       width = 15, 
       dpi = 300)

# Network plots ----
plot(elph_w1_age, 
     vertex.label = "", 
     layout = layout.lgl, 
     vertex.size =  5, 
     edge.width = E(elph_w1_age)$weight*1.5, 
     vertex.color = "gray80", 
     vertex.shape = "square", 
     main = "Elephant - Wave 1")

# Add triangle shape 
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

plot(bab_w1_age,
     vertex.label = "", 
     vertex.color = rep("gray80",length(V(bab_w1_age))),
     vertex.size =  rep(8, length(V(bab_w1_age))),
     edge.width = E(bab_w1_age)$weight*5,
     vertex.shape = "triangle",
     main = "Baboon - Wave 1")

plot(dolphin_w3,
     vertex.label = "", 
     vertex.color = "gray80",
     vertex.size =  5,
     edge.width = E(dolphin_w3)$weight*3,
     main = "Dolphin - Wave 3")

# Efficiency plot ----


rems_df <- read_rds("Data/different_removals_df.rds")

rems_df <- rems_df %>% 
  mutate(
    species = case_when( 
      str_detect(net, "^elephant") ~ "Elephant", 
      str_detect(net, "^dolphin") ~ "Dolphin", 
      str_detect(net, "^bab") ~ "Baboon"
      ), 
    wave = case_when(
      str_detect(net, "w1$") ~ 1, 
      str_detect(net, "w2$") ~ 2, 
      str_detect(net, "w3$") ~ 3, 
      str_detect(net, "w4$") ~ 4, 
      str_detect(net, "w3a") ~ 3, 
      str_detect(net, "w3b") ~ 4,
      str_detect(net, "w5$") ~ 5, 
      str_detect(net, "w6$") ~ 6, 
      str_detect(net, "w4a") ~ 5,
      str_detect(net, "w4b") ~ 6, 
      str_detect(net, "new") ~ 7
    ), 
    wave = case_when(
      species == "Baboon" & 
        wave == 7 ~ 3, 
      species == "Baboon" & 
        wave == 3 ~ 4,
      species == "Baboon" & 
        wave == 4 ~ 5, 
      species == "Baboon" & 
        wave == 5 ~ 6, 
      species == "Baboon" & 
        wave == 6 ~ 7, 
      TRUE ~ wave
    )
  ) 

 
rdl <- rems_df %>% 
  rename(random = avg_eff) %>% 
  pivot_longer(cols = c("centrality", 
                        "degree", 
                        "ages", 
                        "random", 
                        "targeted"), 
               names_to = "removal_type", 
               values_to = "efficiency")  %>% 
  mutate(removal_type = case_when(
    removal_type == "ages" ~ "Age-based",
    removal_type == "centrality" ~ "Centrality",
    removal_type == "degree" ~ "Degree", 
    removal_type == "random" ~ "Random", 
    removal_type == "targeted" ~ "Targeted"
  ), 
  species = factor(species, levels = c("Baboon", "Dolphin", "Elephant")))

eff_plt <- rdl %>%
  filter(removal_type != "Targeted") %>%
  ggplot(aes(x = removals, y = efficiency, color = removal_type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                color = "#999933", 
                alpha = 0.5) +
  scale_color_manual(values = c("#332288",
                                         "#AA4499",
                                         "#44AA99",
                                         "#999933")) +
                                           geom_line() +
  facet_grid( species ~ wave, scales = "free_y") +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Global Efficiency After Removals",
    x = "Percent Nodes Removed",
    y = "Efficiency",
    color = "Removal Type"
  ) + 
  theme(legend.position = "bottom") +
  theme(axis.text = element_text(size = 7))

# Save the plot 
ggsave("Figures/eff_plt.png", 
       height = 10, 
       width = 15, 
       dpi = 300)

# Look at the elephants specifically

eleph_eff_plt <- rdl %>%
  filter(removal_type != "targeted", 
         species == "elephant") %>%
  ggplot(aes(x = removals, y = efficiency, color = removal_type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                color = "#999933") +
  scale_color_manual(values = c("#332288",
                                         "#AA4499",
                                         "#44AA99",
                                         "#999933")) +
                                           geom_line(alpha = 0.7) +
  facet_grid( species ~ wave) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Global efficiency after removals",
    subtitle = "Elephant networks",
    x = "Nodes removed",
    y = "Efficiency",
    color = "Removal"
  ) + 
  theme(legend.position = "top") +
  theme(axis.text = element_text(size = 7))

# Save the plot 
ggsave("Figures/eleph_eff_plt.png", 
       height = 10, 
       width = 15, 
       dpi = 300)

# Long trajectories plot (Only Young Explore) ----

ew1_rand <- readRDS("Results/ew1_rand_oy.rds")
ew2_rand <- readRDS("Results/ew2_rand_oy.rds")
ew3_rand <- readRDS("Results/ew3_rand_oy.rds")
dw1_rand <- readRDS("Results/dw1_rand_oy.rds")
dw2_rand <- readRDS("Results/dw2_rand_oy.rds")
dw3_rand <- readRDS("Results/dw3_rand_oy.rds")
dw4_rand <- readRDS("Results/dw4_rand_oy.rds")
dw5_rand <- readRDS("Results/dw5_rand_oy.rds")
dw6_rand <- readRDS("Results/dw6_rand_oy.rds")
bw1_rand <- read_rds("Results/bw1_rand_oy.rds")
bw2_rand <- read_rds("Results/bw2_rand_oy.rds")
bw3a_rand <- read_rds("Results/bw3a_rand_oy.rds")
bw3b_rand <- read_rds("Results/bw3b_rand_oy.rds")
bw4a_rand <- read_rds("Results/bw4a_rand_oy.rds")
bw4b_rand <- read_rds("Results/bw4b_rand_oy.rds")
bn_rand <- read_rds("Results/bn_rand_oy.rds")

ew1_rand <- ew1_rand %>% 
  mutate(species = "elephant", 
         wave = 1, 
         removal_type = "random")
ew2_rand <- ew2_rand %>% 
  mutate(species = "elephant", 
         wave = 2, 
         removal_type = "random")
ew3_rand <- ew3_rand %>% 
  mutate(species = "elephant", 
         wave = 3, 
         removal_type = "random")

dw1_rand <- dw1_rand %>% 
  mutate(species = "dolphin", 
         wave = 1, 
         removal_type = "random")
dw2_rand <- dw2_rand %>% 
  mutate(species = "dolphin", 
         wave = 2, 
         removal_type = "random")
dw3_rand <- dw3_rand %>% 
  mutate(species = "dolphin", 
         wave = 3, 
         removal_type = "random")
dw4_rand <- dw4_rand %>% 
  mutate(species = "dolphin", 
         wave = 4, 
         removal_type = "random")
dw5_rand <- dw5_rand %>% 
  mutate(species = "dolphin", 
         wave = 5, 
         removal_type = "random")
dw6_rand <- dw6_rand %>% 
  mutate(species = "dolphin", 
         wave = 6, 
         removal_type = "random")

bw1_rand <- bw1_rand %>% 
  mutate(species = "baboon", 
         wave = 1, 
         removal_type = "random")
bw2_rand <- bw2_rand %>% 
  mutate(species = "baboon", 
         wave = 2, 
         removal_type = "random")
bw3a_rand <- bw3a_rand %>% 
  mutate(species = "baboon", 
         wave = 3, 
         removal_type = "random")
bw3b_rand <- bw3b_rand %>% 
  mutate(species = "baboon", 
         wave = 4, 
         removal_type = "random")
bw4a_rand <- bw4a_rand %>% 
  mutate(species = "baboon", 
         wave = 5, 
         removal_type = "random")
bw4b_rand <- bw4b_rand %>% 
  mutate(species = "baboon", 
         wave = 6, 
         removal_type = "random")
bn_rand <- bn_rand %>% 
  mutate(species = "baboon", 
         wave = 7, 
         removal_type = "random")

all_data_r <- rbind(ew1_rand, 
                  ew2_rand, 
                  ew3_rand, 
                  dw1_rand, 
                  dw2_rand, 
                  dw3_rand, 
                  dw4_rand,
                  dw5_rand, 
                  dw6_rand, 
                  bw1_rand, 
                  bw2_rand, 
                  bw3a_rand,
                  bw3b_rand,
                  bw4a_rand,
                  bw4b_rand, 
                  bn_rand)

ew1 <- readRDS("Results/ew1_oy.rds")
ew2 <- readRDS("Results/ew2_oy.rds")
ew3 <- readRDS("Results/ew3_oy.rds")
dw1 <- readRDS("Results/dw1_oy.rds")
dw2 <- readRDS("Results/dw2_oy.rds")
dw3 <- readRDS("Results/dw3_oy.rds")
dw4 <- readRDS("Results/dw4_oy.rds")
dw5 <- readRDS("Results/dw5_oy.rds")
dw6 <- readRDS("Results/dw6_oy.rds")
bw1 <- read_rds("Results/bw1_oy.rds")
bw2 <- read_rds("Results/bw2_oy.rds")
bw3a <- read_rds("Results/bw3a_oy.rds")
bw3b <- read_rds("Results/bw3b_oy.rds")
bw4a <- read_rds("Results/bw4a_oy.rds")
bw4b <- read_rds("Results/bw4b_oy.rds")
bn <- read_rds("Results/bwn_oy.rds")

ew1 <- ew1 %>% 
  mutate(species = "elephant", 
         wave = 1, 
         removal_type = "age-based")
ew2 <- ew2 %>% 
  mutate(species = "elephant", 
         wave = 2, 
         removal_type = "age-based")
ew3 <- ew3 %>% 
  mutate(species = "elephant", 
         wave = 3, 
         removal_type = "age-based")

dw1 <- dw1 %>% 
  mutate(species = "dolphin", 
         wave = 1, 
         removal_type = "age-based")
dw2 <- dw2 %>% 
  mutate(species = "dolphin", 
         wave = 2, 
         removal_type = "age-based")
dw3 <- dw3 %>% 
  mutate(species = "dolphin", 
         wave = 3, 
         removal_type = "age-based")
dw4 <- dw4 %>% 
  mutate(species = "dolphin", 
         wave = 4, 
         removal_type = "age-based")
dw5 <- dw5 %>% 
  mutate(species = "dolphin", 
         wave = 5, 
         removal_type = "age-based")
dw6 <- dw6 %>% 
  mutate(species = "dolphin", 
         wave = 6, 
         removal_type = "age-based")

bw1 <- bw1 %>% 
  mutate(species = "baboon", 
         wave = 1, 
         removal_type = "age-based")
bw2 <- bw2 %>% 
  mutate(species = "baboon", 
         wave = 2, 
         removal_type = "age-based")
bw3a <- bw3a %>% 
  mutate(species = "baboon", 
         wave = 3, 
         removal_type = "age-based")
bw3b <- bw3b %>% 
  mutate(species = "baboon", 
         wave = 4, 
         removal_type = "age-based")
bw4a <- bw4a %>% 
  mutate(species = "baboon", 
         wave = 5, 
         removal_type = "age-based")
bw4b <- bw4b %>% 
  mutate(species = "baboon", 
         wave = 6, 
         removal_type = "age-based")
bn <- bn %>% 
  mutate(species = "baboon", 
         wave = 7, 
         removal_type = "age-based")

all_data_ab <- rbind(ew1, 
                     ew2, 
                     ew3, 
                     dw1, 
                     dw2, 
                     dw3, 
                     dw4,
                     dw5, 
                     dw6, 
                     bw1, 
                     bw2, 
                     bw3a,
                     bw3b,
                     bw4a,
                     bw4b, 
                     bn)

all_data <- rbind(all_data_ab, 
                  all_data_r)

rm(all_data_ab, 
   all_data_r)

srs_diff <- all_data %>% 
  select(species, wave, removal_type, p_exp, gm, remove, time, run) %>% 
  pivot_wider(names_from = remove, 
              values_from = gm)  %>% 
  mutate(diff = nr-r) %>% 
  group_by(species, wave, p_exp, time, removal_type) %>% 
  summarise(avg = median(diff), 
            avg_u = quantile(diff, .75), 
            avg_l = quantile(diff, .25))

baboon_diff <- srs_diff %>% 
  filter(species == "baboon") %>% 
  ggplot(aes(x = time, y = avg, 
             fill = removal_type, 
             col = removal_type)) +
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Baboons", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))

elephant_diff <- srs_diff %>% 
  filter(species == "elephant") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Elephants", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))


dolphin_diff <- srs_diff %>% 
  filter(species == "dolphin") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Dolphins", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))

diff_patches_oy <- baboon_diff + dolphin_diff + elephant_diff + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Trajectories by Removal Type", 
                  subtitle = "Only Young Explore") & theme(legend.position = 'bottom')

# Long trajectories plot (Only Old Explore) ----


ew1_rand <- readRDS("Results/ew1_rand_oo.rds")
ew2_rand <- readRDS("Results/ew2_rand_oo.rds")
ew3_rand <- readRDS("Results/ew3_rand_oo.rds")
dw1_rand <- readRDS("Results/dw1_rand_oo.rds")
dw2_rand <- readRDS("Results/dw2_rand_oo.rds")
dw3_rand <- readRDS("Results/dw3_rand_oo.rds")
dw4_rand <- readRDS("Results/dw4_rand_oo.rds")
dw5_rand <- readRDS("Results/dw5_rand_oo.rds")
dw6_rand <- readRDS("Results/dw6_rand_oo.rds")
bw1_rand <- read_rds("Results/bw1_rand_oo.rds")
bw2_rand <- read_rds("Results/bw2_rand_oo.rds")
bw3a_rand <- read_rds("Results/bw3a_rand_oo.rds")
bw3b_rand <- read_rds("Results/bw3b_rand_oo.rds")
bw4a_rand <- read_rds("Results/bw4a_rand_oo.rds")
bw4b_rand <- read_rds("Results/bw4b_rand_oo.rds")
bn_rand <- read_rds("Results/bn_rand_oo.rds")

ew1_rand <- ew1_rand %>% 
  mutate(species = "elephant", 
         wave = 1, 
         removal_type = "random")
ew2_rand <- ew2_rand %>% 
  mutate(species = "elephant", 
         wave = 2, 
         removal_type = "random")
ew3_rand <- ew3_rand %>% 
  mutate(species = "elephant", 
         wave = 3, 
         removal_type = "random")

dw1_rand <- dw1_rand %>% 
  mutate(species = "dolphin", 
         wave = 1, 
         removal_type = "random")
dw2_rand <- dw2_rand %>% 
  mutate(species = "dolphin", 
         wave = 2, 
         removal_type = "random")
dw3_rand <- dw3_rand %>% 
  mutate(species = "dolphin", 
         wave = 3, 
         removal_type = "random")
dw4_rand <- dw4_rand %>% 
  mutate(species = "dolphin", 
         wave = 4, 
         removal_type = "random")
dw5_rand <- dw5_rand %>% 
  mutate(species = "dolphin", 
         wave = 5, 
         removal_type = "random")
dw6_rand <- dw6_rand %>% 
  mutate(species = "dolphin", 
         wave = 6, 
         removal_type = "random")

bw1_rand <- bw1_rand %>% 
  mutate(species = "baboon", 
         wave = 1, 
         removal_type = "random")
bw2_rand <- bw2_rand %>% 
  mutate(species = "baboon", 
         wave = 2, 
         removal_type = "random")
bw3a_rand <- bw3a_rand %>% 
  mutate(species = "baboon", 
         wave = 3, 
         removal_type = "random")
bw3b_rand <- bw3b_rand %>% 
  mutate(species = "baboon", 
         wave = 4, 
         removal_type = "random")
bw4a_rand <- bw4a_rand %>% 
  mutate(species = "baboon", 
         wave = 5, 
         removal_type = "random")
bw4b_rand <- bw4b_rand %>% 
  mutate(species = "baboon", 
         wave = 6, 
         removal_type = "random")
bn_rand <- bn_rand %>% 
  mutate(species = "baboon", 
         wave = 7, 
         removal_type = "random")

all_data_r <- rbind(ew1_rand, 
                    ew2_rand, 
                    ew3_rand, 
                    dw1_rand, 
                    dw2_rand, 
                    dw3_rand, 
                    dw4_rand,
                    dw5_rand, 
                    dw6_rand, 
                    bw1_rand, 
                    bw2_rand, 
                    bw3a_rand,
                    bw3b_rand,
                    bw4a_rand,
                    bw4b_rand, 
                    bn_rand)

ew1 <- readRDS("Results/ew1_oo.rds")
ew2 <- readRDS("Results/ew2_oo.rds")
ew3 <- readRDS("Results/ew3_oo.rds")
dw1 <- readRDS("Results/dw1_oo.rds")
dw2 <- readRDS("Results/dw2_oo.rds")
dw3 <- readRDS("Results/dw3_oo.rds")
dw4 <- readRDS("Results/dw4_oo.rds")
dw5 <- readRDS("Results/dw5_oo.rds")
dw6 <- readRDS("Results/dw6_oo.rds")
bw1 <- read_rds("Results/bw1_oo.rds")
bw2 <- read_rds("Results/bw2_oo.rds")
bw3a <- read_rds("Results/bw3a_oo.rds")
bw3b <- read_rds("Results/bw3b_oo.rds")
bw4a <- read_rds("Results/bw4a_oo.rds")
bw4b <- read_rds("Results/bw4b_oo.rds")
bn <- read_rds("Results/bn_oo.rds")

ew1 <- ew1 %>% 
  mutate(species = "elephant", 
         wave = 1, 
         removal_type = "age-based")
ew2 <- ew2 %>% 
  mutate(species = "elephant", 
         wave = 2, 
         removal_type = "age-based")
ew3 <- ew3 %>% 
  mutate(species = "elephant", 
         wave = 3, 
         removal_type = "age-based")

dw1 <- dw1 %>% 
  mutate(species = "dolphin", 
         wave = 1, 
         removal_type = "age-based")
dw2 <- dw2 %>% 
  mutate(species = "dolphin", 
         wave = 2, 
         removal_type = "age-based")
dw3 <- dw3 %>% 
  mutate(species = "dolphin", 
         wave = 3, 
         removal_type = "age-based")
dw4 <- dw4 %>% 
  mutate(species = "dolphin", 
         wave = 4, 
         removal_type = "age-based")
dw5 <- dw5 %>% 
  mutate(species = "dolphin", 
         wave = 5, 
         removal_type = "age-based")
dw6 <- dw6 %>% 
  mutate(species = "dolphin", 
         wave = 6, 
         removal_type = "age-based")

bw1 <- bw1 %>% 
  mutate(species = "baboon", 
         wave = 1, 
         removal_type = "age-based")
bw2 <- bw2 %>% 
  mutate(species = "baboon", 
         wave = 2, 
         removal_type = "age-based")
bw3a <- bw3a %>% 
  mutate(species = "baboon", 
         wave = 3, 
         removal_type = "age-based")
bw3b <- bw3b %>% 
  mutate(species = "baboon", 
         wave = 4, 
         removal_type = "age-based")
bw4a <- bw4a %>% 
  mutate(species = "baboon", 
         wave = 5, 
         removal_type = "age-based")
bw4b <- bw4b %>% 
  mutate(species = "baboon", 
         wave = 6, 
         removal_type = "age-based")
bn <- bn %>% 
  mutate(species = "baboon", 
         wave = 7, 
         removal_type = "age-based")

all_data_ab <- rbind(ew1, 
                     ew2, 
                     ew3, 
                     dw1, 
                     dw2, 
                     dw3, 
                     dw4,
                     dw5, 
                     dw6, 
                     bw1, 
                     bw2, 
                     bw3a,
                     bw3b,
                     bw4a,
                     bw4b, 
                     bn)

all_data <- rbind(all_data_ab, 
                  all_data_r)

rm(all_data_ab, 
   all_data_r)

srs_diff <- all_data %>% 
  select(species, wave, removal_type, p_exp, gm, remove, time, run) %>% 
  pivot_wider(names_from = remove, 
              values_from = gm)  %>% 
  mutate(diff = nr-r) %>% 
  group_by(species, wave, p_exp, time, removal_type) %>% 
  summarise(avg = median(diff), 
            avg_u = quantile(diff, .75), 
            avg_l = quantile(diff, .25))

baboon_diff_oo <- srs_diff %>% 
  filter(species == "baboon") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Baboons", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))

elephant_diff_oo <- srs_diff %>% 
  filter(species == "elephant") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Elephants", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))


dolphin_diff_oo <- srs_diff %>% 
  filter(species == "dolphin") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Dolphins", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))

diff_patches_oo <- baboon_diff_oo + dolphin_diff_oo + elephant_diff_oo + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Trajectories by Removal Type", 
                  subtitle = "Only Old Explore") & theme(legend.position = 'bottom')

# Long trajectories plot (All Explore) ----

ew1_rand <- readRDS("Results/ew1_rand_all.rds")
ew2_rand <- readRDS("Results/ew2_rand_all.rds")
ew3_rand <- readRDS("Results/ew3_rand_all.rds")
dw1_rand <- readRDS("Results/dw1_rand_all.rds")
dw2_rand <- readRDS("Results/dw2_rand_all.rds")
dw3_rand <- readRDS("Results/dw3_rand_all.rds")
dw4_rand <- readRDS("Results/dw4_rand_all.rds")
dw5_rand <- readRDS("Results/dw5_rand_all.rds")
dw6_rand <- readRDS("Results/dw6_rand_all.rds")
bw1_rand <- read_rds("Results/bw1_rand_all.rds")
bw2_rand <- read_rds("Results/bw2_rand_all.rds")
bw3a_rand <- read_rds("Results/bw3a_rand_all.rds")
bw3b_rand <- read_rds("Results/bw3b_rand_all.rds")
bw4a_rand <- read_rds("Results/bw4a_rand_all.rds")
bw4b_rand <- read_rds("Results/bw4b_rand_all.rds")
bn_rand <- read_rds("Results/bn_rand_all.rds")

ew1_rand <- ew1_rand %>% 
  mutate(species = "elephant", 
         wave = 1, 
         removal_type = "random")
ew2_rand <- ew2_rand %>% 
  mutate(species = "elephant", 
         wave = 2, 
         removal_type = "random")
ew3_rand <- ew3_rand %>% 
  mutate(species = "elephant", 
         wave = 3, 
         removal_type = "random")

dw1_rand <- dw1_rand %>% 
  mutate(species = "dolphin", 
         wave = 1, 
         removal_type = "random")
dw2_rand <- dw2_rand %>% 
  mutate(species = "dolphin", 
         wave = 2, 
         removal_type = "random")
dw3_rand <- dw3_rand %>% 
  mutate(species = "dolphin", 
         wave = 3, 
         removal_type = "random")
dw4_rand <- dw4_rand %>% 
  mutate(species = "dolphin", 
         wave = 4, 
         removal_type = "random")
dw5_rand <- dw5_rand %>% 
  mutate(species = "dolphin", 
         wave = 5, 
         removal_type = "random")
dw6_rand <- dw6_rand %>% 
  mutate(species = "dolphin", 
         wave = 6, 
         removal_type = "random")

bw1_rand <- bw1_rand %>% 
  mutate(species = "baboon", 
         wave = 1, 
         removal_type = "random")
bw2_rand <- bw2_rand %>% 
  mutate(species = "baboon", 
         wave = 2, 
         removal_type = "random")
bw3a_rand <- bw3a_rand %>% 
  mutate(species = "baboon", 
         wave = 3, 
         removal_type = "random")
bw3b_rand <- bw3b_rand %>% 
  mutate(species = "baboon", 
         wave = 4, 
         removal_type = "random")
bw4a_rand <- bw4a_rand %>% 
  mutate(species = "baboon", 
         wave = 5, 
         removal_type = "random")
bw4b_rand <- bw4b_rand %>% 
  mutate(species = "baboon", 
         wave = 6, 
         removal_type = "random")
bn_rand <- bn_rand %>% 
  mutate(species = "baboon", 
         wave = 7, 
         removal_type = "random")

all_data_r <- rbind(ew1_rand, 
                    ew2_rand, 
                    ew3_rand, 
                    dw1_rand, 
                    dw2_rand, 
                    dw3_rand, 
                    dw4_rand,
                    dw5_rand, 
                    dw6_rand, 
                    bw1_rand, 
                    bw2_rand, 
                    bw3a_rand,
                    bw3b_rand,
                    bw4a_rand,
                    bw4b_rand, 
                    bn_rand)

ew1 <- readRDS("Results/ew1_all.rds")
ew2 <- readRDS("Results/ew2_all.rds")
ew3 <- readRDS("Results/ew3_all.rds")
dw1 <- readRDS("Results/dw1_all.rds")
dw2 <- readRDS("Results/dw2_all.rds")
dw3 <- readRDS("Results/dw3_all.rds")
dw4 <- readRDS("Results/dw4_all.rds")
dw5 <- readRDS("Results/dw5_all.rds")
dw6 <- readRDS("Results/dw6_all.rds")
bw1 <- read_rds("Results/bw1_all.rds")
bw2 <- read_rds("Results/bw2_all.rds")
bw3a <- read_rds("Results/bw3a_all.rds")
bw3b <- read_rds("Results/bw3b_all.rds")
bw4a <- read_rds("Results/bw4a_all.rds")
bw4b <- read_rds("Results/bw4b_all.rds")
bn <- read_rds("Results/bn_all.rds")

ew1 <- ew1 %>% 
  mutate(species = "elephant", 
         wave = 1, 
         removal_type = "age-based")
ew2 <- ew2 %>% 
  mutate(species = "elephant", 
         wave = 2, 
         removal_type = "age-based")
ew3 <- ew3 %>% 
  mutate(species = "elephant", 
         wave = 3, 
         removal_type = "age-based")

dw1 <- dw1 %>% 
  mutate(species = "dolphin", 
         wave = 1, 
         removal_type = "age-based")
dw2 <- dw2 %>% 
  mutate(species = "dolphin", 
         wave = 2, 
         removal_type = "age-based")
dw3 <- dw3 %>% 
  mutate(species = "dolphin", 
         wave = 3, 
         removal_type = "age-based")
dw4 <- dw4 %>% 
  mutate(species = "dolphin", 
         wave = 4, 
         removal_type = "age-based")
dw5 <- dw5 %>% 
  mutate(species = "dolphin", 
         wave = 5, 
         removal_type = "age-based")
dw6 <- dw6 %>% 
  mutate(species = "dolphin", 
         wave = 6, 
         removal_type = "age-based")

bw1 <- bw1 %>% 
  mutate(species = "baboon", 
         wave = 1, 
         removal_type = "age-based")
bw2 <- bw2 %>% 
  mutate(species = "baboon", 
         wave = 2, 
         removal_type = "age-based")
bw3a <- bw3a %>% 
  mutate(species = "baboon", 
         wave = 3, 
         removal_type = "age-based")
bw3b <- bw3b %>% 
  mutate(species = "baboon", 
         wave = 4, 
         removal_type = "age-based")
bw4a <- bw4a %>% 
  mutate(species = "baboon", 
         wave = 5, 
         removal_type = "age-based")
bw4b <- bw4b %>% 
  mutate(species = "baboon", 
         wave = 6, 
         removal_type = "age-based")
bn <- bn %>% 
  mutate(species = "baboon", 
         wave = 7, 
         removal_type = "age-based")

all_data_ab <- rbind(ew1, 
                     ew2, 
                     ew3, 
                     dw1, 
                     dw2, 
                     dw3, 
                     dw4,
                     dw5, 
                     dw6, 
                     bw1, 
                     bw2, 
                     bw3a,
                     bw3b,
                     bw4a,
                     bw4b, 
                     bn)

all_data <- rbind(all_data_ab, 
                  all_data_r)

rm(all_data_ab, 
   all_data_r)

srs_diff <- all_data %>% 
  select(species, wave, removal_type, p_exp, gm, remove, time, run) %>% 
  pivot_wider(names_from = remove, 
              values_from = gm)  %>% 
  mutate(diff = nr-r) %>% 
  group_by(species, wave, p_exp, time, removal_type) %>% 
  summarise(avg = median(diff), 
            avg_u = quantile(diff, .75), 
            avg_l = quantile(diff, .25))

baboon_diff_all <- srs_diff %>% 
  filter(species == "baboon") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_line() + 
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Baboons", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))

elephant_diff_all <- srs_diff %>% 
  filter(species == "elephant") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_line() + 
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Elephants", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))


dolphin_diff_all <- srs_diff %>% 
  filter(species == "dolphin") %>% 
  ggplot(aes(x = time, y = avg, fill = removal_type, 
             col = removal_type)) +
  geom_line() + 
  geom_ribbon(aes(ymin = avg_l, 
                  ymax = avg_u), 
              alpha = 0.2, 
              color = NA) +
  geom_hline(aes(yintercept = 0), col = "black", linetype = "dotted")+
  scale_color_manual(
    values = c("age-based" = "#332288", 
               "random" = "#AA4499"), 
    aesthetics = c("color", "fill")
  ) +
  facet_grid(wave~p_exp)  +
  labs(title = "Dolphins", 
       y = "Difference in Mean Fitness", 
       color = "Removal Type", 
       fill = "Removal Type") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(1, 250, 500),labels = c(1, 50, 100))

diff_patches_all <- baboon_diff_all + dolphin_diff_all + elephant_diff_all + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Trajectories by Removal Type", 
                  subtitle = "All Explore") & theme(legend.position = 'bottom')
