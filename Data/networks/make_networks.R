# Make networks ans save them 

# Get the packages we need 
library(tidyverse)
library(igraph)
library(brainGraph)
library(patchwork)
library(bayestestR)
library(netrankr)
library(ggrepel)
library(corrplot)
library(here)
library(readxl)

here::here()
#### Efficiency Functions ####
inverse_efficiency <- function(g) {
  # Turn it into an adjacency matrix 
  net_mat <- as_adj(g,
                    attr = 'weight', 
                    sparse = F)
  # Get the inverse matrix 
  mat_inv <- net_mat 
  edges <- which(mat_inv > 0)
  mat_inv[edges] <- 1.0001 - mat_inv[edges]
  # Populate the diagonal with 0s
  diag(mat_inv) <- 0
  
  # Create the new graph 
  net_inv <- graph_from_adjacency_matrix(mat_inv, 
                                         weighted = T, 
                                         mode = 'undirected')
  D <- distances(net_inv, 
                 weights = E(net_inv)$weight)
  D <- D + 1 
  diag(D) <- 0
  Nv <- nrow(D)
  Dinv <- 1/D
  eff <- colSums(Dinv * is.finite(Dinv), na.rm = T)/(Nv - 1)
  geff <- sum(eff)/length(eff)
  
  return(geff)
}

decrease_efficiency <- function(g) {
  # Get the original efficiency
  og_geff <- inverse_efficiency(g)
  # A matrix to store the data
  removal_df <- matrix(NA, ncol = 6, nrow = length(V(g)))
  # Inverse network 
  net_mat <- as_adj(g, 
                    attr = "weight",
                    sparse = FALSE)
  # Get the inverse matrix 
  mat_inv <- net_mat 
  edges <- which(mat_inv > 0)
  mat_inv[edges] <- 1.0001 - mat_inv[edges]
  inv_network <- graph_from_adjacency_matrix(mat_inv,
                                             mode = "undirected",
                                             weighted = TRUE)
  for (i in 1:length(V(g))) {
    vert <- V(g)[i]
    deg <- degree(g)[vert]
    ecent <- eigen_centrality(g, weights = E(g)$weight)$vector[vert]
    bcent <- betweenness(g, directed = FALSE, weights = E(inv_network)$weight)[vert]
    net_mat <- as_adj(g, attr = 'weight', sparse = F)
    sum_weigths <- sum(net_mat[vert,], na.rm = T)
    ng <- delete.vertices(g, vert)
    eff <- inverse_efficiency(ng)
    removal_df[i,] <- c(names(vert), 
                        eff-og_geff,
                        deg,
                        sum_weigths, 
                        ecent, 
                        bcent)
  }
  
  removal_df <- data.frame(removal_df)
  names(removal_df) <- c("node_name", "change_efficiency",
                         "degree", "sum_edge_weights", 
                         "eigen_centrality", 
                         "betweenness")
  removal_df <- removal_df %>% 
    mutate_at(.vars = c("change_efficiency",
                        "degree", "sum_edge_weights", 
                        "eigen_centrality", 
                        "betweenness"), 
              as.numeric)
  
  return(removal_df)
}

#### Elephant Data - No Removals ####

# Import the data 
elephant_data <- read_csv("/Users/nrestrepo/Documents/diffusion_animal_networks/Data/dist.matrix.t1.csv")
# Transform data into square matrix
elephant_matrix <- elephant_data %>% 
  select(2:ncol(elephant_data)) %>% 
  as.matrix() 
# Change the names so that nodes have the same ID as in the dataset
node_IDs <- names(elephant_data)[2:ncol(elephant_data)]
colnames(elephant_matrix) <- node_IDs
rownames(elephant_matrix) <- node_IDs
# Build inverse matrix
inv_elephant_matrix <- matrix(1, 97, 97) - elephant_matrix
# Populate the diagonal with 0s
diag(inv_elephant_matrix) <- 0
# Now create the network
elephant_graph_inv <- 
  graph_from_adjacency_matrix(inv_elephant_matrix, mode = "undirected", weighted = T)

# Network for wave 2 elephant data 
# Import the data 
elephant_data_w2 <- read_csv("/Users/nrestrepo/Documents/diffusion_animal_networks/Data/dist.matrix.t2.csv")
# Transform data into square matrix
elephant_matrix_w2 <- elephant_data_w2 %>% 
  select(2:ncol(elephant_data_w2)) %>% 
  as.matrix() 
# Change the names so that nodes have the same ID as in the dataset
node_IDs <- names(elephant_data_w2)[2:ncol(elephant_data_w2)]
colnames(elephant_matrix_w2) <- node_IDs
rownames(elephant_matrix_w2) <- node_IDs
# Build inverse matrix
inv_elephant_matrix_w2 <- matrix(1, 130, 130) - elephant_matrix_w2
# Populate the diagonal with 0s
diag(inv_elephant_matrix_w2) <- 0
# Replace with NAs with 0s 
inv_elephant_matrix_w2 <- replace_na(inv_elephant_matrix_w2, 0)
# Now create the network
elephant_graph_w2<- 
  graph_from_adjacency_matrix(inv_elephant_matrix_w2, mode = "undirected", weighted = T)

# Network elephant data wave 3 
# Import the data 
elephant_data_w3 <- read_csv("/Users/nrestrepo/Documents/diffusion_animal_networks/Data/dist.matrix.t3.csv")
# Transform data into square matrix
elephant_matrix_w3 <- elephant_data_w3 %>% 
  select(2:ncol(elephant_data_w3)) %>% 
  as.matrix() 
# Change the names so that nodes have the same ID as in the dataset
node_IDs <- names(elephant_data_w3)[2:ncol(elephant_data_w3)]
colnames(elephant_matrix_w3) <- node_IDs
rownames(elephant_matrix_w3) <- node_IDs
# Build inverse matrix
inv_elephant_matrix_w3 <- matrix(1, 120, 120) - elephant_matrix_w3
# Populate the diagonal with 0s
diag(inv_elephant_matrix_w3) <- 0
# Replace with NAs with 0s 
inv_elephant_matrix_w3 <- replace_na(inv_elephant_matrix_w3, 0)
# Now create the network
elephant_graph_w3<- 
  graph_from_adjacency_matrix(inv_elephant_matrix_w3, mode = "undirected", weighted = T)



#### Dolphin Data - No Removals ####

# Network dolphin data wave 1
# Import the data
dolphin_edge_lists <- read_csv("/Users/nrestrepo/Documents/diffusion_animal_networks/Data/dolphin_edge_lists.csv")

# Write function to return the graph from a wave
# Function to plot the networks 
dolphin_edgelist <- function(w, t) {
  # Conditional statements for the waves 
  if (w == 1) {
    c <- "T2008"
    title <- "Wave 1"
  } else {
    if(w ==2) {
      c <- "T2010"
      title <- "Wave 2"
    } else {
      if (w==3) {
        c <- "T2012"
        title <- "Wave 3"
      } else {
        if (w == 4) {
          c <- "T2014"
          title <- "Wave 4"
        } else {
          if (w == 5) {
            c <- "T2016"
            title <- "Wave 5"
          } else {
            c <- "T2018"
            title <- "Wave 6"
          }
        }
      }
    }
  }
  # Take the wave 
  edgelist <- dolphin_edge_lists %>% 
    select(1,2, c, 9) %>% 
    filter(!is.na(.[,3]) & .[,3] > t) %>% 
    rename(weight = c, 
           from = ID1, 
           to = ID2)
  
  net <- graph_from_data_frame(edgelist, directed = FALSE)
  return(net)
}
# Import the ID list data 
# Contains removal information 
id_list <- read_csv("/Users/nrestrepo/Documents/diffusion_animal_networks/Data/ID_list.csv")
dolphin_w1 <- dolphin_edgelist(w = 1, 
                               t = 0)

# Dolphin wave 2
dolphin_w2 <- dolphin_edgelist(w = 2, 
                               t = 0)

# Dolphin wave 3
dolphin_w3 <- dolphin_edgelist(w = 3, 
                               t = 0)

# Dolphin wave 4 
dolphin_w4 <- dolphin_edgelist(w = 4, 
                               t = 0)

# Dolphin wave 5 
dolphin_w5 <- dolphin_edgelist(w = 5, 
                               t = 0)

# Dolphin wave 6 
dolphin_w6 <- dolphin_edgelist(w = 6, 
                               t = 0)

#### Baboon Data - No Removals ####
# Importing the data
bab_dat <- read_excel("~/Documents/diffusion_animal_networks/Data/baboon_data/Data_Tables/grooming.xls")

# We'll create the matrix 
# 1) sub-setting data per year (there is an intentional gap in 1990 because of an atypical group fission event)
bab_dat1<- bab_dat%>%
  mutate(date = as.Date(date), format = "%Y-%m-%d")%>%
  filter(date <= "1997-12-31 ")%>%
  dplyr::select(actor, actee)%>%
  group_by(actor, actee)%>%
  count()%>%
  dplyr::rename(weight = n)

bab_dat2<- bab_dat%>%
  mutate(date = as.Date(date), format = "%Y-%m-%d")%>%
  filter(date > "1997-12-31" & date <= "1998-12-31")%>%
  dplyr::select(actor, actee)%>%
  group_by(actor, actee)%>%
  count()%>%
  dplyr::rename(weight = n)

bab_dat3<- bab_dat%>%
  mutate(date = as.Date(date), format = "%Y-%m-%d")%>%
  filter(date > "1999-07-31" & date <= "2000-08-01")%>%
  dplyr::select(actor, actee)%>%
  group_by(actor, actee)%>%
  count()%>%
  dplyr::rename(weight = n)

bab_dat4<- bab_dat%>%
  mutate(date = as.Date(date), format = "%Y-%m-%d")%>%
  filter(date > "2000-07-31" & date <= "2001-08-01")%>%
  dplyr::select(actor, actee)%>%
  group_by(actor, actee)%>%
  count()%>%
  dplyr::rename(weight = n)

# all nodes per dataset
babs1<- unique(c(bab_dat1$actor, bab_dat1$actee))
babs2<- unique(c(bab_dat2$actor, bab_dat2$actee))
babs3<- unique(c(bab_dat3$actor, bab_dat3$actee))
babs4<- unique(c(bab_dat4$actor, bab_dat4$actee))

# directed graphs
bab1_graph<- graph_from_data_frame(d = bab_dat1, vertices = babs1, direct = T) 
bab2_graph<- graph_from_data_frame(d = bab_dat2, vertices = babs2, direct = T) 
bab3_graph<- graph_from_data_frame(d = bab_dat3, vertices = babs3, direct = T) 
bab4_graph<- graph_from_data_frame(d = bab_dat4, vertices = babs4, direct = T) 

# collapse into undirected graph
bab1_graph<- as.undirected(bab1_graph, mode = "collapse",
                           edge.attr.comb = igraph_opt("edge.attr.comb"))
bab2_graph<- as.undirected(bab2_graph, mode = "collapse",
                           edge.attr.comb = igraph_opt("edge.attr.comb"))
bab3_graph<- as.undirected(bab3_graph, mode = "collapse",
                           edge.attr.comb = igraph_opt("edge.attr.comb"))
bab4_graph<- as.undirected(bab4_graph, mode = "collapse",
                           edge.attr.comb = igraph_opt("edge.attr.comb"))
# Standardize the networks 
bab1_mat<- as_adjacency_matrix(bab1_graph, type = "both", sparse = F, attr = "weight")
bab2_mat<- as_adjacency_matrix(bab2_graph, type = "both", sparse = F, attr = "weight")
bab3_mat<- as_adjacency_matrix(bab3_graph, type = "both", sparse = F, attr = "weight")
bab4_mat<- as_adjacency_matrix(bab4_graph, type = "both", sparse = F, attr = "weight")

bab1_mat_z <- bab1_mat/max(bab1_mat)
bab1_graph <- graph_from_adjacency_matrix(bab1_mat_z, 
                                          mode = "undirected", 
                                          weighted = T)

bab2_mat_z <- bab2_mat/max(bab2_mat)
bab2_graph <- graph_from_adjacency_matrix(bab2_mat_z, 
                                          mode = "undirected", 
                                          weighted = T)

bab3_mat_z <- bab3_mat/max(bab3_mat)
bab3_graph <- graph_from_adjacency_matrix(bab3_mat_z, 
                                          mode = "undirected", 
                                          weighted = T)

bab4_mat_z <- bab4_mat/max(bab4_mat)
bab4_graph <- graph_from_adjacency_matrix(bab4_mat_z, 
                                          mode = "undirected", 
                                          weighted = T)

b3cs <- components(bab3_graph)$membership

bab3_a <- delete.vertices(bab3_graph, 
                          names(which(b3cs == 1)))

bab3_b <- delete.vertices(bab3_graph, 
                          names(which(b3cs == 2)))

b4cs <- components(bab4_graph)$membership

bab4_a <- delete.vertices(bab4_graph, 
                          names(which(b4cs == 1)))

bab4_b <- delete.vertices(bab4_graph, 
                          names(which(b4cs == 2)))

#### Functions - Backboning ####
backbone_graph <- function(alpha, simple, g) {
  
  # All edges and their weight 
  e <- cbind(igraph::as_data_frame(g)[, 1:2 ], 
             weight = E(g)$weight)
  
  # Add graph strength and in degree of each source to graph dataframe
  # in
  # Sum all edge weigths for each node 
  w_in <- graph.strength(g,
                         mode = "in")
  w_in <- data.frame(to = names(w_in), 
                     w_in, stringsAsFactors = FALSE)
  # Get the degree for each node 
  k_in <- degree(g, mode = "in")
  k_in <- data.frame(to = names(k_in), 
                     k_in, 
                     stringsAsFactors = FALSE)
  
  # Join them both 
  e_in <- e %>%
    left_join(w_in, by = "to") %>%
    left_join(k_in, by = "to") %>%
    mutate(alpha_in = (1-(weight/w_in))^(k_in-1))
  
  # Same for out degree 
  w_out <- graph.strength(g, mode = "out")
  w_out <- data.frame(from = names(w_out), w_out, stringsAsFactors = FALSE)
  k_out <- degree(g, mode = "out")
  k_out <- data.frame(from = names(k_out), k_out, stringsAsFactors = FALSE)
  
  e_out <- e %>%
    left_join(w_out, by = "from") %>%
    left_join(k_out, by = "from") %>%
    mutate(alpha_out = (1-(weight/w_out))^(k_out-1))
  
  # Join everything and add alpha to the network
  #full
  e_full <- left_join(e_in, 
                      e_out, 
                      by = c("from", "to", "weight"))
  
  e_full <- e_full %>%
    mutate(alpha = ifelse(alpha_in < alpha_out, alpha_in, alpha_out)) %>%
    select(from, to, alpha)
  
  E(g)$alpha <- e_full$alpha
  
  # Finally prune the network 
  graph_pruned <- delete.edges(g, 
                               which(E(g)$alpha >= alpha))
  graph_pruned_deleted <- delete.vertices(graph_pruned, 
                                          which(degree(graph_pruned) == 0))
  
  # Count number of vertices
  num_e <- gsize(graph_pruned)
  num_v <- gorder(graph_pruned_deleted)
  
  res <- c(num_e, num_v)
  
  # return statement 
  if (simple == T) {
    return(res)
  } else {
    return(list(res, graph_pruned))
  }
}
backbone_rough <- function(g, quant) {
  copy_mat <- as_adjacency_matrix(g, 
                                  attr = 'weight')
  copy_mat[copy_mat < quantile(E(g)$weight, quant)] <- 0  
  net <- graph_from_adjacency_matrix(copy_mat, mode = "undirected", weighted = T)
  return(net)
} 
#### Backbone Networks ####
backbone_ew1_25 <- backbone_graph(alpha = 0.4, 
                                  simple = F, 
                                  g = elephant_graph_inv)
enet_w1 <- backbone_ew1_25[[2]]
enet_w1_cut <- backbone_rough(elephant_graph_inv,
                              quant = 0.75)
backbone_ew2_25 <- backbone_graph(alpha = 0.4, 
                                  simple = F, 
                                  g = elephant_graph_w2)
enet_w2 <- backbone_ew2_25[[2]]
enet_w2_cut <- backbone_rough(elephant_graph_w2,
                              quant = 0.75)
backbone_ew3_25 <- backbone_graph(alpha = 0.36, 
                                  simple = F, 
                                  g = elephant_graph_w3)
enet_w3 <- backbone_ew3_25[[2]]
enet_w3_cut <- backbone_rough(elephant_graph_w3,
                              quant = 0.75)
backbone_dw1_25 <- backbone_graph(alpha = 0.25, 
                                  simple = F, 
                                  g = dolphin_w1)
dnet_w1 <- backbone_dw1_25[[2]]
dnet_w1_cut <- backbone_rough(dolphin_w1,
                              quant = 0.75)
backbone_dw2_25 <- backbone_graph(alpha = 0.23, 
                                  simple = F, 
                                  g = dolphin_w2)
dnet_w2 <- backbone_dw2_25[[2]]
dnet_w2_cut <- backbone_rough(dolphin_w2,
                              quant = 0.75)
backbone_dw3_25 <- backbone_graph(alpha = 0.24, 
                                  simple = F, 
                                  g = dolphin_w3)
dnet_w3 <- backbone_dw3_25[[2]]
dnet_w3_cut <- backbone_rough(dolphin_w3,
                              quant = 0.75)
backbone_dw4_25 <- backbone_graph(alpha = 0.26, 
                                  simple = F, 
                                  g = dolphin_w4)
dnet_w4 <- backbone_dw4_25[[2]]
dnet_w4_cut <- backbone_rough(dolphin_w4,
                              quant = 0.75)
backbone_dw5_25 <- backbone_graph(alpha = 0.25, 
                                  simple = F, 
                                  g = dolphin_w5)
dnet_w5 <- backbone_dw5_25[[2]]
dnet_w5_cut <- backbone_rough(dolphin_w5,
                              quant = 0.75)
backbone_dw6_25 <- backbone_graph(alpha = 0.25, 
                                  simple = F, 
                                  g = dolphin_w6)
dnet_w6 <- backbone_dw6_25[[2]]
dnet_w6_cut <- backbone_rough(dolphin_w6,
                              quant = 0.75)

backbone_bw1_25 <- backbone_graph(alpha = 0.2, 
                                  simple = F, 
                                  g = bab1_graph)
bnet_w1 <- backbone_bw1_25[[2]]
bnet_w1_cut <- backbone_rough(bab1_graph,
                              quant = 0.75)

backbone_bw2_25 <- backbone_graph(alpha = 0.32, 
                                  simple = F, 
                                  g = bab2_graph)
bnet_w2 <- backbone_bw2_25[[2]]
bnet_w2_cut <- backbone_rough(bab2_graph,
                              quant = 0.75)

backbone_bw3a_25 <- backbone_graph(alpha = 0.23, 
                                  simple = F, 
                                  g = bab3_a)
bnet_w3a <- backbone_bw3a_25[[2]]
bnet_w3a_cut <- backbone_rough(bab3_a,
                              quant = 0.75)

backbone_bw3b_25 <- backbone_graph(alpha = 0.23, 
                                   simple = F, 
                                   g = bab3_b)
bnet_w3b <- backbone_bw3b_25[[2]]
bnet_w3b_cut <- backbone_rough(bab3_b,
                               quant = 0.75)

backbone_bw4a_25 <- backbone_graph(alpha = 0.25, 
                                   simple = F, 
                                   g = bab4_a)
bnet_w4a <- backbone_bw4a_25[[2]]
bnet_w4a_cut <- backbone_rough(bab4_a,
                               quant = 0.75)

backbone_bw4b_25 <- backbone_graph(alpha = 0.3, 
                                   simple = F, 
                                   g = bab4_b)
bnet_w4b <- backbone_bw4b_25[[2]]
bnet_w4b_cut <- backbone_rough(bab4_b,
                               quant = 0.75)

#### Save networks 
saveRDS(elephant_graph_inv, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w1.rds")
saveRDS(elephant_graph_w2, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w2.rds")
saveRDS(elephant_graph_w3, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w3.rds")
saveRDS(enet_w1, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w1_bb.rds")
saveRDS(enet_w2, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w2_bb.rds")
saveRDS(enet_w3, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w3_bb.rds")
saveRDS(enet_w1_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w1_cut.rds")
saveRDS(enet_w2_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w2_cut.rds")
saveRDS(enet_w3_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/elephant_w3_cut.rds")
saveRDS(dolphin_w1, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w1.rds")
saveRDS(dolphin_w2, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w2.rds")
saveRDS(dolphin_w3, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w3.rds")
saveRDS(dolphin_w4, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w4.rds")
saveRDS(dolphin_w5, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w5.rds")
saveRDS(dolphin_w6, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w6.rds")
saveRDS(dnet_w1, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w1_bb.rds")
saveRDS(dnet_w2, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w2_bb.rds")
saveRDS(dnet_w3, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w3_bb.rds")
saveRDS(dnet_w4, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w4_bb.rds")
saveRDS(dnet_w5, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w5_bb.rds")
saveRDS(dnet_w6, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w6_bb.rds")
saveRDS(dnet_w1_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w1_cut.rds")
saveRDS(dnet_w2_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w2_cut.rds")
saveRDS(dnet_w3_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w3_cut.rds")

saveRDS(dnet_w4_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w4_cut.rds")
saveRDS(dnet_w5_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w5_cut.rds")
saveRDS(dnet_w6_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/dolphin_w6_cut.rds")
saveRDS(bab1_graph , 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w1.rds")
saveRDS(bab2_graph , 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w2.rds")
saveRDS(bab3_a , 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w3_a.rds")
saveRDS(bab3_b , 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w3_b.rds")
saveRDS(bab4_a , 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w4_a.rds")
saveRDS(bab4_b , 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w4_b.rds")
saveRDS(bnet_w1 , 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w1_bb.rds")
saveRDS(bnet_w2, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w2_bb.rds")
saveRDS(bnet_w2, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w2_bb.rds")
saveRDS(bnet_w3a, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w3a_bb.rds")
saveRDS(bnet_w3b, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w3b_bb.rds")
saveRDS(bnet_w4a, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w4a_bb.rds")
saveRDS(bnet_w4b, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w4b_bb.rds")
saveRDS(bnet_w1_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w1_cut.rds")
saveRDS(bnet_w2_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w2_cut.rds")
saveRDS(bnet_w3a_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w3a_cut.rds")
saveRDS(bnet_w3b_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w3b_cut.rds")
saveRDS(bnet_w4a_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w4a_cut.rds")
saveRDS(bnet_w4b_cut, 
        "/Users/nrestrepo/Documents/diffusion_animal_networks/Data/networks/bab_w4b_cut.rds")


















