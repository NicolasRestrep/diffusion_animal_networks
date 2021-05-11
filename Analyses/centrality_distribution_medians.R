library(tidyverse)
library(igraph)
library(bayestestR)

####  Functions ####
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
info_contagion <- function(net, rewire, e = 1, r_max, sim = 1){
  
  # Rewire network if random is set to TRUE
  if(rewire){
    net <- rewire(graph = net, with = keeping_degseq(loops = F, niter = 10^3))
  }
  
  # Get adjacency matrix from network
  adjm <- get.adjacency(net, 
                        sparse = F, 
                        attr = "weight")
  
  # Turn adjacency matrix into boolean (TRUE / FALSE) - if you dont want weights
  # adjm_bool <- adjm > 0
  
  # Set number of individuals based adjacency matrix
  N <- vcount(net)
  
  # Create a vector indicating possession of info and set one entry to TRUE
  info <- rep(FALSE, N)
  info[sample(x = N, size = 1)] <- TRUE
  
  # Create a reporting variable
  proportion <- rep(0, r_max)
  
  # Rounds
  for(r in 1:r_max){
    # In random sequence go through all individuals without info
    for(i in sample(N)){
      # Select i's neighbourhood 
      nei <- adjm[i,] > 0
      # If you dont want to include weights, quote above, unquote below
      #nei <- adjm_bool[i,]
      # Proceed if there is at least one neighbour
      if(sum(nei) > 0){
        # Simple contagion for e = 1 and complex contagion for e = 2
        if(runif(n = 1, min = 0, max = 1) <= (sum(adjm[i,][info])/sum(nei))^e){
          info[i] <- TRUE
        }
      }
    }
    # Record proportion of the population with info
    proportion[r] <- sum(info) / N
    # Increment the round counter
    r <- r + 1
  }
  # Return a tibble with simulation results
  return(tibble(time = 1:r_max, 
                proportion = proportion, 
                time_to_max = which(proportion == max(proportion))[1],
                e = e, 
                network = ifelse(test = rewire, yes = "random", no = "model output"),
                sim = sim))
}
multiple_diff_auc <- function(ev, g, g2, reps, turns) {
  # Place holder 
  values <- rep(NA, reps)
  
  for (i in 1:reps) {
    sim <- info_contagion(g, 
                          rewire = F, 
                          e = ev, 
                          r_max = turns)
    auc_nr <- area_under_curve(x = sim$time, 
                               y = sim$proportion)
    
    # Now with removal
    sim_r <- info_contagion(g2, 
                            rewire = F, 
                            e = ev, 
                            r_max = turns)
    auc_r <- area_under_curve(x = sim_r$time, 
                              y = sim_r$proportion)
    
    diff <- auc_r/auc_nr
    
    values[i] <- diff
  }
  return(values)
}
distribution_of_medians <- function(net, 
                                    removal_type, 
                                    num_iterations) {
  removals <- decrease_efficiency(net)
  # Remove the top 10 in efficiency 
  top_nodes_targeted <- removals %>% 
    arrange(change_efficiency) %>% 
    slice(1:10) %>% 
    pull(node_name)
  # Pull top 10 centrality
  top_nodes_cent<- removals %>% 
    arrange(desc(betweenness)) %>% 
    slice(1:10) %>% 
    pull(node_name)
  
  # Remove the vertex with highest impact 
  net_rems_targeted <- delete.vertices(net, V(net)[top_nodes_targeted])
  net_rems_cent <- delete.vertices(net, V(net)[top_nodes_cent])
  
  medians <- rep(NA, num_iterations)
  for (i in 1:num_iterations) {
    if(removal_type == "targeted") {
      net_rems <- net_rems_targeted
    } else {
      net_rems <- net_rems_cent
    }
    print(paste0("Working on: ",graph_id(net),"-",i))
    auc_diffs <- multiple_diff_auc(ev = 1, 
                                   g = net, 
                                   g2 = net_rems,
                                   reps = 100, 
                                   turns = 350)
    medians[[i]] <- median(auc_diffs)
  }
  return(medians)
}
#### Aalysis ####
elephant_w1 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w1.rds")
elephant_w2 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w2.rds")
elephant_w3 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w3.rds")
elephant_w1_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w1_bb.rds")
elephant_w2_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w2_bb.rds")
elephant_w3_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w3_bb.rds")
elephant_w1_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w1_cut.rds")
elephant_w2_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w2_cut.rds")
elephant_w3_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/elephant_w3_cut.rds")
dolphin_w1 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w1.rds")
dolphin_w2 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w2.rds")
dolphin_w3 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w3.rds")
dolphin_w4 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w4.rds")
dolphin_w5 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w5.rds")
dolphin_w6 <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w6.rds")
dolphin_w1_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w1_bb.rds")
dolphin_w2_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w2_bb.rds")
dolphin_w3_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w3_bb.rds")
dolphin_w4_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w4_bb.rds")
dolphin_w5_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w5_bb.rds")
dolphin_w6_bb <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w6_bb.rds")
dolphin_w1_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w1_cut.rds")
dolphin_w2_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w2_cut.rds")
dolphin_w3_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w3_cut.rds")
dolphin_w4_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w4_cut.rds")
dolphin_w5_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w5_cut.rds")
dolphin_w6_cut <- read_rds("~/Documents/diffusion_animal_networks/Data/networks/dolphin_w6_cut.rds")

nets <- list(elephant_w1, 
             elephant_w2, 
             elephant_w3, 
             elephant_w1_bb, 
             elephant_w2_bb, 
             elephant_w3_bb, 
             elephant_w1_cut, 
             elephant_w2_cut, 
             elephant_w3_cut, 
             dolphin_w1, 
             dolphin_w2, 
             dolphin_w3, 
             dolphin_w4, 
             dolphin_w5, 
             dolphin_w6, 
             dolphin_w1_cut, 
             dolphin_w2_cut, 
             dolphin_w3_cut, 
             dolphin_w4_cut, 
             dolphin_w5_cut, 
             dolphin_w6_cut, 
             dolphin_w1_bb, 
             dolphin_w2_bb, 
             dolphin_w3_bb, 
             dolphin_w4_bb, 
             dolphin_w5_bb, 
             dolphin_w6_bb)

set.seed(76)
cent_medians <- map(nets, 
                        distribution_of_medians, 
                        num_iterations = 50, 
                        removal_type = "centrality")

saveRDS(cent_medians, 
        "cent_medians_aucs.rds")