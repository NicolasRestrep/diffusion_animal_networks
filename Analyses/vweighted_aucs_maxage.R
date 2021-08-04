# Get the packages we need 
library(tidyverse)
library(igraph)

area_under_curve <- function (x, y, method = c("trapezoid", "step", "spline"), ...)
{
  if (length(x) != length(y)) {
    stop("length x must equal length y")
  }
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  switch(match.arg(arg = method, choices = c("trapezoid", "step",
                                             "spline")), trapezoid = sum((rowMeans(cbind(y[-length(y)],
                                                                                         y[-1]))) * (x[-1] - x[-length(x)])), step = sum(y[-length(y)] *
                                                                                                                                           (x[-1] - x[-length(x)])), spline = stats::integrate(stats::splinefun(x,
                                                                                                                                                                                                                y, method = "natural"), lower = min(x), upper = max(x))$value)
}

# Write function to return the graph from a wave
# This time I want age in there 
dolphin_edgelist <- function(w, t, node_file) {
  if (w == 1) {
    c <- "T2008"
    title <- "Wave 1"
    current_time <- as.Date("2008-12-31")
  } else {
    if(w ==2) {
      c <- "T2010"
      title <- "Wave 2"
      current_time <- as.Date("2010-12-31")
    } else {
      if (w==3) {
        c <- "T2012"
        title <- "Wave 3"
        current_time <- as.Date("2012-12-31")
      } else {
        if (w == 4) {
          c <- "T2014"
          title <- "Wave 4"
          current_time <- as.Date("2014-12-31")
        } else {
          if (w == 5) {
            c <- "T2016"
            title <- "Wave 5"
            current_time <- as.Date("2016-12-31")
          } else {
            c <- "T2018"
            title <- "Wave 6"
            current_time <- as.Date("2018-12-31")
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
  node_file$Birth.Date <- as.Date(node_file$Birth.Date)
  node_file$age <- current_time-node_file$Birth.Date
  
  net <- graph_from_data_frame(edgelist, directed = FALSE)
  df_ordered_nodes <- tibble(Dolphin.ID = V(net)$name)
  df_ordered_nodes <- df_ordered_nodes %>% 
    left_join(node_file, by = "Dolphin.ID")
  
  V(net)$age <- df_ordered_nodes$age
  return(net)
}
# Now let's read the data in 
# Edgelist for all networks 
dolphin_edge_lists <- read_csv("~/diffusion_animal_networks/Data/dolphin_edge_lists.csv")
# Attribute file
id_list <- read_csv("~/diffusion_animal_networks/Data/ID_list.csv")

# Wave 1 
dolphin_w1 <- dolphin_edgelist(w = 1, 
                               t = 0, 
                               node_file = id_list)
# Wave 2
dolphin_w2 <- dolphin_edgelist(w = 2, 
                               t = 0, 
                               node_file = id_list)
# Wave 3 
dolphin_w3 <- dolphin_edgelist(w = 3, 
                               t = 0, 
                               node_file = id_list)
# Wave 4 
dolphin_w4 <- dolphin_edgelist(w = 4, 
                               t = 0, 
                               node_file = id_list)
# Wave 5 
dolphin_w5 <- dolphin_edgelist(w = 5, 
                               t = 0, 
                               node_file = id_list)
# Wave 6
dolphin_w6 <- dolphin_edgelist(w = 6, 
                               t = 0, 
                               node_file = id_list)




# Now write function for vertical social learning
info_contagion_vertical <- function(net, rewire, e = 1, r_max, sim = 1){
  
  # Rewire network if random is set to TRUE
  if(rewire){
    net <- rewire(graph = net, with = keeping_degseq(loops = F, niter = 10^3))
  }
  
  # Get adjacency matrix from network
  adjm <- get.adjacency(net, 
                        sparse = F, 
                        attr = "weight")
  mat_inv <- adjm 
  edges <- which(mat_inv > 0)
  mat_inv[edges] <- 1.0001 - mat_inv[edges]
  inv_network <- graph_from_adjacency_matrix(mat_inv,
                                             mode = "undirected",
                                             weighted = TRUE)
  bcent <- betweenness(net, directed = FALSE, weights = E(inv_network)$weight)
  
  # Turn adjacency matrix into boolean (TRUE / FALSE) - if you dont want weights
  # adjm_bool <- adjm > 0
  
  # Set number of individuals based adjacency matrix
  N <- vcount(net)
  ages <- V(net)$age
  elder <- ages > quantile(ages, 0.9)
  seed <- which.max(bcent[elder])
  # Create a vector indicating possession of info and set one entry to TRUE
  info <- rep(FALSE, N)
  info[which.max(ages)] <- TRUE
  
  # Create a reporting variable
  proportion <- rep(0, r_max)
  
  # Rounds
  for(r in 1:r_max){
    # In random sequence go through all individuals without info
    for(i in sample(N)){
      # Select i's neighbourhood 
      ties <- adjm[i,] > 0
      nei <- ties & ages > ages[i]
      # If you dont want to include weights, quote above, unquote below
      #nei <- adjm_bool[i,]
      # Proceed if there is at least one neighbour
      if(sum(nei) > 0){
        # Simple contagion for e = 1 and complex contagion for e = 2
        if(runif(n = 1, min = 0, max = 1) <= (sum(adjm[i,which(nei & info)])/sum(adjm[i,]))^e){
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

# Function for efficiency 
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
multiple_diff_auc <- function(ev, g, g2, reps, turns) {
  # Place holder 
  values <- rep(NA, reps)
  
  for (i in 1:reps) {
    sim <- info_contagion_vertical(g, 
                                   rewire = F, 
                                   e = ev, 
                                   r_max = turns)
    auc_nr <- area_under_curve(x = sim$time, 
                               y = sim$proportion)
    
    # Now with removal
    sim_r <- info_contagion_vertical(g2, 
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
  rmv <- sample(1:vcount(net), 10)
  net_rems_random <- delete.vertices(net,rmv)
  medians <- rep(NA, num_iterations)
  for (i in 1:num_iterations) {
    if(removal_type == "targeted") {
      net_rems <- net_rems_targeted
    } else if(removal_type == "centrality") {
      net_rems <- net_rems_cent
    } else {
      net_rems <- net_rems_random
    }
    print(paste0("Working on: ",graph_id(net),"-",i))
    auc_diffs <- multiple_diff_auc(ev = 1, 
                                   g = net, 
                                   g2 = net_rems,
                                   reps = 100, 
                                   turns = 50)
    medians[[i]] <- median(auc_diffs)
  }
  return(medians)
}

# List of all the dolphin networks 
nets <- list(dolphin_w1, 
             dolphin_w2, 
             dolphin_w3, 
             dolphin_w4, 
             dolphin_w5, 
             dolphin_w6)

set.seed(76)
targeted_medians <- map(nets, 
                        distribution_of_medians, 
                        num_iterations = 25, 
                        removal_type = "targeted")

saveRDS(targeted_medians, 
        "~/diffusion_animal_networks/vweighted_aucs_maxage.rds")
