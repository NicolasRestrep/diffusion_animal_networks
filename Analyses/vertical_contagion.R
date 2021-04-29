# Get the packages we need 
library(tidyverse)
library(igraph)
library(brainGraph)
library(patchwork)
library(bayestestR)
library(netrankr)
library(ggrepel)
library(corrplot)

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
# Edgelist
dolphin_edge_lists <- read_csv("Data/dolphin_edge_lists.csv")
# Attribute file
id_list <- read_csv("Data/ID_list.csv")

# Get the network 
dolphin_w1 <- dolphin_edgelist(w = 1, 
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
      ties <- adjm[i,] > 0
      nei <- ties & V(net)$age > V(net)$age[i]
      # If you dont want to include weights, quote above, unquote below
      #nei <- adjm_bool[i,]
      # Proceed if there is at least one neighbour
      if(sum(nei) > 0){
        # Simple contagion for e = 1 and complex contagion for e = 2
        if(runif(n = 1, min = 0, max = 1) <= (sum(adjm[i,which(nei & info)])/sum(nei))^e){
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

# Does it work? 

results <- info_contagion_vertical(net = dolphin_w1, 
                                   rewire = F, 
                                   e = 1, 
                                   r_max = 500, 
                                   sim = 1)

results %>% 
  ggplot(aes(x = time, 
             y = proportion)) + 
  geom_line()

# What about after some random removal 
rmv <- sample(1:vcount(dolphin_w1), 10)
ng <- delete.vertices(dolphin_w1, rmv)

results <- info_contagion_vertical(net = ng, 
                                   rewire = F, 
                                   e = 1, 
                                   r_max = 500, 
                                   sim = 1)

results %>% 
  ggplot(aes(x = time, 
             y = proportion)) + 
  geom_line()