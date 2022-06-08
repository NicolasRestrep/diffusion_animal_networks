library(igraph)
library(tidyverse)

set.seed(NULL)

net <- read_rds("Data/networks/trial_dolphin.rds")
landscape <- read_rds("Data/trial_landscape.rds")
r_max <- 100
sims <- 5
p_exp <- 0.1
age_thr <- 25



nrems <- round(length(V(net))*0.1)

d <- data.frame(
  node = names(V(net)),
  age = V(net)$age
)

ntbr <- d %>%
  arrange(desc(age)) %>%
  slice(1:nrems) %>%
  pull(node)

#net <- delete.vertices(net, ntbr)


sample_robust <-  function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}
create_initial_positions <- function(net, 
                                     r_max, 
                                     landscape) {
  # Let's do the starting positions 
  # These stay constant across sims 
  # Create an empty dataset for the actors positions
  agents <- matrix(NA,
                   nrow = length(V(net)),
                   ncol = r_max)
  # Grab positions equal to the number of agents 
  
  sampled_postions <- landscape[sort(sample(1:2 ^ 12, length(V(net)), replace = F)), 
                                1:3]
  
  sorted_ages <- tibble(id = V(net), 
                        age = V(net)$age) %>% 
    arrange(desc(age))
  
  sorted_postions <- sampled_postions %>% arrange(desc(fitness))
  
  ages_positions <- cbind(sorted_ages, 
                          sorted_postions) %>% 
    arrange(id) # back to the original order
  
  return(ages_positions)
  
}

for (j in 1:sims) {
  
  print(paste0("working on simulation ", j))
  
  # Get adjacency matrix of the network 
  adjm <- get.adjacency(net,
                        sparse = F,
                        attr = "weight")
  # An empty list to keep record of the results 
  dfs <- list()
  
  # Create an empty dataset for the actors positions
  agents <- matrix(NA,
                   nrow = length(V(net)),
                   ncol = r_max)
  
  # Grab positions equal to the number of agents 
  
  sampled_postions <- landscape[sort(sample(1:2 ^ 12, length(V(net)), replace = F)), 
                                1:3]
  
  sorted_ages <- tibble(id = V(net), 
                        age = V(net)$age) %>% 
    arrange(desc(age))
  
  sorted_postions <- sampled_postions %>% arrange(desc(fitness))
  
  ages_positions <- cbind(sorted_ages, 
                          sorted_postions) %>% 
    arrange(id) # back to the original order
  
  print(ages_positions[1,])
  # Place agents in a position related to age 
  agents[, 1] <-
    ages_positions[,'position_name']
  
  agents_fitness <- matrix(NA,
                           nrow = length(V(net)),
                           ncol = r_max)
  
  for (z in 1:length(V(net))) {
    agents_fitness[z,1] <- landscape[landscape$position_name == agents[z,1], "fitness"]
  }
  
  # Record starting positions
  starting_positions <- agents[, 1]
  N <- vcount(net)
  ages <- V(net)$age
  ages_years <- ages/365
  
  learners <- which(ages_years<age_thr)
  nonlearners <- which(ages_years>=age_thr)
  
  
  maximum <- max(landscape[, "fitness"])
  
  med_fit <- c()
  max_fit <- c()
  min_fit <- c()
  numb_positions <- c()
  learners_med <- c()
  nonlearners_med <- c()
  
  med_fit[1] <- median(agents_fitness[,1])
  max_fit[1] <- max(agents_fitness[,1])
  min_fit[1] <- min(agents_fitness[,1])
  numb_positions[1] <- n_distinct(starting_positions) 
  
  # Rounds
  for (r in 2:r_max) {
    
    # for adults, stay where you are
    agents[nonlearners, r] <- agents[nonlearners, r - 1]
    agents_fitness[nonlearners, r] <- agents_fitness[nonlearners, r-1] 
    # Now, for the learners 
    for (i in 1:length(learners)) {
      
      # Select i's neighbourhood
      ties <- adjm[learners[i],] > 0
      nei <- ties & ages > ages[learners[i]]
      
      # Proceed if there is at least one neighbour
      if (sum(nei) > 0) {
        # Does the agent explore?
        explores <- rbernoulli(1, p = p_exp)
        if (explores == TRUE) {
          what_slot <- sample(1:12, 1)
          num_position <- as.numeric(str_split(agents[learners[i], r-1], "")[[1]])
          num_position[what_slot] <- ifelse(num_position[what_slot] == 1,
                                            0,
                                            1)
          explore_position <- paste(num_position, collapse = "")
          
          
          agents[learners[i], r] <-
            ifelse(landscape[landscape$position_name == explore_position, "fitness"] >
                     landscape[landscape$position_name == agents[learners[i], r - 1], "fitness"],
                   explore_position,
                   agents[learners[i], r -
                            1])
          agents_fitness[learners[i],r] <- landscape[landscape$position_name == agents[learners[i],r], "fitness"]
          try(if(agents_fitness[learners[i],r] < agents_fitness[learners[i],r-1]) stop("someone switched when they shouldn't"))
        } else {
          # If vertical learning, copy best elder
          nei_positions <- agents[nei, r - 1]
          
          best_position <- landscape[landscape$position_name %in% nei_positions,"position_name"][which.max(landscape[landscape$position_name %in% nei_positions, "fitness"])]
          best_elder <- sample_robust(x = c(which(nei_positions==best_position)), size = 1)
          agents[learners[i], r] <-
            ifelse(
              landscape[landscape$position_name == nei_positions[best_elder], "fitness"] >
                landscape[landscape$position_name == agents[learners[i], r - 1], "fitness"],
              nei_positions[best_elder],
              agents[learners[i], r - 1])
          agents_fitness[learners[i],r] <- landscape[landscape$position_name == agents[learners[i],r], "fitness"]
          try(if(agents_fitness[learners[i],r] < agents_fitness[learners[i],r-1]) stop("someone switched when they shouldn't"))
        }
      } else {
        # If none, just stay where you are
        agents[learners[i], r] <- agents[learners[i], r - 1]
        agents_fitness[learners[i],r] <- landscape[landscape$position_name == agents[learners[i],r], "fitness"]
        try(if(agents_fitness[i,r] < agents_fitness[i,r-1]) stop("someone switched when they shouldn't"))
      }
    }
    # Record relative max
    positions <- agents[,r]
    med_fit[r] <- median(agents_fitness[,r])
    max_fit[r] <- max(agents_fitness[,r])
    min_fit[r] <- min(agents_fitness[,r])
    numb_positions[r] <- n_distinct(positions)
    learners_med[r] <- median(agents_fitness[learners,r])
    nonlearners_med[r] <- median(agents_fitness[nonlearners,r])
  }
  # Return a tibble with simulation results
  df <- tibble(
    time = 1:r_max,
    num_agents = vcount(net),
    med = med_fit,
    maxf = max_fit, 
    minf = min_fit,
    numb_positions = numb_positions, 
    learners_med = learners_med, 
    nonlearners_med = nonlearners_med, 
    p_exp = p_exp
  )
  dfs[[j]] <- df
}
df_complete <- map_df(dfs, cbind)
df_complete <- df_complete %>%
  mutate(run = rep(1:sims, each = r_max))
mean(agents_fitness[,100])
