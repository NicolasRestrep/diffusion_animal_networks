#### Packages ####
library(tidyverse)
library(igraph)
library(furrr)

#### Load landscape and net ####
net <- read_rds("~/Documents/animal_networks_simulations/Data/networks/bab_w3_b_age.rds")
landscape <- read_rds("~/Documents/animal_networks_simulations/Data/trial_landscape.rds")

#### Functions ####
sample_robust <-  function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}
geometric_mean <- function (x, na.rm = TRUE)
{
  if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = na.rm))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}

rugged_crawl_vertical <- function(r_max,
                                  net,
                                  landscape,
                                  sims,
                                  p_exp,
                                  age_thr,
                                  perc_remove = 0.1,
                                  removal, 
                                  explore_thr) {
  # An empty list to keep record of the results
  dfs <- list()
  
  for (j in 1:sims) {
    
    # Print a record of where we are
    print(paste0("working on simulation ", j))
    
    # Now we will create the initial positions
    # Grab positions equal to the number of agents
    
    sampled_postions <- landscape[sort(sample(1:2 ^ 12, length(V(net)), replace = F)),
                                  1:3]
    
    # Match ids with ages
    sorted_ages <- tibble(id = V(net),
                          age = V(net)$age) %>%
      arrange(desc(age))
    
    sorted_postions <- sampled_postions %>% arrange(desc(fitness))
    
    ages_positions <- cbind(sorted_ages,
                            sorted_postions) %>%
      arrange(id) # back to the original order
    
    print(ages_positions[1,])
    
    # Now create records for the agents
    
    if (removal == TRUE) {
      
      nrems <- round(length(V(net))*perc_remove)
      
      d <- data.frame(
        node = names(V(net)),
        age = V(net)$age
      )
      
      ntbr <- d %>%
        arrange(desc(age)) %>%
        slice(1:nrems) %>%
        pull(node)
      
      net_rems <- delete.vertices(net, ntbr)
      
      # Get adjacency matrix of the network
      adjm <- get.adjacency(net_rems,
                            sparse = F,
                            attr = "weight")
      
      # Create an empty dataset for the actors positions
      agents <- matrix(NA,
                       nrow = length(V(net_rems)),
                       ncol = r_max)
      
      ages_positions <- ages_positions %>%
        arrange(desc(age)) %>%
        slice(-c(1:nrems)) %>%
        arrange(id)
      
      
      # Place agents in a position related to age
      agents[, 1] <-
        ages_positions[,'position_name']
      
      agents_fitness <- matrix(NA,
                               nrow = length(V(net_rems)),
                               ncol = r_max)
      
      for (z in 1:length(V(net_rems))) {
        agents_fitness[z,1] <- landscape[landscape$position_name == agents[z,1], "fitness"]
      }
      
      # Record starting positions
      starting_positions <- agents[, 1]
      N <- vcount(net_rems)
      ages <- V(net_rems)$age
      
      learners <- which(ages<age_thr)
      nonlearners <- which(ages>=age_thr)
      
    } else {
      
      # Get adjacency matrix of the network
      adjm <- get.adjacency(net,
                            sparse = F,
                            attr = "weight")
      
      # Create an empty dataset for the actors positions
      agents <- matrix(NA,
                       nrow = length(V(net)),
                       ncol = r_max)
      
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
      
      learners <- which(ages<age_thr)
      nonlearners <- which(ages>=age_thr)
    }
    
    maximum <- max(landscape[, "fitness"])
    
    med_fit <- vector("double", r_max)
    mn_fit <- vector("double", r_max)
    gm_fit <- vector("double", r_max)
    numb_positions <- vector("double", r_max)
    learners_med <- vector("double", r_max)
    nonlearners_med <- vector("double", r_max)
    
    med_fit[1] <- median(agents_fitness[,1])
    mn_fit[1] <- mean(agents_fitness[,1])
    gm_fit[1] <- geometric_mean(agents_fitness[,1])
    numb_positions[1] <- n_distinct(starting_positions)
    
    # Rounds
    for (r in 2:r_max) {
      
      # Now, for the learners
      for (i in 1:length(learners)) {
        
        # Select i's neighbourhood
        ties <- adjm[learners[i],] > 0
        nei <- ties & ages > ages[learners[i]]
        
        # Does the agent explore?
        explores <- rbernoulli(1, p = p_exp)
        
        if (explores == TRUE & ages[learners[i]] < explore_thr) {
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
          if (sum(nei) > 0) {
            # If vertical learning, copy best elder
            nei_positions <- agents[nei, r - 1]
            
            best_position <- landscape[landscape$position_name %in% nei_positions,"position_name"][which.max(landscape[landscape$position_name %in% nei_positions, "fitness"])]
            best_elder <- sample_robust(x = c(which(nei_positions==best_position)), size = 1)
            socially_learns <- rbernoulli(1, p = adjm[i,names(nei[nei==T][best_elder])])
            agents[learners[i], r] <-
              ifelse(
                (landscape[landscape$position_name == nei_positions[best_elder], "fitness"] >
                   landscape[landscape$position_name == agents[learners[i], r - 1], "fitness"]) &
                  socially_learns == T,
                nei_positions[best_elder],
                agents[learners[i], r - 1])
            agents_fitness[learners[i],r] <- landscape[landscape$position_name == agents[learners[i],r], "fitness"]
            try(if(agents_fitness[learners[i],r] < agents_fitness[learners[i],r-1]) stop("someone switched when they shouldn't"))
          } else {
            # If none, just stay where you are
            agents[learners[i], r] <- agents[learners[i], r - 1]
            agents_fitness[learners[i],r] <- landscape[landscape$position_name == agents[learners[i],r], "fitness"]
            try(if(agents_fitness[i,r] < agents_fitness[i,r-1]) stop("someone switched when they shouldn't"))
          }
        }
      }
      # Record relative max
      positions <- agents[,r]
      med_fit[r] <- median(agents_fitness[,r])
      mn_fit[r] <- mean(agents_fitness[,r])
      gm_fit[r] <- geometric_mean(agents_fitness[,r])
      numb_positions[r] <- n_distinct(positions)
      learners_med[r] <- median(agents_fitness[learners,r])
      nonlearners_med[r] <- median(agents_fitness[nonlearners,r])
    }
    # Return a tibble with simulation results
    df <- tibble(
      time = 1:r_max,
      num_agents = vcount(net),
      med = med_fit,
      mn = mn_fit,
      gm = gm_fit,
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
  return(df_complete)
}

#### Parameter exploration ####

pars <- expand_grid(p_exp = seq(from = 0, to = 0.2, by = 0.05),
                    removal = c(TRUE, FALSE))
V(net)$age <- V(net)$age/365
plan(multisession, 
     workers = 10)
options <- furrr_options(seed =202312)
df <- future_pmap_dfr(pars,
              rugged_crawl_vertical,
              r_max = 500,
              net = net,
              landscape = landscape,
              sims = 100,
              age_thr = 1e3,
              perc_remove = 0.1, 
              explore_thr = 14, 
              .options = options)

#### Clean and Save ####

df <- df %>%
  mutate(remove = rep(rep(c("r", "nr"), each = 50000),times = 5),
         remove = as.factor(remove))

saveRDS(df,
        "~/Documents/animal_networks_simulations/Results/bw3b_oy.rds")
