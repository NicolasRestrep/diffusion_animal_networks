################################################
### Calculate Efficiency after Removal Types ###
################################################

#### Packages ####
library(tidyverse)
library(igraph)

#### Dolphin Data ####
# Dolphin Networks
dolphin_w1 <- read_rds("Data/networks/dolphin_w1.rds")
dolphin_w2 <- read_rds("Data/networks/dolphin_w2.rds")
dolphin_w3 <- read_rds("Data/networks/dolphin_w3.rds")
dolphin_w4 <- read_rds("Data/networks/dolphin_w4.rds")
dolphin_w5 <- read_rds("Data/networks/dolphin_w5.rds")
dolphin_w6 <- read_rds("Data/networks/dolphin_w6.rds")

#### Baboon Data ####
# Read in the networks
bab_w1_age <- read_rds("Data/networks/bab_w1_age.rds")
bab_w2_age <- read_rds("Data/networks/bab_w2_age.rds")
bab_w3_a_age <- read_rds("Data/networks/bab_w3_a_age.rds")
bab_w3_b_age <- read_rds("Data/networks/bab_w3_b_age.rds")
bab_w4_a_age <- read_rds("Data/networks/bab_w4_a_age.rds")
bab_w4_b_age <- read_rds("Data/networks/bab_w4_b_age.rds")

#### Elephant data ####

# Elephant Networks
elph_w1_age <- read_rds("Data/networks/elph_w1_age.rds")
elph_w2_age <- read_rds("Data/networks/elph_w2_age.rds")
elph_w3_age <- read_rds("Data/networks/elph_w3_age.rds")

#### Functions ####
inverse_efficiency <- function(g) {
  # Turn it into an adjacency matrix
  net_mat <- as_adj(g,
                    attr = 'weight',
                    sparse = F)
  # Get the inverse matrix
  mat_inv <- net_mat
  edges <- which(mat_inv > 0)
  mat_inv[edges] <- 1.000001 - mat_inv[edges]
  # Populate the diagonal with 0s
  diag(mat_inv) <- 0
  
  # Create the new graph
  net_inv <- graph_from_adjacency_matrix(mat_inv,
                                         weighted = T,
                                         mode = 'undirected')
  D <- distances(net_inv,
                 weights = E(net_inv)$weight)
  Nv <- nrow(D)
  edges <- which(D > 0)
  D[edges] <- D[edges] + 1
  diag(D) <- 0
  Dinv <- 1/D
  eff <- colSums(Dinv * is.finite(Dinv), na.rm = T)/(Nv - 1)
  geff <- sum(eff)/length(eff)
  
  return(geff)
}
decrease_efficiency <- function(g,
                                w = 1,
                                species = "elephant") {
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
  # Build an edgelist to find family ties
  edgelist <- get.data.frame(g)
  
  if ( species == "elephant") {
    
    # See which lines fullfil the requirements for kinship ties
    kinship_ties <- rep(NA, nrow(edgelist))
    
    for (i in 1:nrow(edgelist)) {
      pat_from <- edgelist$from[i]
      if (str_detect(edgelist$to[i], "\\.") != TRUE) {
        kinship_ties[i] <- 0
      } else {
        pat_to <- sub("\\..*", "", edgelist$to[i])
        if (pat_from==pat_to) {
          kinship_ties[i] <- 1
        } else {
          kinship_ties[i] <- 0
        }
      }
    }
    
    edgelist$kinship <- kinship_ties
    
    mothers <- edgelist %>%
      filter(kinship==1) %>%
      pull(from)
    
    children <- edgelist %>%
      filter(kinship==1) %>%
      pull(to)
    
    removal_df <- removal_df %>%
      mutate(mothers = if_else(node_name %in% mothers, 1, 0),
             children = if_else(node_name %in% children, 1, 0))
  } else if(species == "dolphin") {
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
      filter(!is.na(.[,3]) & .[,3] > 0) %>%
      rename(weight = c,
             from = ID1,
             to = ID2)
    related_df <- edgelist %>%
      group_by(from) %>%
      summarize(total_relatedness = sum(relatedness_coef>0, na.rm = T)) %>%
      select(from, total_relatedness) %>%
      rename(node_name = from)
    
    removal_df <- removal_df %>%
      left_join(related_df, by = "node_name")
    
    sex_df <- id_list_dolphin %>%
      select(1, w+1, 12) %>%
      rename(node_name = Dolphin.ID) %>%
      mutate(sex_binary = case_when(Sex == "MALE" ~ 0,
                                    Sex == "FEMALE" ~ 1)) %>%
      select(node_name, sex_binary)
    
    removal_df <- removal_df %>%
      left_join(sex_df, by = "node_name") %>%
      mutate(perc_related = total_relatedness/degree) %>% select(-total_relatedness)
    
  } else if (species == "baboon") {
    all_ages <- all_ages %>%
      rename(node_name = sname) %>%
      mutate(current_time = as.Date("2021-08-28")-as.Date(birth),
             current_time = as.numeric(current_time),
             sex_binary = if_else(sex == "F", 0, 1),
             sex_binary = as.numeric(sex_binary)) %>%
      select(node_name, sex_binary, current_time)
    
    removal_df <- removal_df %>%
      left_join(all_ages, by = "node_name")
  }
  
  return(removal_df)
}

efficiency_after_random_removal <- function(g,
                                            reps,
                                            num_removals,
                                            decrease_data,
                                            type_removal) {
  # placeholder
  efficiencies <- rep(NA, reps)
  nrems <- round(length(V(g))*(num_removals*1/100))
  if (type_removal == "random") {
    print(paste0(graph_id(g), " ", nrems, " "))
    for (i in 1:reps) {
      rmv <- sample(1:vcount(g), nrems)
      ng <- delete.vertices(g,rmv)
      eff <- inverse_efficiency(ng)
      efficiencies[i] <- eff
    }
    avg_eff <- mean(efficiencies)
    sd_eff <- sd(efficiencies)
    upper <- avg_eff + 1.96*sd_eff
    lower <- avg_eff - 1.96*sd_eff
    return(c(avg_eff = avg_eff,
             upper = upper,
             lower = lower,
             removals = num_removals))
  } else if (type_removal == "age") {
    print(paste0(graph_id(g), " ", nrems, " "))
    d <- data.frame(
      node = names(V(g)),
      age = V(g)$age
    )
    ntbr <- d %>%
      arrange(desc(age)) %>%
      slice(1:nrems) %>%
      pull(node)
    ng <- delete.vertices(g, ntbr)
    eff <- inverse_efficiency(ng)
    return(eff)
  } else if (type_removal == "centrality") {
    print(paste0(graph_id(g), " ", num_removals, " "))
    ntbr <- decrease_data %>%
      arrange(desc(betweenness)) %>%
      slice(1:nrems) %>%
      pull(node_name)
    ng <- delete.vertices(g, ntbr)
    eff <- inverse_efficiency(ng)
    return(eff)
  } else if (type_removal == "degree") {
    print(paste0(graph_id(g), " ", num_removals, " "))
    ntbr <- decrease_data %>%
      arrange(desc(degree)) %>%
      slice(1:nrems) %>%
      pull(node_name)
    ng <- delete.vertices(g, ntbr)
    eff <- inverse_efficiency(ng)
    return(eff)
  } else if (type_removal == "targeted") {
    print(paste0(graph_id(g), " ", num_removals, " "))
    ntbr <- decrease_data %>%
      arrange(change_efficiency) %>%
      slice(1:nrems) %>%
      pull(node_name)
    ng <- delete.vertices(g, ntbr)
    eff <- inverse_efficiency(ng)
  }
}

all_nets <- list(elph_w1_age,
                 elph_w2_age,
                 elph_w3_age,
                 dolphin_w1,
                 dolphin_w2,
                 dolphin_w3,
                 dolphin_w4,
                 dolphin_w5,
                 dolphin_w6,
                 bab_w1_age,
                 bab_w2_age,
                 bab_w3_a_age,
                 bab_w3_b_age,
                 bab_w4_a_age,
                 bab_w4_b_age)

nets_rep <- rep(all_nets,
                times = 20)



number_removals <- rep(seq(from = 1, to = 20),
                       each = 15)

random_rems_df <- map2_df(.x = nets_rep,
                          .y = number_removals,
                          efficiency_after_random_removal,
                          reps = 500,
                          type_removal = "random")
net_info <- rep(c("elephant_w1",
                  "elephant_w2",
                  "elephant_w3",
                  "dolphin_w1",
                  "dolphin_w2",
                  "dolphin_w3",
                  "dolphin_w4",
                  "dolphin_w5",
                  "dolphin_w6",
                  "bab_w1",
                  "bab_w2",
                  "bab_w3a",
                  "bab_w3b",
                  "bab_w4a",
                  "bab_w4b"), 20)
random_rems_df$net <- net_info

re1 <- decrease_efficiency(elph_w1_age)
re2 <- decrease_efficiency(elph_w2_age)
re3 <- decrease_efficiency(elph_w3_age)
rd1 <- decrease_efficiency(dolphin_w1)
rd2 <- decrease_efficiency(dolphin_w2)
rd3 <- decrease_efficiency(dolphin_w3)
rd4 <- decrease_efficiency(dolphin_w4)
rd5 <- decrease_efficiency(dolphin_w5)
rd6 <- decrease_efficiency(dolphin_w6)
rb1 <- decrease_efficiency(bab_w1_age)
rb2 <- decrease_efficiency(bab_w2_age)
rb3a <- decrease_efficiency(bab_w3_a_age)
rb3b <- decrease_efficiency(bab_w3_b_age)
rb4a <- decrease_efficiency(bab_w4_a_age)
rb4b <- decrease_efficiency(bab_w4_b_age)

all_rems <- list(re1,
                 re2,
                 re3,
                 rd1,
                 rd2,
                 rd3,
                 rd4,
                 rd5,
                 rd6,
                 rb1,
                 rb2,
                 rb3a,
                 rb3b,
                 rb4a,
                 rb4b)
vars_to_run <- tibble(
  g = nets_rep,
  reps = 1,
  num_removals = number_removals,
  decrease_data = rep(all_rems, 20),
  type_removal = "centrality"
)

centrality_rems <- pmap_dbl(vars_to_run,
                            efficiency_after_random_removal)

random_rems_df$centrality <- centrality_rems


vars_to_run <- tibble(
  g = nets_rep,
  reps = 1,
  num_removals = number_removals,
  decrease_data = rep(all_rems, 20),
  type_removal = "targeted"
)

targeted_rems <- pmap_dbl(vars_to_run,
                          efficiency_after_random_removal)

random_rems_df$targeted <- targeted_rems

vars_to_run <- tibble(
  g = nets_rep,
  reps = 1,
  num_removals = number_removals,
  decrease_data = rep(all_rems, 20),
  type_removal = "age"
)

age_rems <- pmap_dbl(vars_to_run,
                     efficiency_after_random_removal)

random_rems_df$ages <- age_rems

vars_to_run <- tibble(
  g = nets_rep,
  reps = 1,
  num_removals = number_removals,
  decrease_data = rep(all_rems, 20),
  type_removal = "degree"
)

deg_rems <- pmap_dbl(vars_to_run,
                     efficiency_after_random_removal)

random_rems_df$degree <- deg_rems


saveRDS(random_rems_df,
        "Data/different_removals_df.rds")