#### Packages ####
library(igraph)
library(tidyverse)

#### Build Networks ####
#### Dolphin Data ####
# Import the data
dolphin_edge_lists <- read_csv("~/Documents/diffusion_animal_networks/Data/dolphin_edge_lists_complete.csv")
# Import the ID list data
dolphin_id_list <- read_csv("~/Documents/diffusion_animal_networks/Data/ID_list_complete.csv")

# Write function to return the graph from a wave
dolphin_edgelist_females <- function(w, t, node_file) {
  # Conditional statements for the waves
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
  
  df_ordered_nodes <- tibble(Dolphin.ID = V(net)$name)
  df_ordered_nodes <- df_ordered_nodes %>%
    left_join(node_file, by = "Dolphin.ID")
  
  V(net)$sex <- df_ordered_nodes$Sex
  jfems <- delete.vertices(net,
                           V(net)$sex != "FEMALE")
  
  return(jfems)
}

# Wave 1
dolphin_w1 <- dolphin_edgelist_females(w = 1,
                                       t = 0,
                                       node_file = dolphin_id_list)
# Wave 2
dolphin_w2 <- dolphin_edgelist_females(w = 2,
                                       t = 0,
                                       node_file = dolphin_id_list)

# Wave 3
dolphin_w3 <- dolphin_edgelist_females(w = 3,
                                       t = 0,
                                       node_file = dolphin_id_list)

# Wave 4
dolphin_w4 <- dolphin_edgelist_females(w = 4,
                                       t = 0,
                                       node_file = dolphin_id_list)

# Wave 5
dolphin_w5 <- dolphin_edgelist_females(w = 5,
                                       t = 0,
                                       node_file = dolphin_id_list)

# Wave 6
dolphin_w6 <- dolphin_edgelist_females(w = 6,
                                       t = 0,
                                       node_file = dolphin_id_list)


#### Baboon Data ####
# Read in the networks
bab_w1 <- read_rds("~/Documents/diffusion_animal_networks/shelf/bab_w1.rds")
bab_w2 <- read_rds("~/Documents/diffusion_animal_networks/shelf/bab_w2.rds")
bab_w3_a <- read_rds("~/Documents/diffusion_animal_networks/shelf/bab_w3_a.rds")
bab_w3_b <- read_rds("~/Documents/diffusion_animal_networks/shelf/bab_w3_b.rds")
bab_w4_a <- read_rds("~/Documents/diffusion_animal_networks/shelf/bab_w4_a.rds")
bab_w4_b <- read_rds("~/Documents/diffusion_animal_networks/shelf/bab_w4_b.rds")
biography <- read_csv("~/Documents/diffusion_animal_networks/Data/baboon_data/Data_Tables/biography.csv")
missing_bab_data <-  read_csv("~/Documents/diffusion_animal_networks/Data/baboon_data/Data_Tables/missing_baboon_info.csv")
missing_bab_data <- missing_bab_data %>%
  select(sname, sex, birth, earliestbirth, mom, dad) %>%
  rename(statdate = earliestbirth)

# Get the ages & sexes for baboons
bab_data <- rbind(biography,
                  missing_bab_data)

# Now I am just getting the female
# Make a list of the networks
baboon_nets_list <- list(bab_w1,
                         bab_w2,
                         bab_w3_a,
                         bab_w3_b,
                         bab_w4_a,
                         bab_w4_b)

# Go through networks and add info about sex
for (i in 1:length(baboon_nets_list)) {
  df_ordered_nodes <- tibble(sname = V(baboon_nets_list[[i]])$name)
  df_ordered_nodes <- df_ordered_nodes %>%
    left_join(bab_data, by = "sname")
  
  V(baboon_nets_list[[i]])$sex <- df_ordered_nodes$sex
  baboon_nets_list[[i]] <- delete.vertices(baboon_nets_list[[i]],
                                           V(baboon_nets_list[[i]])$sex != "F")
  
}

# Take out the only-female nets
bab_w1_fem <- baboon_nets_list[[1]]
bab_w2_fem <- baboon_nets_list[[2]]
bab_w3_a_fem <- baboon_nets_list[[3]]
bab_w3_b_fem <- baboon_nets_list[[4]]
bab_w4_a_fem <- baboon_nets_list[[5]]
bab_w4_b_fem <- baboon_nets_list[[6]]

# Now let's get the ages
biography <- read_csv("~/Documents/diffusion_animal_networks/Data/baboon_data/Data_Tables/biography.csv")
missing_baboon_ages <- read_csv("~/Documents/diffusion_animal_networks/Data/baboon_data/Data_Tables/missing_baboon_ages.csv")
all_ages <- full_join(biography,
                      missing_baboon_ages %>% select(sname, birth))
all_ages$birth <- as.Date(all_ages$birth)

baboon_nets_list <- list(bab_w1_fem,
                         bab_w2_fem,
                         bab_w3_a_fem,
                         bab_w3_b_fem,
                         bab_w4_a_fem,
                         bab_w4_b_fem)

# Let's add ages to these vertices
for (i in 1:length(baboon_nets_list)) {
  if (i == 1) {
    current_time <- as.Date("1997-12-31")
  } else if (i == 2) {
    current_time <- as.Date("1998-12-31")
  } else if (i == 3 | i == 4) {
    current_time <- as.Date("2000-08-01")
  } else {
    current_time <- as.Date("2001-08-16")
  }
  
  wave_ages <- all_ages %>%
    mutate(age = current_time - birth)
  df_ordered_nodes <- tibble(sname = V(baboon_nets_list[[i]])$name)
  df_ordered_nodes <- df_ordered_nodes %>%
    left_join(wave_ages, by = "sname")
  
  V(baboon_nets_list[[i]])$age <- df_ordered_nodes$age
}

bab_w1_age <- baboon_nets_list[[1]]
bab_w2_age <- baboon_nets_list[[2]]
bab_w3_a_age <- baboon_nets_list[[3]]
bab_w3_b_age <- baboon_nets_list[[4]]
bab_w4_a_age <- baboon_nets_list[[5]]
bab_w4_b_age <- baboon_nets_list[[6]]



#### Elephant data ####

elephant_w1 <- read_rds("~/Documents/diffusion_animal_networks/shelf/elephant_w1.rds")
elephant_w2 <- read_rds("~/Documents/diffusion_animal_networks/shelf/elephant_w2.rds")
elephant_w3 <- read_rds("~/Documents/diffusion_animal_networks/shelf/elephant_w3.rds")

elephant_ages <- read_csv("~/Documents/diffusion_animal_networks/Data/SamburuAges_Nico.csv") %>%
  select(2:3)
elephant_ages2 <- read_csv("~/Documents/diffusion_animal_networks/Data/elephant_ages_1.csv")
elephant_ages3 <- read_csv("~/Documents/diffusion_animal_networks/Data/elephant_ages_2.csv")

colnames(elephant_ages) <- c('nodes', 'birth')
colnames(elephant_ages2) <- c('nodes', 'birth')
colnames(elephant_ages3) <- c('nodes', 'birth')

elephant_ages <- rbind(elephant_ages,
                       elephant_ages2,
                       elephant_ages3)

elephant_ages <- elephant_ages %>%
  mutate(age = 2021- as.numeric(birth))

dea <- distinct(elephant_ages,
                nodes,
                .keep_all = T)

# Make a list of the networks
elephant_nets_list <- list(elephant_w1,
                           elephant_w2,
                           elephant_w3)

for (i in 1:length(elephant_nets_list)) {
  
  edf_ordered_nodes <- tibble(nodes = V(elephant_nets_list[[i]])$name)
  edf_ordered_nodes <- edf_ordered_nodes %>%
    left_join(dea, by = "nodes")
  V(elephant_nets_list[[i]])$age <- edf_ordered_nodes$age
}

elph_w1_age <-elephant_nets_list[[1]]
elph_w2_age <- elephant_nets_list[[2]]
elph_w3_age <- elephant_nets_list[[3]]

#### Save the networks #### 

saveRDS(dolphin_w1, 
        "~/Documents/diffusion_animal_networks/Data/networks/dolphin_w1.rds")
saveRDS(dolphin_w2, 
        "~/Documents/diffusion_animal_networks/Data/networks/dolphin_w2.rds")
saveRDS(dolphin_w3, 
        "~/Documents/diffusion_animal_networks/Data/networks/dolphin_w3.rds")
saveRDS(dolphin_w4, 
        "~/Documents/diffusion_animal_networks/Data/networks/dolphin_w4.rds")
saveRDS(dolphin_w5, 
        "~/Documents/diffusion_animal_networks/Data/networks/dolphin_w5.rds")
saveRDS(dolphin_w6, 
        "~/Documents/diffusion_animal_networks/Data/networks/dolphin_w6.rds")

saveRDS(elph_w1_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/elph_w1_age.rds")
saveRDS(elph_w2_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/elph_w2_age.rds")
saveRDS(elph_w3_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/elph_w3_age.rds")


saveRDS(bab_w1_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/bab_w1_age.rds")
saveRDS(bab_w2_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/bab_w2_age.rds")
saveRDS(bab_w3_a_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/bab_w3_a_age.rds")
saveRDS(bab_w3_b_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/bab_w3_b_age.rds")
saveRDS(bab_w4_a_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/bab_w4_a_age.rds")
saveRDS(bab_w4_b_age, 
        "~/Documents/diffusion_animal_networks/Data/networks/bab_w4_b_age.rds")




