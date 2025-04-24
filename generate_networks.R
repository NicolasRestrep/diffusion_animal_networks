#####################################
### Create networks and save them ###
#####################################

#### Packages ####
library(igraph)
library(tidyverse)
library(readxl)

#### Build Networks ####
#### Dolphin Networks ####
# Import the data
dolphin_edge_lists <- read_csv("Data/dolphin_edge_lists_complete.csv")
# Import the ID list data
dolphin_id_list <- read_csv("Data/ID_list_complete.csv")

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
# Importing the data
bab_dat <- read_excel("Data/baboon_data/Data_Tables/grooming.xls")

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

# Now the turn adjacency matrices 
bab_w1 <- graph_from_adjacency_matrix(bab1_mat_z, 
                                          mode = "undirected", 
                                          weighted = T)

bab2_mat_z <- bab2_mat/max(bab2_mat)
bab_w2 <- graph_from_adjacency_matrix(bab2_mat_z, 
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

# Nets with multiple components get split
b3cs <- components(bab3_graph)$membership

bab_w3_a <- delete.vertices(bab3_graph, 
                          names(which(b3cs == 1)))

bab_w3_b <- delete.vertices(bab3_graph, 
                          names(which(b3cs == 2)))

b4cs <- components(bab4_graph)$membership

bab_w4_a <- delete.vertices(bab4_graph, 
                          names(which(b4cs == 1)))

bab_w4_b <- delete.vertices(bab4_graph, 
                          names(which(b4cs == 2)))

# This is the period that is wonky
bab_dat_new <- bab_dat%>%
  mutate(date = as.Date(date), format = "%Y-%m-%d")%>%
  filter(date > "1998-12-31" & date <= "1999-07-31")%>%
  dplyr::select(actor, actee)%>%
  group_by(actor, actee)%>%
  count()%>%
  dplyr::rename(weight = n)

babs_new<- unique(c(bab_dat_new$actor, bab_dat_new$actee))

# Make directed graphs 
bab_new_graph<- graph_from_data_frame(d = bab_dat_new, vertices = babs_new, direct = T) 


# Make it undirected 
bab_new_graph<- as.undirected(bab_new_graph, mode = "collapse",
                              edge.attr.comb = igraph_opt("edge.attr.comb"))

bab_new_mat<- as_adjacency_matrix(bab_new_graph, type = "both", sparse = F, attr = "weight")

bab_new_mat_z <- bab_new_mat/max(bab_new_mat)

bab_new <- graph_from_adjacency_matrix(bab_new_mat_z, 
                                             mode = "undirected", 
                                             weighted = T)

# Read in data for ages 
biography <- read_csv("Data/baboon_data/Data_Tables/biography.csv")
missing_bab_data <-  read_csv("Data/baboon_data/Data_Tables/missing_baboon_info.csv")
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
                         bab_w4_b, 
                         bab_new)

# Go through networks and add info about sex
# And only keep females
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
bab_new_fem <- baboon_nets_list[[7]]

# Now let's get the ages
biography <- read_csv("Data/baboon_data/Data_Tables/biography.csv")
missing_baboon_ages <- read_csv("Data/baboon_data/Data_Tables/missing_baboon_ages.csv")
all_ages <- full_join(biography,
                      missing_baboon_ages %>% select(sname, birth))
all_ages$birth <- as.Date(all_ages$birth)

baboon_nets_list <- list(bab_w1_fem,
                         bab_w2_fem,
                         bab_w3_a_fem,
                         bab_w3_b_fem,
                         bab_w4_a_fem,
                         bab_w4_b_fem, 
                         bab_new_fem)

# Let's add ages to these vertices
for (i in 1:length(baboon_nets_list)) {
  if (i == 1) {
    current_time <- as.Date("1997-12-31")
  } else if (i == 2) {
    current_time <- as.Date("1998-12-31")
  } else if (i == 3 | i == 4) {
    current_time <- as.Date("2000-08-01")
  } else if ( i == 5 | i == 6 ){
    current_time <- as.Date("2001-08-16")
  } else {
    current_time <- as.Date("1999-07-31")
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
bab_new_age <- baboon_nets_list[[7]]



#### Elephant data ####

# Import the data 
elephant_data <- read_csv("Data/dist.matrix.t1.csv")
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
elephant_w1 <- 
  graph_from_adjacency_matrix(inv_elephant_matrix, mode = "undirected", weighted = T)

# Network for wave 2 elephant data 
# Import the data 
elephant_data_w2 <- read_csv("Data/dist.matrix.t2.csv")
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
elephant_w2<- 
  graph_from_adjacency_matrix(inv_elephant_matrix_w2, mode = "undirected", weighted = T)

# Network elephant data wave 3 
# Import the data 
elephant_data_w3 <- read_csv("Data/dist.matrix.t3.csv")
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
elephant_w3<- 
  graph_from_adjacency_matrix(inv_elephant_matrix_w3, mode = "undirected", weighted = T)


elephant_ages <- read_csv("Data/SamburuAges_Nico.csv") %>%
  select(2:3)
elephant_ages2 <- read_csv("Data/elephant_ages_1.csv")
elephant_ages3 <- read_csv("Data/elephant_ages_2.csv")

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
        "Data/networks/dolphin_w1.rds")
saveRDS(dolphin_w2, 
        "Data/networks/dolphin_w2.rds")
saveRDS(dolphin_w3, 
        "Data/networks/dolphin_w3.rds")
saveRDS(dolphin_w4, 
        "Data/networks/dolphin_w4.rds")
saveRDS(dolphin_w5, 
        "Data/networks/dolphin_w5.rds")
saveRDS(dolphin_w6, 
        "Data/networks/dolphin_w6.rds")

saveRDS(elph_w1_age, 
        "Data/networks/elph_w1_age.rds")
saveRDS(elph_w2_age, 
        "Data/networks/elph_w2_age.rds")
saveRDS(elph_w3_age, 
        "Data/networks/elph_w3_age.rds")


saveRDS(bab_w1_age, 
        "Data/networks/bab_w1_age.rds")
saveRDS(bab_w2_age, 
        "Data/networks/bab_w2_age.rds")
saveRDS(bab_w3_a_age, 
        "Data/networks/bab_w3_a_age.rds")
saveRDS(bab_w3_b_age, 
        "Data/networks/bab_w3_b_age.rds")
saveRDS(bab_w4_a_age, 
        "Data/networks/bab_w4_a_age.rds")
saveRDS(bab_w4_b_age, 
        "Data/networks/bab_w4_b_age.rds")
saveRDS(bab_new_age, 
        "Data/networks/bab_new_age.rds")




