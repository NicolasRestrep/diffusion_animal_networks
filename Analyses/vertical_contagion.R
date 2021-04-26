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

