#### Most influential nodes ####
library(tidyverse)

most_influential_nodes <- function(path) {
d <- read_csv(path)
inf_nodes <- d %>% 
  arrange(change_efficiency) %>% 
  slice(1:10) %>% 
  pull(node_name)
return(inf_nodes)
}

list_paths <- c("Data/elephant_removal_w1.csv",
                "Data/elephant_removal_w2.csv", 
                "Data/elephant_removal_w3.csv", 
                "Data/dolphin_removal_w1.csv", 
                "Data/dolphin_removal_w2.csv",
                "Data/dolphin_removal_w3.csv", 
                "Data/dolphin_removal_w4.csv",
                "Data/dolphin_removal_w5.csv", 
                "Data/dolphin_removal_w6.csv")

list_nodes <- map(list_paths, 
    most_influential_nodes)

df_inf <- matrix(NA, nrow = 9, ncol = 12)

df_inf[,1] <- c(rep("elephant", 3), rep("dolphin", 6))
df_inf[,2] <- c(1:3, 1:6)
for (i in 1:9) {
  df_inf[i, 3:12] <- list_nodes[[i]]
}

df_inf <- data.frame(df_inf)
names(df_inf) <- c("species", "wave", as.character(1:10))
write.csv(df_inf, "most_influential_nodes.csv")
