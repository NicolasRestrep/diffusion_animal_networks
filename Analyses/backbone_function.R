# Build the function 
backbone_graph <- function(alpha, simple, g) {
  
# All edges and their weight 
e <- cbind(igraph::as_data_frame(g)[, 1:2 ], 
           weight = E(g)$weight)

# Add graph strength and in degree of each source to graph dataframe
# in
# Sum all edge weigths for each node 
w_in <- graph.strength(g,
                       mode = "in")
w_in <- data.frame(to = names(w_in), 
                   w_in, stringsAsFactors = FALSE)
# Get the degree for each node 
k_in <- degree(g, mode = "in")
k_in <- data.frame(to = names(k_in), 
                   k_in, 
                   stringsAsFactors = FALSE)

# Join them both 
e_in <- e %>%
  left_join(w_in, by = "to") %>%
  left_join(k_in, by = "to") %>%
  mutate(alpha_in = (1-(weight/w_in))^(k_in-1))

# Same for out degree 
w_out <- graph.strength(g, mode = "out")
w_out <- data.frame(from = names(w_out), w_out, stringsAsFactors = FALSE)
k_out <- degree(g, mode = "out")
k_out <- data.frame(from = names(k_out), k_out, stringsAsFactors = FALSE)

e_out <- e %>%
  left_join(w_out, by = "from") %>%
  left_join(k_out, by = "from") %>%
  mutate(alpha_out = (1-(weight/w_out))^(k_out-1))

# Join everything and add alpha to the network
#full
e_full <- left_join(e_in, 
                    e_out, 
                    by = c("from", "to", "weight"))

e_full <- e_full %>%
  mutate(alpha = ifelse(alpha_in < alpha_out, alpha_in, alpha_out)) %>%
  select(from, to, alpha)

E(g)$alpha <- e_full$alpha

# Finally prune the network 
graph_pruned <- delete.edges(g, 
                             which(E(g)$alpha >= alpha))
graph_pruned_deleted <- delete.vertices(graph_pruned, 
                                        which(degree(graph_pruned) == 0))

# Count number of vertices
num_e <- gsize(graph_pruned)
num_v <- gorder(graph_pruned_deleted)

res <- c(num_e, num_v)

# return statement 
if (simple == T) {
  return(res)
} else {
  return(list(res, graph_pruned))
}
}

# Now iterate trough different values of alpha 
alphas <- seq(from = 0, to = 1, by = 0.05)

result_e <- vector("list", length(alphas))
result_v <- vector("list", length(alphas))
i = 1

for (alpha in alphas){
  temp_res <- backbone_graph(alpha = alpha, 
                             simple = T, 
                             g = dolphin_w6)  
  result_e[i]<- temp_res[1]
  result_v[i]<- temp_res[2]
  i <- i + 1
}

res_df <- as.data.frame(result_e)
colnames(res_df) <- alphas
res_df <- res_df %>% 
  gather() %>% 
  mutate(pct_edges = value / (gsize(dolphin_w6))) %>% 
  rename(alpha = key, num_edges = value)

ggplot(res_df, aes(alpha, num_edges)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label= round(pct_edges,2)), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal()

