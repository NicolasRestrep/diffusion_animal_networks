
net <- dolphin_w1
e <- 1
r_max <- 500
sim <- 1

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
    nei <- adjm[i,] > 0
    # If you dont want to include weights, quote above, unquote below
    #nei <- adjm_bool[i,]
    # Proceed if there is at least one neighbour
    if(sum(nei) > 0){
      # Simple contagion for e = 1 and complex contagion for e = 2
      if(runif(n = 1, min = 0, max = 1) <= (sum(adjm[i,][info])/sum(nei))^e){
        info[i] <- TRUE
      }
    }
  }
  # Record proportion of the population with info
  proportion[r] <- sum(info) / N
  # Increment the round counter
  r <- r + 1
}
