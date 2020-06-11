#### Master thesis Floor Nauta
#code for creating true networks;
#simulating random, scale free and small world networks;
#example plot;
#getting data ready for analysis

#loading packages
library(dplyr)
library(modelr)
library(parSim)
library(tidyr)
library(bootnet)
library(qgraph)
library(igraph)

# create true networks ----

#make a function to generate edge weights
make_weighted_adj_matrix <- function(Graph){
  Graph <- as.matrix(get.adjacency(Graph)) * rnorm(20^2, 0.05, 0.1) #normal distribution mean=0.05, sd=0.1
  Graph[lower.tri(Graph)] <- t(Graph)[lower.tri(Graph)]
  Graph <- round(Graph, 8)
  Graph
} 


##RANDOM##
random_temp <- matrix(NA, 1, 20)
p_conn <- c( 0.023, 0.1, 0.3, 0.6)

set.seed(111)

for (i in 1:4) {
  graph_random <- erdos.renyi.game(20, p_conn[i], type = "gnp")
  #use function to get weighted adjacency matrix
  random_temp <- rbind(random_temp, make_weighted_adj_matrix(graph_random))
}

random_temp <- na.omit(random_temp) #remove first row with NA's

#make a data frame with network type, connectivity level and all weights
weights_random <- data.frame(
  network_type = "random",
  conn_level = c(rep(0.023, 20), rep(0.1, 20), rep(0.3, 20), rep(0.6, 20)),
  w = rep(c(1:20), 4),
  weights = random_temp
)



##SCALE FREE##
scalefree_temp <- matrix(NA, 1, 20)
m <-  c(1, 3, 7)

set.seed(222)

for (i in 1:3) {
  graph_scalefree <- barabasi.game(20, power = 1, m = m[i], directed = FALSE) 
  scalefree_temp <- rbind(scalefree_temp, make_weighted_adj_matrix(graph_scalefree))
}
scalefree_temp <- na.omit(scalefree_temp) #remove first row with NA's

#make a data frame with network type, connectivity level and all weights
weights_scalefree <- data.frame(
  network_type = "scalefree",
  conn_level = c(rep(1, 20), rep(3, 20), rep(7, 20)),
  w = rep(c(1:20), 3),
  weights = scalefree_temp
)


##SMALL WORLD##
smallworld_temp <- matrix(NA, 1, 20)
n_neighbours <-  c(1, 3, 6) #nr of neighbours

set.seed(333)

for (i in 1:3) {
  graph_smallworld <- watts.strogatz.game(1, 20, n_neighbours[i], 0.1)
  smallworld_temp <-rbind(smallworld_temp,  make_weighted_adj_matrix(graph_smallworld))
}
smallworld_temp <- na.omit(smallworld_temp) #remove first row with NA's

#make a data frame with network type, connectivity level and all weights
weights_smallworld <- data.frame(
  network_type = "smallworld",
  conn_level = c(rep(1, 20), rep(3, 20), rep(6, 20)),
  w = rep(c(1:20), 3),
  weights = smallworld_temp
)

##ALL NETWORK TYPES##
#combine all network types & save
true_networks <- rbind(weights_random, weights_scalefree, weights_smallworld)
save(true_networks, file = 'true_networks.rdata')
load('true_networks.rdata')



# simulation RANDOM network ----

library(parSim)
library(bootnet)
load('true_networks.rdata')


parSim(
  
  #conditions
  n = c(100, 250, 500, 1000, 2500),
  p_conn = c(0.023, 0.1, 0.3, 0.6),
  
  #setup
  reps = 1000,
  nCores = 1,
  write = TRUE,
  name = "sim_results_random",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    source("compare_results.R")
    load("true_networks.rdata")
    
    #data
    true_weights <- true_networks %>%
      filter(conn_level == p_conn) %>%
      select(weights.1:weights.20) 
    true_weights <- as.matrix(true_weights)
    Sigma <- cov2cor(solve(diag(ncol(true_weights)) - true_weights)) #invert correlation to covariance
    #generate data from multivariate normal distribution
    Data_sim <- mvtnorm::rmvnorm(n, sigma = Sigma) #sigma = covariance matrix
    
    #network based on simulation data
    Network_sim <- estimateNetwork(Data_sim, default = "EBICglasso", refit = TRUE) 
    est_weights <- getWmat(Network_sim)
    
    #outcome measures
    results <- compare_results(est_weights, true_weights)
    
    #results
    return(results)
    
    
  }
  
  
)


results_random <- read.table("sim_results_random.txt", header = TRUE)



# simulation SCALE FREE network ----

library(parSim)
library(bootnet)
load('true_networks.rdata')


parSim(
  
  #conditions
  n = c(100, 250, 500, 1000, 2500),
  m = c(1, 3, 7),
  
  #setup
  reps = 1000,
  nCores = 1,
  write = TRUE,
  name = "sim_results_scalefree",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    source("compare_results.R")
    load("true_networks.rdata")
    
    #data
    true_weights <- true_networks %>%
      filter(network_type == "scalefree", conn_level == m) %>%
      select(weights.1:weights.20) 
    true_weights <- as.matrix(true_weights)
    Sigma <- cov2cor(solve(diag(ncol(true_weights)) - true_weights)) #invert correlation to covariance
    #generate data from multivariate normal distribution
    Data_sim <- mvtnorm::rmvnorm(n, sigma = Sigma) #sigma = covariance matrix
    
    #network based on simulation data
    Network_sim <- estimateNetwork(Data_sim, default = "EBICglasso", refit = TRUE) 
    est_weights <- getWmat(Network_sim)
    
    #outcome measures
    results <- compare_results(est_weights, true_weights)
    
    #results
    return(results)
    
    
  }
  
  
)


results_scalefree <- read.table("sim_results_scalefree.txt", header = TRUE)




# simulation SMALL WORLD network ----

library(parSim)
library(bootnet)
load('true_networks.rdata')


parSim(
  
  #conditions
  n = c(100, 250, 500, 1000, 2500),
  n_neighbours = c(1, 3, 6),
  
  #setup
  reps = 1000,
  nCores = 1,
  write = TRUE,
  name = "sim_results_smallworld",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    source("compare_results.R")
    load("true_networks.rdata")
    
    #data
    true_weights <- true_networks %>%
      filter(network_type == "smallworld", conn_level == n_neighbours) %>%
      select(weights.1:weights.20) 
    true_weights <- as.matrix(true_weights)
    Sigma <- cov2cor(solve(diag(ncol(true_weights)) - true_weights)) #invert correlation to covariance
    #generate data from multivariate normal distribution
    Data_sim <- mvtnorm::rmvnorm(n, sigma = Sigma) #sigma = covariance matrix
    
    #network based on simulation data
    Network_sim <- estimateNetwork(Data_sim, default = "EBICglasso", refit = TRUE) 
    est_weights <- getWmat(Network_sim)
    
    #outcome measures
    results <- compare_results(est_weights, true_weights)
    
    #results
    return(results)
    
    
  }
  
  
)


results_smallworld <- read.table("sim_results_smallworld.txt", header = TRUE)



#PLOT what the networks look like ----
load('true_networks.rdata')


#RANDOM
#unweighted
qgraph(true_networks[21:40, 4:23], vsize = 4, esize = 7, layout = "spring", 
       labels = FALSE, theme = "colorblind", weighted = FALSE) # = p_conn 0.100
#weighted
qgraph(true_networks[21:40, 4:23], vsize = 4, esize = 10, layout = "spring", 
       labels = FALSE, theme = "colorblind") # = p_conn 0.100


#SCALE FREE
#unweighted
qgraph(true_networks[81:100, 4:23], vsize = 4, esize = 7, layout = "spring", 
       labels = FALSE, theme = "colorblind", weighted = FALSE) # = m 1
#weighted
qgraph(true_networks[81:100, 4:23], vsize = 4, esize = 10, layout = "spring", 
       labels = FALSE, theme = "colorblind") # = m 1


#SMALL WORLD
#unweighted
qgraph(true_networks[141:160, 4:23], vsize = 4, esize = 7, layout = "spring", 
       labels = FALSE, theme = "colorblind", weighted = FALSE) # = n_neighbours 1
#weighted
qgraph(true_networks[141:160, 4:23], vsize = 4, esize = 10, layout = "spring", 
       labels = FALSE, theme = "colorblind") # = n_neighbours 1




#getting data ready for analysis ----

results_random2 <- mutate(results_random, network_type = "random",
                          conn_level = ifelse(p_conn == 0.023, 1,
                                              ifelse(p_conn == 0.1, 2,
                                                     ifelse(p_conn == 0.3, 3, 4)))) %>%
  select(network_type, conn_level, n, sensitivity:MAE)

results_scalefree2 <- mutate(results_scalefree, network_type = "scalefree",
                             conn_level = ifelse(m == 1, 2,
                                                 ifelse(m == 3, 3, 4))) %>%
  select(network_type, conn_level, n, sensitivity:MAE)

results_smallworld2 <- mutate(results_smallworld, network_type = "smallworld",
                              conn_level = ifelse(n_neighbours == 1, 2,
                                                  ifelse(n_neighbours == 3, 3, 4))) %>%
  select(network_type, conn_level, n, sensitivity:MAE)


results_total <- rbind(results_random2, results_scalefree2, results_smallworld2)
results_total$n <- as.factor(results_total$n)
results_total$conn_level <- as.factor(results_total$conn_level)

# save(results_total, file = 'results_total.rdata')
# load('results_total.rdata')


#releveling data for analyses
results_releveled <- results_total
results_releveled$n  <- relevel(results_releveled$n, ref = "2500") # set 2500 as reference
results_releveled$conn_level  <- relevel(results_releveled$conn_level, ref = 4) # set level 4 (= 0.6 for random) as reference
# save(results_releveled, file = 'results_releveled.rdata')
# load('results_releveled.rdata')

