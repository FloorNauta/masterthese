#### Master thesis Floor Nauta
#code for creating Figure 7 with estimated edge weights

library(bootnet)
library(dplyr)
library(qgraph)

load('true_networks.rdata')


#simulation ----

parSim(
  
  #conditions
  n = c(500, 1000, 2500),
  p_conn = c(0.1, 0.3, 0.6),
  
  #setup
  reps = 10,
  nCores = 1,
  write = TRUE,
  name = "TEST_weights",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
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
    
    results <- list(est_weights = as.vector(est_weights),
                    true_weights = as.vector(true_weights))
    return(results)
  }
)


TEST_weights_results <- read.table("TEST_weights.txt", header = TRUE)



# plot ----

TEST_weights_results <- transform(TEST_weights_results, 
                                  n = factor(n, levels = c(500, 1000, 2500),
                                             labels = c("n = 500", "n = 1000", "n = 2500")),
                                  p_conn = factor(p_conn, levels = c(0.1, 0.3, 0.6),
                                  labels = c("Low", "Medium", "High")))


ggplot(TEST_weights_results, aes(x = est_weights, y = true_weights, colour = as.factor(p_conn))) +
  geom_point(aes(shape = as.factor(n))) +
  ylim(-0.4, 0.4) +
  xlim(-0.4, 0.4) +
  labs(x = "Estimated edge weights", y = "True edge weights", colour = "Connectivity\nlevel", shape = "Sample size") +
  theme_bw() +
  facet_grid(n ~ p_conn)


#plot showing false positives
TEST_weights_results <- mutate(TEST_weights_results,
                               false_pos = if_else(true_weights == 0 & est_weights != 0, 1, 0))

TEST_weights_results <- transform(TEST_weights_results, 
                                  false_pos = factor(false_pos, levels = c(0, 1),
                                             labels = c("No", "Yes")))

ggplot(TEST_weights_results, aes(x = est_weights, y = true_weights, colour = as.factor(false_pos))) +
  geom_point() +
  ylim(-0.4, 0.4) +
  xlim(-0.4, 0.4) +
  labs(x = "Estimated edge weights", y = "True edge weights", colour = "False positive?") +
  scale_colour_manual(values = c(1, 2)) +
  theme_bw() +
  facet_grid(n ~ p_conn)
