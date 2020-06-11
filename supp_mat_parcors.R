#### Master thesis Floor Nauta
#code supplementary material;
#simulation study on size partial correlations
#plots

#loading packages
library(dplyr)
library(modelr)
library(parSim)
library(tidyr)
library(bootnet)
library(qgraph)
library(igraph)
library(BDgraph)

#TEST1 size parcors rnorm m = 0, sd = 0.15 ----
# parcors between approximately -0.4 and 0.4 (like williams), but normally distributed

# create true networks

#make a function to generate edge weights
make_weighted_adj_matrix <- function(Graph){
  Graph <- as.matrix(get.adjacency(Graph)) * rnorm(20^2, 0, 0.15) #normal distribution mean=0, sd=0.15
  Graph[lower.tri(Graph)] <- t(Graph)[lower.tri(Graph)]
  Graph <- round(Graph, 8)
  Graph
} 


##RANDOM##
random_temp <- matrix(NA, 1, 20)
p_conn <- c( 0.023, 0.1, 0.3, 0.6)

set.seed(222)

for (i in 1:4) {
  graph_random <- erdos.renyi.game(20, p_conn[i], type = "gnp")
  #use function to get weighted adjacency matrix
  random_temp <- rbind(random_temp, make_weighted_adj_matrix(graph_random))
}

random_temp <- na.omit(random_temp) #remove first row with NA's

#test min & max parcors
min(random_temp)
max(random_temp)


#make a data frame with network type, connectivity level and all weights
weights_random1 <- data.frame(
  network_type = "random",
  conn_level = c(rep(0.023, 20), rep(0.1, 20), rep(0.3, 20), rep(0.6, 20)),
  w = rep(c(1:20), 4),
  weights = random_temp
)


# simulation RANDOM network

parSim(
  
  #conditions
  n = c(100, 250, 500, 1000, 2500),
  p_conn = c(0.023, 0.1, 0.3, 0.6),
  
  #setup
  reps = 100,
  nCores = 1,
  write = TRUE,
  name = "TEST1_parcors",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    source("compare_results.R")
    
    
    #data
    true_weights <- weights_random1 %>%
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


TEST1_results <- read.table("TEST1_parcors.txt", header = TRUE)

#PLOT TEST1----

#sensitivity
ggplot(TEST1_results, aes(x = factor(n), y = sensitivity, fill = factor(p_conn))) +
  geom_boxplot() +
  ylim(0, 1)

#specificity
ggplot(TEST1_results, aes(x = factor(n), y = specificity, fill = factor(p_conn))) +
  geom_boxplot() +
  ylim(0, 1)

#specificity is low for high conn_level, still not as low as williams


#make a nice plot
#average outcome 
outcome_TEST1 <- TEST1_results %>%
  group_by(n, p_conn) %>%
  summarise(average_sen = mean(sensitivity),
            sd_sen = sd(sensitivity),
            average_spec = mean(specificity),
            sd_spec = sd(specificity))#,
            # average_cor = mean(correlation),
            # sd_cor = sd(correlation),
            # average_MAE = mean(MAE),
            # sd_MAE = sd(MAE))
#outcome_TEST1

average_TEST1 <- outcome_TEST1 %>% 
  pivot_longer(c('average_sen', 'average_spec'), 
               names_to = "outcome_measure", values_to = "outcome_value")

sd_TEST1 <- outcome_TEST1 %>%
  pivot_longer(c('sd_sen', 'sd_spec'), 
               names_to = "sd_measure", values_to = "sd_value")

average_sd_TEST1 <- cbind(as.data.frame(average_TEST1), sd = sd_TEST1$sd_value)
average_sd_TEST1 <- mutate(average_sd_TEST1, ymax = outcome_value + sd,
                          ymin = outcome_value - sd)

#values between 0 and 1
average_sd_TEST1$ymax[average_sd_TEST1$ymax > 1] <- 1
average_sd_TEST1$ymin[average_sd_TEST1$ymin < 0] <- 0

average_sd_TEST1 <- transform(average_sd_TEST1, 
                             outcome_measure = factor(outcome_measure, levels = c('average_sen', 'average_spec'), 
                                                      labels = c("Sensitivity", "Specificity")),
                             p_conn = factor(p_conn, levels = c(0.023, 0.1, 0.3, 0.6),
                                                 labels = c("Lowest", "Low", "Medium", "High")))


ggplot(average_sd_TEST1, aes(x = as.factor(n), y = outcome_value, colour = as.factor(p_conn), group = as.factor(p_conn))) +
  # geom_col(position = "dodge", colour = 1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = as.factor(p_conn)), alpha = 0.15) +  
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  labs(x = "Sample size", y = "", colour = "Connectivity\nlevel", fill = "Connectivity\nlevel") +
  #  scale_colour_manual(values = c(3, 4, 5, 6)) +
  theme_bw() +
  facet_grid(. ~ outcome_measure)


# TEST1 weights + plot----
#simulation saving edge weights
parSim(
  
  #conditions
  n = c(500, 1000, 2500),
  p_conn = c(0.1, 0.3, 0.6),
  
  #setup
  reps = 10,
  nCores = 1,
  write = TRUE,
  name = "TEST1_weights",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)

    
    #data
    true_weights <- weights_random1 %>%
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


TEST1_weights_results <- read.table("TEST1_weights.txt", header = TRUE)


#make false_pos variable
TEST1_weights_results <- mutate(TEST1_weights_results,
                               false_pos = if_else(true_weights == 0 & est_weights != 0, 1, 0))

#relabel variables
TEST1_weights_results <- transform(TEST1_weights_results, 
                                  n = factor(n, levels = c(500, 1000, 2500),
                                             labels = c("n = 500", "n = 1000", "n = 2500")),
                                  p_conn = factor(p_conn, levels = c(0.1, 0.3, 0.6),
                                                  labels = c("Low", "Medium", "High")),
                                  false_pos = factor(false_pos, levels = c(0, 1),
                                                     labels = c("No", "Yes")))

#plot
ggplot(TEST1_weights_results, aes(x = est_weights, y = true_weights, colour = as.factor(false_pos))) +
  geom_point() +
  ylim(-0.5, 0.5) +
  xlim(-0.5, 0.5) +
  labs(x = "Estimated edge weights", y = "True edge weights", colour = "False positive?") +
  scale_colour_manual(values = c(1, 2)) +
  theme_bw() +
  facet_grid(n ~ p_conn)


# TEST2 g-wishart distribution 20 df ----
#use rgwish() function with 20 degrees of freedom

library(BDgraph)

# create true networks
##RANDOM##
random_temp2 <- matrix(NA, 1, 20)
p_conn <- c( 0.023, 0.1, 0.3, 0.6)

set.seed(222)

for (i in 1:2) {
  graph_random2 <- erdos.renyi.game(20, p_conn[i], type = "gnp")
  adj <- as.matrix(get.adjacency(graph_random2))
  weights <- as.matrix(rgwish(n = 1, adj = adj, b = 20, D = diag(20, 20)))
  diag(weights) <- 0 #for some reason rgwish doesn't make diagonal 0
  random_temp2 <- rbind(random_temp2, round(weights, 8))
}

set.seed(500)
graph_random23 <- erdos.renyi.game(20, 0.3, type = "gnp")
adj <- as.matrix(get.adjacency(graph_random23))
weights <- as.matrix(rgwish(n = 1, adj = adj, b = 20, D = diag(20, 20)))
diag(weights) <- 0 #for some reason rgwish doesn't make diagonal 0
random_temp2 <- rbind(random_temp2, round(weights, 8))

set.seed(82)#82 min -.52 max .54     101 min -.56 max .64
graph_random26 <- erdos.renyi.game(20, 0.6, type = "gnp")
adj <- as.matrix(get.adjacency(graph_random26))
weights <- as.matrix(rgwish(n = 1, adj = adj, b = 20, D = diag(20, 20)))
diag(weights) <- 0 #for some reason rgwish doesn't make diagonal 0
random_temp2 <- rbind(random_temp2, round(weights, 8))

random_temp2 <- na.omit(random_temp2) #remove first row with NA's


#test min & max parcors
min(random_temp2)
max(random_temp2)
mean(random_temp2)
sd(random_temp2)

hist(random_temp2[1:20, ])
hist(random_temp2[21:40, ])
hist(random_temp2[41:60, ])
hist(random_temp2[61:80, ])
# qgraph(random_temp2[61:80, ], layout = "spring")

#make a data frame with network type, connectivity level and all weights
weights_random2 <- data.frame(
  network_type = "random",
  conn_level = c(rep(0.023, 20), rep(0.1, 20), rep(0.3, 20), rep(0.6, 20)),
  w = rep(c(1:20), 4),
  weights = random_temp2
)


# simulation RANDOM network

parSim(
  
  #conditions
  n = c(100, 250, 500, 1000, 2500),
  p_conn = c(0.023, 0.1, 0.3, 0.6),
  
  #setup
  reps = 100,
  nCores = 1,
  write = TRUE,
  name = "TEST2_parcors",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    source("compare_results.R")
    
    
    #data
    true_weights <- weights_random2 %>%
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


TEST2_results <- read.table("TEST2_parcors.txt", header = TRUE)

#PLOT TEST 2 ----

#sensitivity
ggplot(TEST2_results, aes(x = factor(n), y = sensitivity, fill = factor(p_conn))) +
  geom_boxplot() +
  ylim(0, 1)

#specificity
ggplot(TEST2_results, aes(x = factor(n), y = specificity, fill = factor(p_conn))) +
  geom_boxplot() +
  ylim(0, 1)



#specificity is extremely low (between .6 and .4, lower than williams) for conn_level 0.6

#make a nice plot
#average outcome 
outcome_TEST2 <- TEST2_results %>%
  group_by(n, p_conn) %>%
  summarise(average_sen = mean(sensitivity, na.rm = TRUE),
            sd_sen = sd(sensitivity, na.rm = TRUE),
            average_spec = mean(specificity, na.rm = TRUE),
            sd_spec = sd(specificity, na.rm = TRUE))


average_TEST2 <- outcome_TEST2 %>% 
  pivot_longer(c('average_sen', 'average_spec'), 
               names_to = "outcome_measure", values_to = "outcome_value")

sd_TEST2 <- outcome_TEST2 %>%
  pivot_longer(c('sd_sen', 'sd_spec'), 
               names_to = "sd_measure", values_to = "sd_value")

average_sd_TEST2 <- cbind(as.data.frame(average_TEST2), sd = sd_TEST2$sd_value)
average_sd_TEST2 <- mutate(average_sd_TEST2, ymax = outcome_value + sd,
                           ymin = outcome_value - sd)

#values between 0 and 1
average_sd_TEST2$ymax[average_sd_TEST2$ymax > 1] <- 1
average_sd_TEST2$ymin[average_sd_TEST2$ymin < 0] <- 0

average_sd_TEST2 <- transform(average_sd_TEST2, 
                              outcome_measure = factor(outcome_measure, levels = c('average_sen', 'average_spec'), 
                                                       labels = c("Sensitivity", "Specificity")),
                              p_conn = factor(p_conn, levels = c(0.023, 0.1, 0.3, 0.6),
                                              labels = c("Lowest", "Low", "Medium", "High")))


ggplot(average_sd_TEST2, aes(x = as.factor(n), y = outcome_value, colour = as.factor(p_conn), group = as.factor(p_conn))) +
  # geom_col(position = "dodge", colour = 1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = as.factor(p_conn)), alpha = 0.15) +  
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  labs(x = "Sample size", y = "", colour = "Connectivity\nlevel", fill = "Connectivity\nlevel") +
  #  scale_colour_manual(values = c(3, 4, 5, 6)) +
  theme_bw() +
  facet_grid(. ~ outcome_measure)


#look at nr of errors
TEST2_error <- filter(TEST2_results, error == TRUE)

ggplot(TEST2_error, aes(x = factor(n), fill = factor(p_conn))) +
  geom_bar() +
  ylim(0, 100)




# TEST2 weights + plot----
#simulation saving edge weights
parSim(
  
  #conditions
  n = c(500, 1000, 2500),
  p_conn = c(0.1, 0.3, 0.6),
  
  #setup
  reps = 10,
  nCores = 1,
  write = TRUE,
  name = "TEST2_weights",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    
    
    #data
    true_weights <- weights_random2 %>%
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


TEST2_weights_results <- read.table("TEST2_weights.txt", header = TRUE)


#make false_pos variable
TEST2_weights_results <- mutate(TEST2_weights_results,
                                false_pos = if_else(true_weights == 0 & est_weights != 0, 1, 0))

#relabel variables
TEST2_weights_results <- transform(TEST2_weights_results, 
                                   n = factor(n, levels = c(500, 1000, 2500),
                                              labels = c("n = 500", "n = 1000", "n = 2500")),
                                   p_conn = factor(p_conn, levels = c(0.1, 0.3, 0.6),
                                                   labels = c("Low", "Medium", "High")),
                                   false_pos = factor(false_pos, levels = c(0, 1),
                                                      labels = c("No", "Yes")))

#plot
ggplot(TEST2_weights_results, aes(x = est_weights, y = true_weights, colour = as.factor(false_pos))) +
  geom_point() +
  ylim(-1, 1) +
  xlim(-1, 1) +
  labs(x = "Estimated edge weights", y = "True edge weights", colour = "False positive?") +
  scale_colour_manual(values = c(1, 2)) +
  theme_bw() +
  facet_grid(n ~ p_conn)


#look at nr of errors
TEST2w_error <- filter(TEST2_weights_results, error == TRUE)

ggplot(TEST2w_error, aes(x = factor(n), fill = factor(p_conn))) +
  geom_bar() +
  ylim(0, 10)



#TEST3 size parcors rnorm m = 0, sd = 0.19 ----
# parcors between approximately -0.5 and 0.5

# create true networks

#make a function to generate edge weights
make_weighted_adj_matrix <- function(Graph){
  Graph <- as.matrix(get.adjacency(Graph)) * rnorm(20^2, 0, 0.19) #normal distribution mean=0, sd=0.19
  Graph[lower.tri(Graph)] <- t(Graph)[lower.tri(Graph)]
  Graph <- round(Graph, 8)
  Graph
} 


##RANDOM##
random_temp <- matrix(NA, 1, 20)
p_conn <- c( 0.023, 0.1, 0.3, 0.6)

set.seed(111)

for (i in 1:2) {
  graph_random <- erdos.renyi.game(20, p_conn[i], type = "gnp")
  #use function to get weighted adjacency matrix
  random_temp <- rbind(random_temp, make_weighted_adj_matrix(graph_random))
}

set.seed(101)
graph_random3 <- erdos.renyi.game(20, 0.3, type = "gnp")
#use function to get weighted adjacency matrix
random_temp <- rbind(random_temp, make_weighted_adj_matrix(graph_random3))

set.seed(101)
graph_random6 <- erdos.renyi.game(20, 0.6, type = "gnp")
#use function to get weighted adjacency matrix
random_temp <- rbind(random_temp, make_weighted_adj_matrix(graph_random6))

random_temp <- na.omit(random_temp) #remove first row with NA's

#test min & max parcors
min(random_temp)
max(random_temp)


#make a data frame with network type, connectivity level and all weights
weights_random3 <- data.frame(
  network_type = "random",
  conn_level = c(rep(0.023, 20), rep(0.1, 20), rep(0.3, 20), rep(0.6, 20)),
  w = rep(c(1:20), 4),
  weights = random_temp
)


# simulation RANDOM network

parSim(
  
  #conditions
  n = c(100, 250, 500, 1000, 2500),
  p_conn = c(0.023, 0.1, 0.3, 0.6),
  
  #setup
  reps = 100,
  nCores = 1,
  write = TRUE,
  name = "TEST3_parcors",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    source("compare_results.R")
    
    
    #data
    true_weights <- weights_random3 %>%
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


TEST3_results <- read.table("TEST3_parcors.txt", header = TRUE)





#PLOT TEST 3----
#sensitivity
ggplot(TEST3_results, aes(x = factor(n), y = sensitivity, fill = factor(p_conn))) +
  geom_boxplot() +
  ylim(0, 1)

#specificity
ggplot(TEST3_results, aes(x = factor(n), y = specificity, fill = factor(p_conn))) +
  geom_boxplot() +
  ylim(0, 1)

#specificity between .6 and .5 for conn_level 0.6

#make a nice plot
#average outcome 
outcome_TEST3 <- TEST3_results %>%
  group_by(n, p_conn) %>%
  summarise(average_sen = mean(sensitivity, na.rm = TRUE),
            sd_sen = sd(sensitivity, na.rm = TRUE),
            average_spec = mean(specificity, na.rm = TRUE),
            sd_spec = sd(specificity, na.rm = TRUE))#,
# average_cor = mean(correlation),
# sd_cor = sd(correlation),
# average_MAE = mean(MAE),
# sd_MAE = sd(MAE))
#outcome_TEST1

average_TEST3 <- outcome_TEST3 %>% 
  pivot_longer(c('average_sen', 'average_spec'), 
               names_to = "outcome_measure", values_to = "outcome_value")

sd_TEST3 <- outcome_TEST3 %>%
  pivot_longer(c('sd_sen', 'sd_spec'), 
               names_to = "sd_measure", values_to = "sd_value")

average_sd_TEST3 <- cbind(as.data.frame(average_TEST3), sd = sd_TEST3$sd_value)
average_sd_TEST3 <- mutate(average_sd_TEST3, ymax = outcome_value + sd,
                           ymin = outcome_value - sd)

#values between 0 and 1
average_sd_TEST3$ymax[average_sd_TEST3$ymax > 1] <- 1
average_sd_TEST3$ymin[average_sd_TEST3$ymin < 0] <- 0

average_sd_TEST3 <- transform(average_sd_TEST3, 
                              outcome_measure = factor(outcome_measure, levels = c('average_sen', 'average_spec'), 
                                                       labels = c("Sensitivity", "Specificity")),
                              p_conn = factor(p_conn, levels = c(0.023, 0.1, 0.3, 0.6),
                                              labels = c("Lowest", "Low", "Medium", "High")))


ggplot(average_sd_TEST3, aes(x = as.factor(n), y = outcome_value, colour = as.factor(p_conn), group = as.factor(p_conn))) +
  # geom_col(position = "dodge", colour = 1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = as.factor(p_conn)), alpha = 0.15) +  
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  labs(x = "Sample size", y = "", colour = "Connectivity\nlevel", fill = "Connectivity\nlevel") +
  #  scale_colour_manual(values = c(3, 4, 5, 6)) +
  theme_bw() +
  facet_grid(. ~ outcome_measure)

#look at nr of errors
TEST3_error <- filter(TEST3_results, error == TRUE)

ggplot(TEST3_error, aes(x = factor(n), fill = factor(p_conn))) +
  geom_bar() +
  ylim(0, 100)





# TEST3 weights + plot----
#simulation saving edge weights
parSim(
  
  #conditions
  n = c(500, 1000, 2500),
  p_conn = c(0.1, 0.3, 0.6),
  
  #setup
  reps = 10,
  nCores = 1,
  write = TRUE,
  name = "TEST3_weights",
  
  #simulation code
  expression = {
    
    #packages
    library(bootnet)
    library(dplyr)
    library(qgraph)
    
    
    #data
    true_weights <- weights_random3 %>%
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


TEST3_weights_results <- read.table("TEST3_weights.txt", header = TRUE)


#make false_pos variable
TEST3_weights_results <- mutate(TEST3_weights_results,
                                false_pos = if_else(true_weights == 0 & est_weights != 0, 1, 0))

#relabel variables
TEST3_weights_results <- transform(TEST3_weights_results, 
                                   n = factor(n, levels = c(500, 1000, 2500),
                                              labels = c("n = 500", "n = 1000", "n = 2500")),
                                   p_conn = factor(p_conn, levels = c(0.1, 0.3, 0.6),
                                                   labels = c("Low", "Medium", "High")),
                                   false_pos = factor(false_pos, levels = c(0, 1),
                                                      labels = c("No", "Yes")))

#plot
ggplot(TEST3_weights_results, aes(x = est_weights, y = true_weights, colour = as.factor(false_pos))) +
  geom_point() +
  ylim(-1, 1) +
  xlim(-1, 1) +
  labs(x = "Estimated edge weights", y = "True edge weights", colour = "False positive?") +
  scale_colour_manual(values = c(1, 2)) +
  theme_bw() +
  facet_grid(n ~ p_conn)



#look at nr of errors
TEST3w_error <- filter(TEST3_weights_results, error == TRUE)

ggplot(TEST3w_error, aes(x = factor(n), fill = factor(p_conn))) +
  geom_bar() +
  ylim(0, 10)
