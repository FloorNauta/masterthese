#### Master thesis Floor Nauta
#function used in simulation




compare_results <- function(est_weights, true_weights) {
  cor0 <- function(x, y, ...) {
    if (sum(!is.na(x)) < 2 ||
        sum(!is.na(y)) < 2 || 
        sd(x, na.rm = TRUE) == 0 |
        sd(y, na.rm = TRUE) == 0) {
      return(0)
    } else {
      return(cor(as.vector(x[upper.tri(x)]),as.vector(y[upper.tri(y)])))
    }
  }

  
  MAE <- function(x, y){ 
    mean( abs(x - y), na.rm = TRUE)
  }
  
  
  #true positives
  true_pos <- sum(true_weights != 0 & est_weights != 0)
  
  #true negatives
  true_neg <- sum(true_weights == 0 & est_weights == 0)
  
  #false positives
  false_pos <- sum(true_weights == 0 & est_weights != 0)
  
  #false negatives
  false_neg <- sum(true_weights != 0 & est_weights == 0)
  
  out <- list()
  
  #sensitivity
  out$sensitivity <- true_pos / (true_pos + false_neg)
  
  #specificity
  out$specificity <- true_neg / (true_neg + false_pos)
  
  #correlation
  out$correlation <- cor0(est_weights, true_weights)

  #MAE
  out$MAE <- MAE(est_weights, true_weights)

  #mean and sd of adjacency matrix
  out$M <- mean(est_weights)
  out$SD <- sd(est_weights)
  
  return(out)
}
