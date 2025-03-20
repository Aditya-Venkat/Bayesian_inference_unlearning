library(mvtnorm)
library(pracma)

set.seed(1)
#simulate the data 
theta <- rnorm(1e2, 2,1)
prior <- function(theta){
  return(dnorm())
}
#posterior of param
post_samp <- function(num){
  
  n <- length(S)
  mu <- (2 + sum(S))/(n + 1)
  sigma2 <- 1/(n+1)
  
  return(rnorm(num, mu, sigma2))
  
}

post_samp_unlearned <- function(num){
  n <- length(S_rem)
  mu <- (2 + sum(S_rem))/(n + 1)
  sigma2 <- 1/(n+1)
  
  return(rnorm(num, mu, sigma2))
}


S <- rnorm(1e3, 2, 1)
i <- sample(1:length(S) , 50)
S2 <- S[i]
S_rem <- S[-i]



#unlearning function
influence <- function(){
  exp_theta <-post_samp(1e3)
  
  num <- sapply(S2 , function(i){sum(i- exp_theta)/length(exp_theta)})
  num <- sum(num)
  denom <- -(length(S)+1)
  
  
  
  return(-num/denom)
  
  
}

library(entropy)
final_theta = post_samp(1e3) - influence()


density_final <- density(final_theta)
density_retrained <- density(post_samp_unlearned(1e3))


plot(density_final, col = "blue", lwd = 2, main = "Comparison of Models", xlab = "Theta", ylab = "Density")
lines(density_retrained, col = "red", lwd = 2)

legend("topright", legend = c("Unlearned Model", "Retrained Model"), col = c("blue", "red"), lwd = 2)
target <- post_samp_unlearned(length(final_theta))

bins <- seq(min(final_theta, target), max(final_theta, target), length.out = 50)

P_hist <- hist(final_theta, breaks = bins, plot = FALSE)$counts
Q_hist <- hist(target, breaks = bins, plot = FALSE)$counts

P_prob <- P_hist / sum(P_hist)
Q_prob <- Q_hist / sum(Q_hist)

KL.empirical(P_prob + 1e-6, Q_prob + 1e-6)
