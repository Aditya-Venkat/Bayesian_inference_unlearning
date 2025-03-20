library(mvtnorm)
set.seed(1)

u <- rmvnorm(4 , mean = c(0,0) , sigma = diag(2))

sample_means <- sample(1:4, 2e3,replace = TRUE)
sample_means <- u[sample_means,]

zs <- t(mapply(function(mu1, mu2) rmvnorm(1, mean = c(mu1, mu2), sigma = diag(2)), 
               sample_means[,1], sample_means[,2]))


#HMC fitting
library(rstan)

data_list <- list(
  N = nrow(zs),  # Number of data points
  K = 4,         # Number of clusters (since we used 4 means u)
  zs = zs        # Data points
)

fit <- stan(file = "gmm_model.stan", data = data_list, iter = 1000, chains = 4)


posterior_samples <- extract(fit)


log_lik <- function(z, theta, mu, sigma) {
  z <- matrix(z , nrow = 1, ncol = 2)  # Converts to a numeric vector if needed
  
  loglik <- sum(sapply(1:4, function(i) {
    dmvnorm(z, mean = mu[i, , drop= FALSE], sigma = sigma[i, , ], log = FALSE) * theta[i]
  }))
  
  return(loglik)
}


log_mu_prior <- function(mu){
  mu <- matrix(mu, nrow = 1 , ncol = 2)
  return(log(dmvnorm(mu, rep(0,2),diag(25,2))))
}

library(MCMCpack) 

log_sigma_prior <- function(sigma){
  df <- 3  
  scale_matrix <- diag(2)  
  
  return(log(dwish(W = sigma, v = df, S = scale_matrix)))  
}
posterior <- function(data,theta , mu, sigma){
  
  sum_log_lik <- sum(sapply(1:dim(data)[1],function(i)log_lik(data[i,], theta , mu , sigma)))
  sum_log_mu_prior <- sum(sapply(1:4 , function(i)log_mu_prior(mu[i,])))
  sum_log_sigma_prior <- sum(sapply(1:4 ,function(i) log_sigma_prior(sigma[i,,])))
  
  return(sum_log_lik + sum_log_mu_prior + sum_log_sigma_prior)
  
  
}

library(numDeriv)


z <- zs[1, ]  
theta <- posterior_samples$theta[1, ] 
mu <- posterior_samples$mu[1, , ] 
sigma <- posterior_samples$Sigma[1, , , ]  
posterior(zs, theta, mu ,sigma)
library(numDeriv)

flatten_sigma <- function(sigma) {
  return(as.vector(sapply(1:4, function(i) c(sigma[i, 1, 1], sigma[i, 1, 2], sigma[i, 2, 2]))))
}

reconstruct_sigma <- function(sigma_vec) {
  sigma <- array(0, dim = c(4, 2, 2))
  
  for (i in 1:4) {
    idx <- (i - 1) * 3 + 1  
    sigma[i, 1, 1] <- sigma_vec[idx]
    sigma[i, 1, 2] <- sigma_vec[idx + 1]
    sigma[i, 2, 1] <- sigma_vec[idx + 1]
    sigma[i, 2, 2] <- sigma_vec[idx + 2]
  }
  return(sigma)
}


grad_log_lik <- function(theta, mu, sigma, z) {
  sigma_flattened <- flatten_sigma(sigma)  
  
  param_vector <- c(theta, as.vector(mu), sigma_flattened)
  
  log_lik_wrapper <- function(params) {
    theta_new <- params[1:4]
    mu_new <- matrix(params[5:12], nrow = 4, ncol = 2)
    sigma_new <- reconstruct_sigma(params[13:24])  
    
    return(log_lik(z, theta_new, mu_new, sigma_new))
  }
  
  grad_vector <- grad(log_lik_wrapper, param_vector)
  return(grad_vector)
}


hessian_posterior <- function(data, theta, mu, sigma) {
  sigma_flattened <- flatten_sigma(sigma) 
  param_vector <- c(theta, as.vector(mu), sigma_flattened)
  
  posterior_wrapper <- function(params) {
    theta_new <- params[1:4]
    mu_new <- matrix(params[5:12], nrow = 4, ncol = 2)
    sigma_new <- reconstruct_sigma(params[13:24]) 
    
    return(posterior(data, theta_new, mu_new, sigma_new))
  }
  
  hessian_matrix <- hessian(posterior_wrapper, param_vector)
  
  return(hessian_matrix)
}


ind <- sample(1:2e3, 1e2,replace = FALSE)
S2 <- zs[ind,]
S_red <- zs[-ind,]

influence <- function(data){
  
  theta <- posterior_samples$theta
  mu <- posterior_samples$mu
  sigma <- posterior_samples$Sigma
  
  exp_grad <- function(z){
    sum(sapply(1:dim(theta)[1] , function(i)grad_log_lik(theta[i,] , mu[i,,] , sigma[i,,,] , z) ))/dim(theta)[1]
  }
  
  num <- sum(sapply(1:dim(data)[1], function(i) exp_grad(data[i,])))
  
  denom <- sum(sapply(1:dim(theta)[1] , function(i)hessian_posterior(zs, theta[i,] , mu[i,,] , sigma[i,,,] ) ))/dim(theta)[1]
  
  out <- -solve(denom)%*%num
  return(out)
  
}
I <- influence(S2)

#Retraining model

data_list <- list(
  N = nrow(S_red),  
  K = 4,         
  zs = S_red        
)
fit2 <- stan(file = "gmm_model.stan", data = data_list, iter = 500, chains = 4)
save(fit2, file = "gmm_model_retrained.stan")
