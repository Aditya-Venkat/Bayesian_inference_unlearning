library(mvtnorm)
library(numDeriv)


u <- rmvnorm(4 , mean = c(0,0) , sigma = diag(2))

sample_means <- sample(1:4, 2e3,replace = TRUE)
sample_means <- u[sample_means,]

zs <- t(mapply(function(mu1, mu2) rmvnorm(1, mean = c(mu1, mu2), sigma = diag(2)), 
               sample_means[,1], sample_means[,2]))
# ===== Parameter utilities =====
flatten_params <- function(theta_raw, mu, sigma) {
  sigma_flat <- unlist(lapply(1:4, function(i) {
    s <- sigma[i,,]
    c(s[1,1], s[1,2], s[2,2])
  }))
  c(theta_raw, as.vector(mu), sigma_flat)
}

unflatten_params <- function(beta) {
  theta_raw <- beta[1:4]
  mu <- matrix(beta[5:12], ncol = 2, byrow = TRUE)
  sigma_flat <- matrix(beta[13:24], nrow = 4, byrow = TRUE)
  sigma <- array(0, dim = c(4,2,2))
  for (i in 1:4) {
    sigma[i,,] <- matrix(c(sigma_flat[i,1], sigma_flat[i,2], 
                           sigma_flat[i,2], sigma_flat[i,3]), nrow = 2)
  }
  theta <- softmax(theta_raw)
  list(theta = theta, mu = mu, sigma = sigma, theta_raw = theta_raw)
}

softmax <- function(x) {
  exp_x <- exp(x - max(x))
  exp_x / sum(exp_x)
}

# ===== Log-posterior and gradient =====
log_lik <- function(z, theta, mu, sigma) {
  z <- matrix(z, nrow = 1)
  loglik <- sum(sapply(1:4, function(i) {
    dmvnorm(z, mean = mu[i,], sigma = sigma[i,,], log = FALSE) * theta[i]
  }))
  return(log(loglik))
}

log_mu_prior <- function(mu) {
  mu <- matrix(mu, nrow = 1)
  log(dmvnorm(mu, mean = c(0, 0), sigma = diag(25, 2)))
}

log_sigma_prior <- function(sigma) {
  df <- 3
  scale_matrix <- diag(2)
  log(dwish(W = sigma, v = df, S = scale_matrix))
}

posterior <- function(beta, data) {
  p <- unflatten_params(beta)
  theta <- p$theta
  mu <- p$mu
  sigma <- p$sigma

  ll <- sum(apply(data, 1, function(z) log_lik(z, theta, mu, sigma)))
  mu_prior <- sum(sapply(1:4, function(i) log_mu_prior(mu[i,])))
  sigma_prior <- sum(sapply(1:4, function(i) log_sigma_prior(sigma[i,,])))
  return(ll + mu_prior + sigma_prior)
}

grad_posterior <- function(beta, data) {
  grad(func = function(b) posterior(b, data), x = beta)
}

# ===== L-BFGS memory-based Langevin Sampler =====
estimate_xi_eta <- function(gk, zk, s_list, y_list, M, gamma_k = 1.0) {
  # Initialize B0, S0, G0
  d <- length(gk)
  Bk0 <- gamma_k * diag(d)
  Sk0 <- (1 / sqrt(gamma_k)) * diag(d)
  Gk0 <- (1 / gamma_k) * diag(d)
  r <- 0.5
  
  q <- gk
  k <- length(s_list)
  
  alpha <- numeric(k)
  
  # First loop (alpha)
  for (i in k:min(k - M + 1, 0)) {
    si <- s_list[[i]]
    yi <- y_list[[i]]
    
    if(t(si)%*%yi < r*t(si)%*%Bk0%*%si){
      thetai <- ((1-r)*t(si)%*%si)/(t(si)%*%Bk0%*%si - t(si)%*%Bk0)
    }
    else{
      thetai <- 1
    }
    
    yibar <- yi*thetai + (1-yi)*(1-thetai)
    
    
    alpha[i] <- t(si) %*% q / (t(yibar)%*%si)
    q <- q - alpha[i] * yibar
  }
  
  # Middle section: T matrix & a vector
  a <- numeric(k)
  if (k > 0) {
    a[1] <- sum(Bk0 %*% s_list[[max(1, k - M + 1)]])
    T <- matrix(0, nrow = k, ncol = k)
    for (j in 1:min(k, M)) {
      T[1, j] <- sum(Bk0 %*% s_list[[k - M + j]])
    }
    
    for (i in 2:k) {
      for (j in i:k) {
        yi_prev <- y_list[[i - 1]]
        si_prev <- s_list[[i - 1]]
        
        if(t(si_prev)%*%yi_prev < r*t(si_prev)%*%Bk0%*%si_prev){
          thetai_prev <- ((1-r)*t(si_prev)%*%si_prev)/(t(si_prev)%*%Bk0%*%si_prev - t(si_prev)%*%Bk0)
        }
        else{
          thetai_prev <- 1
        }
        yi_prev_bar <- yi_prev*thetai_prev + (1-yi)*(1-thetai_prev)
        
        sj <- s_list[[j]]
        ai_prev <- a[i-1]
        T[i, j] <- T[i - 1, j] +
          ((t(yi_prev_bar) %*% sj) / (t(si_prev) %*% yi_prev_bar))* yi_prev[j] -
          ((t(ai_prev) %*% sj) / (t(si_prev) %*% ai_prev)) * ai_prev[j]
      }
      a[i] <- T[i, i]
    }
    
  }
  
  xi <- Gk0 %*% q
  eta <- Sk0 %*% zk
  
  # Second loop
  for (i in (min(k-M+1):k)) {
    si <- s_list[[i]]
    yi <- y_list[[i]]
    if(t(si)%*%yi < r*t(si)%*%Bk0%*%si){
      thetai <- ((1-r)*t(si)%*%si)/(t(si)%*%Bk0%*%si - t(si)%*%Bk0)
    }
    else{
      thetai <- 1
    }
    
    yibar <- yi*thetai + (1-yi)*(1-thetai)
    p <- si/(t(si)%*%yibar)
    beta_i <- (t(yi_bar) %*% p) / (t(yi_bar) %*% si)
    xi <- xi + (alpha[i] - beta_i) * si
    
    term1 <- (t(yibar) %*% eta) / (t(si)%*%yi_bar)
    norm_term <- sqrt(t(si)%*%yibar) * sqrt(t(a[i]) %*% si)
    eta <- eta - term1 * si - ((t(a[i]) %*% eta )/ norm_term) * si
  }
  
  return(list(xi = xi, eta = eta))
}


optimize_with_lbfgs_and_pruning <- function(beta_init, M, p, tau, num_iterations, grad_fn, prune = TRUE) {
  beta <- beta_init
  d <- length(beta)
  Gk <- diag(d)
  Sk <- diag(d)
  omega <- rep(0.9, num_iterations + 1)  # You can customize omega schedule
  
  for (k in 1:num_iterations) {
    gk <- grad_fn(beta)
    zk <- rnorm(d)
    
    # Step 3: L-BFGS approximation
    xi_eta <- estimate_xi_eta(gk = gk, zk = zk, M = M, gamma_k = 1.0)
    xi_tilde <- xi_eta$xi
    eta_tilde <- xi_eta$eta
    
    # Step 4: Moving average approximations
    Gk_gk <- (1 - omega[k]) * (Gk %*% gk) + omega[k + 1] * xi_tilde
    Sk_zk <- (1 - omega[k]) * (Sk %*% zk) + omega[k + 1] * eta_tilde
    
    # Step 5: Normalize
    xi_k <- Gk_gk / sqrt(sum(Gk_gk^2))
    eta_k <- Sk_zk / sqrt(sum(Sk_zk^2))
    
    # Step 6: Parameter update
    epsilon_k <- 1.0 # You can define a custom step-size scheduler
    beta <- beta + epsilon_k * xi_k + sqrt(2 * epsilon_k / tau) * eta_k
    
    # Step 7–8: Pruning
    if (prune) {
      cutoff <- quantile(abs(beta), probs = p)
      beta[abs(beta) < cutoff] <- 0
      # You could increase `p` slightly if sparse rate increase is desired
    }
  }
  
  return(beta)
}




