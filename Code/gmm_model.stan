data {
  int<lower=1> N;        // Number of data points
  int<lower=1> K;        // Number of components (clusters)
  vector[2] zs[N];       // Observed data (2D points)
}

parameters {
  simplex[K] theta;      // Mixture weights (sum to 1)
  vector[2] mu[K];       // Cluster means
  cov_matrix[2] Sigma[K];// Covariance matrices
}

model {
  // Priors
  for (k in 1:K) {
    mu[k] ~ normal(0, 5);              // Prior on means
    Sigma[k] ~ wishart(3, diag_matrix(rep_vector(1.0, 2))); // Prior on covariance
  }
  
  // Likelihood
  for (n in 1:N) {
    vector[K] log_prob;
    for (k in 1:K) {
      log_prob[k] = log(theta[k]) + multi_normal_lpdf(zs[n] | mu[k], Sigma[k]);
    }
    target += log_sum_exp(log_prob);
  }
}
