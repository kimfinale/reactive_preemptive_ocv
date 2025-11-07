
data {
  int<lower=1> N;
  int<lower=1> G;
  int<lower=1,upper=G> g[N];
  vector<lower=0>[N] t_pk;
  int<lower=1> I_pk[N];
  vector<lower=0>[N] Dur;
}
parameters {
  real mu_t0;         real mu_I0;          real mu_D0;
  real<lower=0> sigma_mu_t;   real<lower=0> sigma_mu_I;  real<lower=0> sigma_mu_D;
  vector[G] mu_t_raw; vector[G] mu_I_raw;  vector[G] mu_D_raw;

  real<lower=0> sigma_t;
  real<lower=0> phi_I;
  real<lower=0> shape_D;        real<lower=0> rate_D;

  cholesky_factor_corr[3] L_corr;
}
transformed parameters {
  vector[G] mu_t = mu_t0 + sigma_mu_t * mu_t_raw;
  vector[G] mu_I = mu_I0 + sigma_mu_I * mu_I_raw;
  vector[G] mu_D = mu_D0 + sigma_mu_D * mu_D_raw;
}
model {
  // Hyper-priors
  mu_t0 ~ normal(log(5), 1.5);
  mu_I0 ~ normal(log(20), 2);
  mu_D0 ~ normal(log(15), 1.5);
  sigma_mu_t ~ exponential(1);
  sigma_mu_I ~ exponential(1);
  sigma_mu_D ~ exponential(1);
  mu_t_raw ~ normal(0,1);
  mu_I_raw ~ normal(0,1);
  mu_D_raw ~ normal(0,1);

  sigma_t  ~ exponential(1);
  phi_I    ~ exponential(1);
  shape_D  ~ gamma(2,0.1);
  rate_D   ~ gamma(2,0.1);
  L_corr   ~ lkj_corr_cholesky(2);

  // Likelihood with Gaussian copula
  for (n in 1:N) {
    int g_n = g[n];
    real mu_tn = mu_t[g_n];
    real mu_In = exp(mu_I[g_n]);
    real mu_Dn = exp(mu_D[g_n]);
    
    // Dur (duration) must be larger than t_pk (peak week)
    if (Dur[n] >= t_pk[n]) { 
      // Marginal log-densities
      target += lognormal_lpdf(t_pk[n] | mu_tn, sigma_t);
      target += neg_binomial_2_lpmf(I_pk[n] | mu_In, 1/phi_I);
      target += gamma_lpdf(Dur[n] | shape_D, rate_D);
  
      // Copula part
      real u1 = lognormal_cdf(t_pk[n], mu_tn, sigma_t);
      real F_lo = neg_binomial_2_cdf(I_pk[n]-1, mu_In, 1/phi_I);
      real F_hi = neg_binomial_2_cdf(I_pk[n], mu_In, 1/phi_I);
      real u2 = 0.5*(F_lo + F_hi);             // mid-CDF for discrete NB
      real u3 = gamma_cdf(Dur[n], shape_D, rate_D);
  
      vector[3] z;
      z[1] = inv_Phi(u1);
      z[2] = inv_Phi(u2);
      z[3] = inv_Phi(u3);
  
      target += multi_normal_cholesky_lpdf(z | rep_vector(0,3), L_corr)
                - normal_lpdf(z[1] | 0,1) - normal_lpdf(z[2] | 0,1)
                - normal_lpdf(z[3] | 0,1);
    } else {
      // Assign zero probability to illegal cases
      target += negative_infinity();
    }
  }
}


