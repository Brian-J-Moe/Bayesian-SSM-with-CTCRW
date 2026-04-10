/*=============================================================================
 * Bayesian CTCRW State-Space Model
 *=============================================================================*/


functions {

  real filter_segment_1d(int idx_start, int n_obs, 
                         array[] real obs_1d,
                         vector dt_vec, vector hpe_vec,
                         real beta_s, real sigma_s) {
    real ll = 0;
    
    // State: [position, velocity]
    real z1 = obs_1d[idx_start];
    real z2 = 0.0;
    
    // 2x2 covariance as 3 scalars (symmetric): [[p11, p12], [p12, p22]]
    real p11 = 4.0;
    real p12 = 0.0;
    real p22 = 0.5;
    
    real beta2 = square(beta_s);
    real sigma2 = square(sigma_s);
    
    for (i in 2:n_obs) {
      int idx = idx_start + i - 1;
      real dti = dt_vec[idx];
      real bdt = beta_s * dti;
      
      // ---- Transition parameters ----
      // T2 = [[1, a], [0, b]]
      real a;
      real b;
      
      if (bdt < 0.01) {
        a = dti;
        b = 1.0;
      } else if (bdt > 20) {
        a = 1.0 / beta_s;
        b = 0.0;
      } else {
        real exp_bdt = exp(-bdt);
        a = (1.0 - exp_bdt) / beta_s;
        b = exp_bdt;
      }
      
      // ---- Process noise ----
      // Q2 = [[q11, q12], [q12, q22]]
      real q11;
      real q12;
      real q22;
      
      if (bdt < 0.01) {
        real dt2 = square(dti);
        q11 = sigma2 * dt2 * dti / 3.0;
        q12 = sigma2 * dt2 / 2.0;
        q22 = sigma2 * dti;
      } else if (bdt > 20) {
        q11 = (sigma2 / beta2) * (dti - 1.5 / beta_s);
        q12 = sigma2 / (2.0 * beta2);
        q22 = sigma2 / (2.0 * beta_s);
      } else {
        real exp_bdt = exp(-bdt);
        real exp_2bdt = exp(-2.0 * bdt);
        real inv_beta = 1.0 / beta_s;
        real inv_beta2 = 1.0 / beta2;
        q11 = sigma2 * inv_beta2 * (dti - 2.0 * (1.0 - exp_bdt) * inv_beta 
              + (1.0 - exp_2bdt) * 0.5 * inv_beta);
        q12 = sigma2 * 0.5 * inv_beta2 * (1.0 - 2.0 * exp_bdt + exp_2bdt);
        q22 = sigma2 * 0.5 * inv_beta * (1.0 - exp_2bdt);
      }
      
      // ---- Predict state ----
      // z_pred = T2 * z_filt = [z1 + a*z2, b*z2]
      real zp1 = z1 + a * z2;
      real zp2 = b * z2;
      
      // ---- Predict covariance ----
      // Pp = T2 * P * T2' + Q2  (inlined 2x2 multiplication)
      // T = [[1,a],[0,b]], P = [[p11,p12],[p12,p22]]
      // Pp11 = p11 + 2*a*p12 + a^2*p22 + q11
      // Pp12 = b*(p12 + a*p22) + q12
      // Pp22 = b^2*p22 + q22
      real a2 = square(a);
      real b2 = square(b);
      real pp11 = p11 + 2.0 * a * p12 + a2 * p22 + q11;
      real pp12 = b * (p12 + a * p22) + q12;
      real pp22 = b2 * p22 + q22;
      
      // ---- Innovation ----
      real h2 = square(hpe_vec[idx]);
      real F_val = pp11 + h2;
      real F_inv = 1.0 / F_val;
      real innov = obs_1d[idx] - zp1;
      
      // ---- Kalman gain ----
      // K = Pp * Z' / F = [pp11, pp12]' / F  (Z = [1, 0])
      real k1 = pp11 * F_inv;
      real k2 = pp12 * F_inv;
      
      // ---- Update state ----
      z1 = zp1 + k1 * innov;
      z2 = zp2 + k2 * innov;
      
      // ---- Update covariance ----
      // Simplified from Joseph form for scalar observation:
      //   p11 = Pp11 * h^2 / F
      //   p12 = Pp12 * h^2 / F
      //   p22 = Pp22 - Pp12^2 / F
      real h2_F_inv = h2 * F_inv;
      p11 = pp11 * h2_F_inv;
      p12 = pp12 * h2_F_inv;
      p22 = pp22 - pp12 * pp12 * F_inv;
      
      // ---- Log-likelihood (scalar normal, not multivariate) ----
      ll += normal_lpdf(obs_1d[idx] | zp1, sqrt(F_val));
    }
    
    return ll;
  }
  
  
  /*---------------------------------------------------------------------------
   * partial_sum: Compute log-likelihood for a slice of segments
   *
   * Used by reduce_sum to distribute segments across threads.
   * Each segment is processed independently (both x and y dimensions).
   *---------------------------------------------------------------------------*/
  real partial_sum_lpmf(array[] int seg_slice, int start, int end,
                        array[] int seg_start, array[] int seg_length, 
                        array[] int seg_id,
                        array[] real obs_x, array[] real obs_y,
                        vector dt_vec, vector hpe_vec,
                        array[] real beta, array[] real sigma) {
    real ll = 0;
    
    for (k in 1:size(seg_slice)) {
      int s = seg_slice[k];
      int ind = seg_id[s];
      int idx_s = seg_start[s];
      int n_s = seg_length[s];
      
      // Filter x dimension
      ll += filter_segment_1d(idx_s, n_s, obs_x, dt_vec, hpe_vec,
                              beta[ind], sigma[ind]);
      // Filter y dimension
      ll += filter_segment_1d(idx_s, n_s, obs_y, dt_vec, hpe_vec,
                              beta[ind], sigma[ind]);
    }
    
    return ll;
  }
}


/*--------------------------------- Data ------------------------------------*/

data {
  int<lower=1> N_total;     // total observations
  int<lower=1> N_indiv;     // number of individuals
  int<lower=1> N_segments;  // number of segments
  
  array[N_segments] int<lower=1> seg_start;   // starting index of each segment
  array[N_segments] int<lower=1> seg_length;  // observations per segment
  array[N_segments] int<lower=1> seg_id;      // individual id for each segment
  
  // Observations stored as separate 1D arrays (avoids array[] vector overhead)
  array[N_total] real obs_x;           // x positions (standardized)
  array[N_total] real obs_y;           // y positions (standardized)
  vector<lower=0>[N_total] dt;         // time intervals (minutes)
  vector<lower=0>[N_total] HPEm;       // predicted horizontal position error
  
  // Priors
  real<lower=0> beta_median;
  real<lower=0> beta_sd;
  real<lower=0> sigma_median;
  real<lower=0> sigma_sd;
  real<lower=0> tau_rate;
  
  // Threading control
  int<lower=1> grainsize;
}


/*--------------------------- Transformed Data ------------------------------*/

transformed data {
  // Segment index array for reduce_sum slicing
  array[N_segments] int seg_indices;
  for (s in 1:N_segments) seg_indices[s] = s;
}


/*------------------------------ Parameters ---------------------------------*/

parameters {
  // Population-level hyperparameters (log scale)
  real<lower=log(1e-4), upper=log(1e2)> log_mu_beta;
  real<lower=log(1e-4), upper=log(1e4)> log_mu_sigma;
  real<lower=0, upper=2> tau_beta;
  real<lower=0, upper=2> tau_sigma;
  
  // Individual-level raw deviates (non-centered parameterization)
  array[N_indiv] real beta_raw;
  array[N_indiv] real sigma_raw;
}


/*------------------------- Transformed Parameters --------------------------*/

transformed parameters {
  real<lower=0> mu_beta  = exp(log_mu_beta);
  real<lower=0> mu_sigma = exp(log_mu_sigma);
  
  array[N_indiv] real<lower=0> beta;
  array[N_indiv] real<lower=0> sigma;
  
  for (i in 1:N_indiv) {
    beta[i]  = exp(log_mu_beta  + tau_beta  * beta_raw[i]);
    sigma[i] = exp(log_mu_sigma + tau_sigma * sigma_raw[i]);
  }
}


/*--------------------------------- Model -----------------------------------*/

model {
  // Hyperpriors
  log_mu_beta  ~ normal(log(beta_median),  beta_sd);
  log_mu_sigma ~ normal(log(sigma_median), sigma_sd);
  
  tau_beta  ~ exponential(tau_rate);
  tau_sigma ~ exponential(tau_rate);
  
  // Individual raw deviates
  beta_raw  ~ std_normal();
  sigma_raw ~ std_normal();
  
  // Likelihood via parallelized Kalman filter
  target += reduce_sum(partial_sum_lpmf, seg_indices, grainsize,
                       seg_start, seg_length, seg_id,
                       obs_x, obs_y, dt, HPEm,
                       beta, sigma);
}


/*--------------------------- Generated Quantities --------------------------*/
/* Parameters only - smoother runs in R post-processing */

generated quantities {
  // Individual-level log-likelihoods (for LOO-CV if needed)
  // Computed per-individual by summing over their segments
  array[N_indiv] real log_lik;
  
  {
    // Initialize
    for (i in 1:N_indiv) log_lik[i] = 0;
    
    for (s in 1:N_segments) {
      int ind = seg_id[s];
      int idx_s = seg_start[s];
      int n_s = seg_length[s];
      
      log_lik[ind] += filter_segment_1d(idx_s, n_s, obs_x, dt, HPEm,
                                        beta[ind], sigma[ind]);
      log_lik[ind] += filter_segment_1d(idx_s, n_s, obs_y, dt, HPEm,
                                        beta[ind], sigma[ind]);
    }
  }
}
