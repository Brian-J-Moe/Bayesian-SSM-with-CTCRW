/*=============================================================================
 * Bayesian CTCRW State-Space Model with a: 
 * linear temperature restorative force
 * asymptotic depth preference force as a function of age
 *=============================================================================*/


functions {
  
  real filter_segment_1d(int idx_start, int n_obs, 
                         array[] real age_t,
                         array[] real obs_1d,
                         vector dt_vec, vector hpe_vec,
                         real beta_s, real sigma_s,
                         vector temp_z, vector temp_grad_z,
                         vector depth_z, vector depth_grad_z,
                         real T_opt_z, real alpha, 
                         real d_inf_z, real d_0, 
                         real kappa_d, real omega) {
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
      
      real f_temp = alpha * (T_opt_z - temp_z[idx]) * temp_grad_z[idx];
      
      real D_pref = d_0 + (d_inf_z - d_0) * (1 - exp(-kappa_d * age_t[idx]));
      real f_depth = omega * (D_pref - depth_z[idx]) * depth_grad_z[idx];
      
      real f = f_temp + f_depth;
      
      // ---- Process noise ----
      // Q2 = [[q11, q12], [q12, q22]]
      real q11;
      real q12;
      real q22;
      
      real a;
      real b;
      real c;
      
      real dti = dt_vec[idx];
      real bdt = beta_s * dti;
      real dt2 = square(dti);
      real dt3 = dt2 * dti;
      real bdt2 = square(bdt);
      real bdt3 = bdt * bdt2;
      
      if (bdt < 0.01) {
        
        a = dti * (1 - bdt / 2.0 + bdt2 / 6.0 - bdt3 / 24);
        b = 1.0 - bdt + bdt2 / 2.0 - bdt3 / 6.0;
        c = dt2 * (0.5 - bdt / 6.0 + bdt2 / 24.0 - bdt3 / 120.0);

        q11 = sigma2 * dt3 * (1.0/3.0 - bdt/4.0 + 7.0 * bdt2/60.0 - bdt3/24.0);
        q12 = sigma2 * dt2 * (0.5     - bdt/2.0 + 7.0 * bdt2/24.0 - bdt3/8.0);
        q22 = sigma2 * dti * (1.0     - bdt     + 2.0 * bdt2/3.0  - bdt3/3.0);
        
      } else {
        real em1 = expm1(-bdt);      
        real em2 = expm1(-2.0 * bdt);

        a = -em1 / beta_s;  
        b = exp(-bdt);
        c = (bdt + em1) / beta2;  
        
        q12 = 0.5 * sigma2 * square(a);
        q22 = -0.5 * sigma2 * em2 / beta_s;
        q11 = sigma2 * dt3 * (bdt + 2.0 * em1 - 0.5 * em2) / bdt3;
      }
      
      // ---- Predict state ----
      // z_pred = T2 * z_filt = [z1 + a*z2, b*z2]
      real zp1 = z1 + a * z2 + c * f;
      real zp2 = b * z2 + a * f;
      
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
   *---------------------------------------------------------------------------*/
  real partial_sum_lpmf(array[] int seg_slice, int start, int end,
                        array[] int seg_start, array[] int seg_length, 
                        array[] int seg_id, array[] real age_t,
                        array[] real obs_x, array[] real obs_y,
                        vector dt_vec, vector hpe_vec,
                        array[] real beta, array[] real sigma, 
                        vector temp_z, vector temp_dX_z, vector temp_dY_z,
                        vector depth_z, vector depth_dX_z, vector depth_dY_z,
                        real T_opt_z, real alpha, real d_inf_z, 
                        real d_0, real kappa_d, real omega) {
    real ll = 0;
    
    for (k in 1:size(seg_slice)) {
      int s = seg_slice[k];
      int ind = seg_id[s];
      int idx_s = seg_start[s];
      int n_s = seg_length[s];
     
      
      // Filter 
      ll += filter_segment_1d(idx_s, n_s, age_t, obs_x, 
                              dt_vec, hpe_vec,
                              beta[ind], sigma[ind], 
                              temp_z, temp_dX_z,
                              depth_z, depth_dX_z,
                              T_opt_z, alpha, d_inf_z, 
                              d_0, kappa_d, omega);
      ll += filter_segment_1d(idx_s, n_s, age_t, obs_y, 
                              dt_vec, hpe_vec,
                              beta[ind], sigma[ind], 
                              temp_z, temp_dY_z,
                              depth_z, depth_dY_z,
                              T_opt_z, alpha, d_inf_z, 
                              d_0, kappa_d, omega);
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
  array[N_total] real age_t;
  array[N_total] real obs_x;           // x positions (standardized)
  array[N_total] real obs_y;           // y positions (standardized)
  vector<lower=0>[N_total] dt;         // time intervals (minutes)
  vector<lower=0>[N_total] HPEm;       // predicted horizontal position error
  
  // Observed annual mean and sd temperatures from receivers
  real temp_mean;
  real temp_sd;
  
  // Position based temp gradients
  vector[N_total] temp_dX;
  vector[N_total] temp_dY;
  vector[N_total] temp;              // estimated temperature for each position
  
  
  // Observed mean and sd of depths (from the raster)
  real depth_mean;
  real depth_sd;
  
  // Position based depth gradients
  vector[N_total] depth_dX;
  vector[N_total] depth_dY;
  vector[N_total] depth;
  
  real depth_max; // prior for d_infty
  real depth_yoy; // yoy preferred depth
  
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
  
  // z-score temp data
  vector[N_total] temp_z = (temp - temp_mean) / temp_sd;
  
  vector[2 * N_total] combined_temp_grad;
  combined_temp_grad = append_row(temp_dX, temp_dY);
  
  real temp_grad_sd = sd(combined_temp_grad);
  vector[N_total] temp_dX_z = temp_dX / temp_grad_sd;
  vector[N_total] temp_dY_z = temp_dY / temp_grad_sd;
  
  
  // z-score depth data
  vector[N_total] depth_z = (depth - depth_mean) / depth_sd;
  
  vector[2 * N_total] combined_depth_grad;
  combined_depth_grad = append_row(depth_dX, depth_dY);
  
  real depth_grad_sd = sd(combined_depth_grad);
  vector[N_total] depth_dX_z = depth_dX / depth_grad_sd;
  vector[N_total] depth_dY_z = depth_dY / depth_grad_sd;
  
  real depth_max_z = (depth_max - depth_mean) / depth_sd;
  real d_0 = (depth_yoy - depth_mean) / depth_sd;
  
  // Segment index array for reduce_sum slicing
  array[N_segments] int seg_indices;
  array[N_total] int ind_for_obs;
  
  for (s in 1:N_segments) seg_indices[s] = s;
  
  for (s in 1:N_segments) {
    int id = seg_id[s];
    int start = seg_start[s];
    int n = seg_length[s];
    
    for (i in 0:(n-1)) ind_for_obs[start + i] = id;
    
  }
  
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
  
  real<lower=-2, upper=4> T_opt_z;  // thermal optimum for restoring drift
  
  real d_inf_z; // asymptotic depth (i.e. depths prefered by largest individuals)
  real log_kappa_d;   // rate of depth transition
  
  real log_alpha; // strength of thermoregulatory response
  real log_omega; // strength of ontogenetic depth shift response          
  
}


/*------------------------- Transformed Parameters --------------------------*/

transformed parameters {
  real T_opt = (T_opt_z * temp_sd) + temp_mean;
  
  real d_inf = (d_inf_z * depth_sd) + depth_mean;
  
  
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> omega = exp(log_omega);
  real<lower=0> kappa_d = exp(log_kappa_d);
  
  array[N_indiv] real<lower=0> beta;
  array[N_indiv] real<lower=0> sigma;
  
  real<lower=0> mu_beta  = exp(log_mu_beta);
  real<lower=0> mu_sigma = exp(log_mu_sigma);
  
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
  
  // Temperature parameters
  T_opt_z ~ normal(0, 1.5);
  log_alpha ~ normal(0, 1.5);
  
  // Depth parameters
  d_inf_z ~ normal(depth_max_z, 0.5);
  log_kappa_d ~ normal(log(0.5), 0.5);
  log_omega ~ normal(0, 1);

  // Likelihood via parallelized Kalman filter
  target += reduce_sum(partial_sum_lpmf, 
                       seg_indices, grainsize,
                       seg_start, seg_length, 
                       seg_id, age_t, obs_x, obs_y, 
                       dt, HPEm, beta, sigma, 
                       temp_z, temp_dX_z, temp_dY_z,
                       depth_z, depth_dX_z, depth_dY_z,
                       T_opt_z, alpha, d_inf_z, 
                       d_0, kappa_d, omega);
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
      
      log_lik[ind] += filter_segment_1d(idx_s, n_s, age_t, 
                                        obs_x, dt, HPEm,
                                        beta[ind], sigma[ind],
                                        temp_z,  temp_dX_z, 
                                        depth_z, depth_dX_z,
                                        T_opt_z, alpha,
                                        d_inf_z, d_0,
                                        kappa_d, omega);
                                        
      log_lik[ind] += filter_segment_1d(idx_s, n_s, age_t, 
                                        obs_y, dt, HPEm,
                                        beta[ind], sigma[ind],
                                        temp_z,  temp_dY_z, 
                                        depth_z, depth_dY_z,
                                        T_opt_z, alpha,
                                        d_inf_z, d_0,
                                        kappa_d, omega);
    }
  }
}
