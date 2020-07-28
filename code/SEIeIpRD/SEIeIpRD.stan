functions{
  
  // SEIR system of ODEs
  real[] SEIR(real t,       // time
              real[] state, // state (S, E, I_early, I_prog, N_IeIp, N_IpD)
              real[] theta, // parameters (beta, gamma, nu, mu, eta)
              data real[] x_r,   // data (real) 
              data int[] x_i) {  // data (integer)
      
      // elementary transitions 
      real dN_SE   = 
        theta[1] * state[1] * (state[3] + 0.8 * state[4]); // beta*S*(I_e+phi_s*I_s + phi_c*I_c)
      real dN_EIe  = theta[2] * state[2]; // gamma * E
      real dN_IeIp = theta[3] * state[3]; // nu * I_e
      real dN_IpR  = theta[4] * state[4]; // mu * I_p
      real dN_IpD  = theta[5] * state[4]; // eta * I_p
      
      // changes in compartment counts
      real dNdt[6];

      // changes in S, E, Ie, Ip, N_IeIp, N_IpD
      dNdt[1] = -dN_SE;         // dS = -beta * S * (Ie + phi*Ip)
      dNdt[2] = dN_SE - dN_EIe; // dE = beta * S * (Ie + phi*Ip) - gamma * E
      dNdt[3] = dN_EIe - dN_IeIp; // dIe = gamma * E - nu * Ie
      dNdt[4] = dN_IeIp - dN_IpR - dN_IpD; // dIp = nu * Ie - (mu + eta) Ip
      dNdt[5] = dN_IeIp; // incidence of progressed infections from early stage
      dNdt[6] = dN_IpD;  // incidence of deaths from progressed infections

      return dNdt;
  }

  // solve the SEIR ODEs, return the path
  // phi: common to all units
  // theta: vary from unit to unit
  // x_r: real valued data, e.g., observation times
  // x_i: integer valued data, e.g., n of timepoints, n cases, n tests
  vector seir_curves(vector phi,   
                     vector theta,
                     data real[] x_r,
                     data int[] x_i) {
                       
    // dimension of x_i                   
    int i_elem = size(x_i);
    
    // number of observation times
    int n_ode_times = x_i[i_elem];   // number of ODE times
    int n_obs_times = n_ode_times-1; // number of observations
    
    // population size
    real popsize_r = x_i[i_elem-1];
    real log_popsize = log(popsize_r);
    
    // array of parameters to ODEs = (beta, gamma, nu, mu, eta)
    real params[5] = to_array_1d(phi[1:5]); 
    
    // container for epi curves
    real epi_curves[n_obs_times, 10];
    
    // container for S, E, Ie, Ip, N_IeIc, N_IcD
    real init_state[6];
    
    // grab the emission distribution parameters
    real alpha_incid_0 = phi[6];
    real alpha_incid_1 = phi[7];
    real kappa_incid   = phi[8];
    real log_rho_death = phi[9];
    real phi_death     = phi[10];
    
    // grab initial state
    init_state[1:4] = to_array_1d(phi[11:14]);
    init_state[5:6] = rep_array(0.0, 2);
    
    // integrate
    epi_curves[, 1:6] = 
        integrate_ode_bdf(
              SEIR, 
              init_state, 
              x_r[1], 
              x_r[2:n_ode_times], 
              params, 
              x_r, x_i,
              1e-10, 1e-10, 1e6);
              
    // Calculate the emission parameters for incident cases and deaths
    for(t in 1:n_obs_times) {
      // incidence
      epi_curves[t, 7] = 
        t == 1 ? 
        inv_logit(alpha_incid_0 + alpha_incid_1 * logit(epi_curves[t, 5] + 1e-8)) :
        inv_logit(alpha_incid_0 + alpha_incid_1 * logit(epi_curves[t, 5] - epi_curves[t-1, 5] + 1e-8)); 
          
      // deaths
      epi_curves[t, 8] = 
        t == 1 ?
        log_rho_death + log_popsize + log(epi_curves[t, 6] + 1e-8) :
        log_rho_death + log_popsize + log(epi_curves[t, 6] - epi_curves[t-1, 6] + 1e-8);
        
      // overdispersion parameters
      epi_curves[t,9] = kappa_incid;
      epi_curves[t,10] = phi_death;
    }
    
    // return epi curves in row major order
    return to_vector(to_array_1d(epi_curves));
  }
  
  // return the data log-likelihoods
  // phi: common to all units
  // theta: vary from unit to unit
  // x_r: ODE times and test rates
  // x_i: integer valued data, e.g., start and end indices, numbers of tests
  vector seir_loglik(vector phi,   
                     vector theta,
                     data real[] x_r,
                     data int[] x_i) {
              
    // vector for log-likelihood
    vector[1] loglik; 
    
    // dimension of x_i                   
    int i_elem = size(x_i);
    
    // number of observation times
    int n_ode_times = x_i[i_elem];   // number of ODE times
    int n_obs_times = n_ode_times-1; // number of observations
    
    // population size
    real popsize_r = x_i[i_elem - 1];
    real log_popsize = log(popsize_r);
    
    // container for S, E, Ie, Ip, N_IeIp, N_IpD
    real init_state[6];
    
    // array of parameters to ODEs = (beta, gamma, nu, mu, eta)
    real params[5] = to_array_1d(phi[1:5]); 
    
    // container for epi curves
    real epi_curves[n_obs_times, 6];
    vector[n_obs_times] rho_incid;      
    vector[n_obs_times] log_mu_death; 
    
    // grab the emission distribution parameters
    real alpha_incid_0 = phi[6];
    real alpha_incid_1 = phi[7];
    real kappa_incid   = phi[8];
    real log_rho_death = phi[9];
    real phi_death     = phi[10];
    
    // grab initial state
    init_state[1:4] = to_array_1d(phi[11:14]);
    init_state[5:6] = rep_array(0.0, 2);
    
    // integrate - also has a function signature
    epi_curves = 
        integrate_ode_bdf(
              SEIR, 
              init_state, 
              x_r[1], 
              x_r[2:n_ode_times], 
              params, 
              x_r, x_i,
              1e-10, 1e-10, 1e6);
              
    // Calculate the adjusted binomial probs for incident cases and deaths
    for(t in 1:n_obs_times) {
      
      // incidence
      rho_incid[t] = 
        t == 1 ? 
        inv_logit(alpha_incid_0 + alpha_incid_1 * logit(epi_curves[t, 5] + 1e-8)) : 
        inv_logit(alpha_incid_0 + alpha_incid_1 * logit(epi_curves[t, 5] - epi_curves[t-1, 5] + 1e-8)); 
          
      // deaths
      log_mu_death[t] = 
        t == 1 ? 
        log_rho_death + log_popsize + log(epi_curves[t, 6] + 1e-8) : 
        log_rho_death + log_popsize + log(epi_curves[t, 6] - epi_curves[t-1, 6] + 1e-8);
    }
    
    // n_times number of case counts in x_i starting at slot 2
    // n_times number of tests in x_i starting at slot 2 + n_obs_times
    // n_times number of deaths in x_i starting at (2 + 2*n_obstimes)
    loglik[1] = 
      (beta_binomial_lpmf(
        x_i[1:n_obs_times] |
          x_i[(1 + n_obs_times):(2 * n_obs_times)],
          kappa_incid * rho_incid,
          kappa_incid * (1 - rho_incid)) +
       neg_binomial_2_log_lpmf(
        x_i[(1 + 2 * n_obs_times):(3 * n_obs_times)] | 
          log_mu_death, phi_death));
                  
    return loglik;
  }
  
  // convert initial fractions to concentrations
  // i.e., frac_carrs, frac_carrs_infec, frac_infec_prog -> S, E, I_m, I_p
  vector init_fracs2concs(vector init_fracs) {
    
    // vector of initial concentrations
    vector[4] init_concs; 
    
    // S, E, Im, Is
    init_concs[1] = 1- init_fracs[1]; 
    init_concs[2] = exp(log(init_fracs[1]) + log1m(init_fracs[2]));
    init_concs[3] = 
      exp(log(init_fracs[1]) + log(init_fracs[2]) + log1m(init_fracs[3])); 
    init_concs[4] = 
      exp(log(init_fracs[1]) + log(init_fracs[2]) + log(init_fracs[3])); 
      
    // return concentrations
    return init_concs; 
  }
}

data {
  int<lower=0, upper=1> priors_only;

  // dimensions of inference objects
  int n_states; // number of states
  int n_obs;    // number of observations
  int num_per;  // (max) number of time points per unit
  
  // dimensions of posterior predictive objects
  int n_obs_pp;   // number of observation times
  int num_per_pp; // max number of time points per unit
  
  int tests[n_obs_pp];    // test counts
  // int popsizes[n_obs_pp]; // state population sizes 

  // integer and real data for map_rect - ODE times includes t0
  // x_i = {case counts, test counts, deaths, popsize, n_ode_times}
  // x_r = {ODE times}
  int x_i[n_states, 2 + 3 * (num_per - 1)];
  real x_r[n_states, num_per];
  
  int x_i_pp[n_states, 2 + 3 * (num_per_pp - 1)];
  real x_r_pp[n_states, num_per_pp];
  
  real log_R0_normal[2];
  real latent_dur_lognormal[2];
  real early_dur_lognormal[2];
  real prog_dur_lognormal[2];
  real IFR_beta[2];
  real frac_carrs_beta[2];
  real frac_carrs_infec_beta[2];
  real frac_infec_early_beta[2];
  real alpha_incid_0_normal[2];
  real alpha_incid_1_beta[2];
  real sqrt_kappa_inv_incid_exponential;
  real rho_death_beta[2];
  real sqrt_phi_inv_death_exponential;
}

transformed data {
  
  // relative infectiousness in clinical and subclinical stages
  // don't forget to change in ODEs if changing here
  real rel_infec_prog = 0.8;
  
  // incides for emission distribution parameters
  int phi_death_inds[n_obs_pp];
  int kappa_incid_inds[n_obs_pp];
  int log_mu_death_inds[n_obs_pp];
  int rho_incid_inds[n_obs_pp];
  
  for(i in 1:n_obs_pp) {
    phi_death_inds[i]    = i * 10;
    kappa_incid_inds[i]  = i * 10 - 1;
    log_mu_death_inds[i] = i * 10 - 2;
    rho_incid_inds[i]    = i * 10 - 3;
  }
}

parameters {

  // parameters controlling the dynamics
  real log_R0;               // R0 = beta / infec_dur
  real<lower=0> latent_dur;  // 1 / gamma
  real<lower=0> early_dur;   // 1 / nu
  real<lower=0> prog_dur;    // 1 / mu
  real<lower=0,upper=1> IFR; // infection fatality rate
  
  // emission parameters for incidence and death
  real<lower=0> alpha_incid_0;         // test prob higher than incid
  real<lower=0,upper=1> alpha_incid_1; // slope parameter for incid
  real<lower=0> sqrt_kappa_inv_incid;  // 1/sqrt(alpha+beta)
  
  real<lower=0,upper=1> rho_death;    // deaths underreported
  real<lower=0> sqrt_phi_inv_death;   // 1/sqrt(phi)

  // initial state concentrations
  real<lower=0,upper=1> frac_carrs;        // initial prevalence
  real<lower=0,upper=1> frac_carrs_infec;  // fraction of carriers that are infectious
  real<lower=0,upper=1> frac_infec_early;  // fraction of infecteds that are in the early stage
}

transformed parameters {

  // parameters common to all units
  vector[14] phi;
  
  // theta contains parameters that vary among states
  vector[0] theta[n_states];
  
  // rates for transitions from I_severe to recovered and dead
  real prog2rec_rate   = (1 - IFR) / prog_dur;
  real prog2death_rate = IFR / prog_dur;
  
  // log expected infectious period duration 
  // (accounting for relative infectiousness of later state infection)
  real log_infec_dur = log(early_dur + rel_infec_prog * prog_dur);
  
  // initial compartment fractions, susc, carrs_infec, infec_early, prog_clin
  vector[3] init_fracs = 
    [frac_carrs, frac_carrs_infec, frac_infec_early]'; 
  
  // rates of state transition
  phi[1] = exp(log_R0 - log_infec_dur); // beta
  phi[2] = 1 / latent_dur;  // E  -> Ie, gamma
  phi[3] = 1 / early_dur;   // Ie -> Ip, nu
  phi[4] = prog2rec_rate;   // Ip -> R,  mu
  phi[5] = prog2death_rate; // Ip -> D,  eta
  
  // emission distribution parameters
  phi[6] = alpha_incid_0;
  phi[7] = alpha_incid_1;
  phi[8] = inv_square(sqrt_kappa_inv_incid);
  phi[9] = log(rho_death);
  phi[10] = inv_square(sqrt_phi_inv_death); 
  
  // initial compartment fractions
  phi[11:14] = init_fracs2concs(init_fracs);
}

model {
  
  // priors
  log_R0     ~ normal(log_R0_normal[1], log_R0_normal[2]);
  latent_dur ~ lognormal(latent_dur_lognormal[1], latent_dur_lognormal[2]);
  early_dur  ~ lognormal(early_dur_lognormal[1], early_dur_lognormal[2]);
  prog_dur   ~ lognormal(prog_dur_lognormal[1], prog_dur_lognormal[2]);
  
  // fraction progressing to clinical infection and death (+jacobian)
  IFR ~ beta(IFR_beta[1], IFR_beta[2]); 
  
  // initial concentrations + jacobian adjustment
  frac_carrs       ~ beta(frac_carrs_beta[1], frac_carrs_beta[2]);
  frac_carrs_infec ~ beta(frac_carrs_infec_beta[1], frac_carrs_infec_beta[2]);
  frac_infec_early ~ beta(frac_infec_early_beta[1], frac_infec_early_beta[2]); 
  
  // underreporting params
  alpha_incid_0        ~ normal(alpha_incid_0_normal[1], alpha_incid_0_normal[2]) T[0, ];
  alpha_incid_1        ~ beta(alpha_incid_1_beta[1], alpha_incid_1_beta[2]); 
  sqrt_kappa_inv_incid ~ exponential(sqrt_kappa_inv_incid_exponential);
  
  rho_death ~ beta(rho_death_beta[1], rho_death_beta[2]); 
  sqrt_phi_inv_death ~ exponential(sqrt_phi_inv_death_exponential);
  
  // likelihood
  // vector map_rect((vector, vector, real[], int[]):vector f,
  //                  vector phi, vector[] thetas,
  //                  data real[ , ] x_rs, data int[ , ] x_is);
  // in each iteration of map_rect; 
  // vector[k] colvec; a column vector
  // row_vector[k] rowvec; a row vector
  // theta[i] is a vector of length n_elem1
  // x_rs[i] is an array of length n_elem2
  // vector res = f(phi, theta[i], x_rs[i], x_is[i]);
   if (!priors_only) {
  target += sum(map_rect(seir_loglik, phi, theta, x_r, x_i)); 
   }
}

generated quantities {
  
  // epi curves (S, E, I, N_EI, Y_pred)
  real cases_pp[n_obs_pp];
  real deaths_pp[n_obs_pp]; 
  vector[10 * n_obs_pp] epi_curves;
  
  // solve for epi curves
  epi_curves = map_rect(seir_curves, phi, theta, x_r_pp, x_i_pp);
  
  // simulate from the posterior predictive
  for(j in 1:n_obs_pp) {
    
    // incidence
    cases_pp[j] = 
      beta_binomial_rng(
        tests[j], 
        epi_curves[rho_incid_inds[j]] * epi_curves[kappa_incid_inds[j]], 
        (1 - epi_curves[rho_incid_inds[j]]) * epi_curves[kappa_incid_inds[j]]);
        
    // deaths
    deaths_pp[j] = 
      neg_binomial_2_log_rng(
        epi_curves[log_mu_death_inds[j]],
        epi_curves[phi_death_inds[j]]);
  }
}

