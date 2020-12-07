find_initial_conditions <- function(city_name, folder_name, prev_folder_name) {
  oc_city_pop <- read_csv("data/oc_city_pop.csv")

  full_county_popsize <- 3175692
  initial_county_popsize <- read_rds(here("code", model_name, loc_name, folder_name, "model_objects.rds"))$popsize

  full_city_popsize <- filter(oc_city_pop, city == city_name)$population
  city_D_plus_R <- round(filter(oc_city_pop, city == city_name)$prop_cases * (full_county_popsize - initial_county_popsize))
  initial_city_popsize <- full_city_popsize - city_D_plus_R
  folder_name

  model_objects <- read_rds(here("code", model_name, loc_name, folder_name, "model_objects.rds"))

  n <- 2000
  target <-
    tibble(frac_carrs = rbeta(n, model_objects$frac_carrs_beta[1], model_objects$frac_carrs_beta[2]),
           frac_carrs_infec = rbeta(n, model_objects$frac_carrs_infec_beta[1], model_objects$frac_carrs_infec_beta[2]),
           frac_infec_early = rbeta(n, model_objects$frac_infec_early_beta[1], model_objects$frac_infec_early_beta[2]),
           IFR = sample(unname(unlist(rstan::extract(read_rds(here("code", model_name, loc_name, prev_folder_name, "oc_post.rds")), "IFR"))), n, replace = F)) %>%
    mutate(EIeIp = frac_carrs * model_objects$popsize,
           S = (1 - frac_carrs) * model_objects$popsize,
           IeIP = frac_carrs_infec * EIeIp,
           E = (1 - frac_carrs_infec) * EIeIp,
           Ie = frac_infec_early * IeIP,
           Ip = (1 - frac_infec_early) * IeIP) %>%
    mutate(R = (1 - IFR) * (full_county_popsize - initial_county_popsize)) %>%
    select(S, E, Ie, Ip, R) %>%
    `*`(full_city_popsize / full_county_popsize) %>%
    round() %>%
    mutate(., D = full_city_popsize - rowSums(.))

  init_states <- round(colMeans(target))

  C <- optimize(f = function(C) sum(extraDistr::ddirmnom(x = target, size = full_city_popsize, alpha = init_states / C, log = T)), lower = 1, upper = 10000, maximum = T)$maximum

  list(init_states = init_states, C = C)
}

create_city_data <- function(city_name, start_date, end_date, time_interval_in_days){
  oc_city_data <- read_csv("data/oc_city_data.csv")

  lump_oc_data <-
    function(oc_data,
             time_interval_in_days,
             first_day = "0000-01-01",
             last_day = "9999-12-31") {
      oc_data %>%
        filter(date >= lubridate::ymd(first_day),
               date <= lubridate::ymd(last_day)) %>%
        group_by(lump = as.integer(floor((max(date) - date) / time_interval_in_days))) %>%
        filter(n() == time_interval_in_days) %>%
        dplyr::summarize(start_date = min(date),
                         end_date = max(date),
                         cases = sum(cases),
                         tests = sum(tests),
                         deaths = sum(deaths)) %>%
        dplyr::select(-lump) %>%
        arrange(start_date)
    }

  lumped_data <- oc_city_data %>%
    filter(city == city_name) %>%
    rename(date = posted_date) %>%
    rename_all(~str_remove(., "new_")) %>%
    lump_oc_data(first_day = start_date, last_day = end_date, time_interval_in_days = time_interval_in_days)

  dat <- lumped_data %>%
    mutate(time = as.numeric((lumped_data$end_date - min(lumped_data$start_date) + 1L) / 7)) %>%
    select(time, cases, tests, deaths, date = end_date)

  dat
}

fit_city_model <- function(dat, init_states, C) {
  popsize <- sum(init_states)
  obs_times <- dat$time

  strata <- NULL # no strata
  compartments <- c("S", "E", "Ie", "Ip", "R", "D")

  rates <-
    list(rate(rate = "beta * (Ie + 0.8 * Ip)", # individual level rate (unlumped)
              from = "S",        # source compartment
              to   = "E",        # destination compartment
              incidence = F),    # compute incidence of S2I transitions, required for simulating incidence data
         rate(rate = "gamma",
              from = "E",
              to = "Ie",
              incidence = F),
         rate(rate = "nu_early",
              from = "Ie",
              to = "Ip",
              incidence = T),
         rate(rate = "mu_rec",
              from = "Ip",
              to = "R",
              incidence = F),
         rate(rate = "mu_death",       # individual level rate
              from = "Ip",        # source compartment
              to   = "D",        # destination compartment
              incidence = TRUE)) # compute incidence of I2R transitions (not required for simulating data)

  state_initializer <-
    list(stem_initializer(
      init_states = init_states, # must match compartment names
      fixed = F,
      prior = init_states / C,
      dist = "dirmultinom"
    )) # initial state fixed for simulation, we'll change this later

  parameters <- numeric(10); names(parameters) <- c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa")
  constants <- c(t0 = 0)
  tcovar <- data.frame(time = obs_times,
                       tests = dat$tests)
  tmax <- max(tcovar$time)

  # list of emission distribution lists (analogous to rate specification)
  emissions <-
    list(emission(meas_var = "cases", # transition or compartment being measured (S->I transitions)
                  distribution    = "betabinomial",        # emission distribution
                  emission_params =
                    c("tests",
                      "kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)",
                      "kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)"), # distribution pars, here overdispersion and mean
                  incidence       = TRUE,                  # is the data incidence
                  obstimes        = obs_times), # vector of observation times
         emission(meas_var = "deaths",
                  distribution = "negbinomial",
                  emission_params = c("phi_death",
                                      "rho_death * Ip2D"),
                  incidence = T,
                  obstimes        = obs_times)) # vector of observation times)
  # list of emission distribution lists (analogous to rate specification)

  dynamics <-
    stem_dynamics(
      rates = rates,
      parameters = parameters,
      state_initializer = state_initializer,
      compartments = compartments,
      constants = constants,
      tcovar = tcovar,
      tmax = tmax,
      compile_ode = T,   # compile ODE functions
      compile_rates = F, # compile MJP functions for Gillespie simulation
      compile_lna = T,   # compile LNA functions
      messages = F       # don't print messages
    )



  measurement_process <-
    stem_measure(emissions = emissions,
                 dynamics = dynamics,
                 data = dat %>%
                   mutate(t = obs_times) %>%
                   select(t, cases, deaths))

  stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)

  # Build Priors ------------------------------------------------------------
  to_estimation_scale = function(params_nat) {
    c(R0_est = log(params_nat[["beta"]]) + log(popsize) + log(1 / params_nat[["nu_early"]] + 0.8 / (params_nat[["mu_rec"]] + params_nat[["mu_death"]])), # log(R0)
      dur_latent_est = log(params_nat[["gamma"]]), # -log(dur_latent)
      dur_early_est = log(params_nat[["nu_early"]]), # -log(dur_early)
      dur_progress_est = log(params_nat[["mu_rec"]] + params_nat[["mu_death"]]), # -log(dur_progress)
      ifr_est = log(params_nat[["mu_death"]]) - log(params_nat[["mu_rec"]]), # logit(ifr)
      rho_death_est = logit(params_nat[["rho_death"]]), # logit(rho_death)
      phi_death_est = -0.5 * log(params_nat[["phi_death"]]), # -0.5 * log(phi_death)
      alpha0_est = log(params_nat[["alpha0"]]), # log(alpha0)
      alpha1_est = logit(params_nat[["alpha1"]]), # logit(alpha1)
      kappa_est = -0.5 * log(params_nat[["kappa"]])) # -0.5 * log(kappa)
  }


  from_estimation_scale = function(params_est) {
    c(beta = exp(params_est[["R0_est"]] - log(popsize) - log(exp(-params_est[["dur_early_est"]]) + 0.8 * exp(-params_est[["dur_progress_est"]]))),
      gamma = exp(params_est[["dur_latent_est"]]),
      nu_early = exp(params_est[["dur_early_est"]]),
      mu_rec = exp(params_est[["dur_progress_est"]]) / (1 + exp(params_est[["ifr_est"]])),
      mu_death = exp(params_est[["dur_progress_est"]]) / (1 + exp(-params_est[["ifr_est"]])),
      rho_death = expit(params_est[["rho_death_est"]]),
      phi_death = exp(-2 * params_est[["phi_death_est"]]),
      alpha0 = exp(params_est[["alpha0_est"]]),
      alpha1 = expit(params_est[["alpha1_est"]]),
      kappa = exp(-2 * params_est[["kappa_est"]]))
  }


  logprior =
    function(params_est) {
      sum(dnorm(params_est["R0_est"], -0.2554128198465173693599, 0.7, log = TRUE), # log(R0)
          dnorm(-params_est["dur_latent_est"], 0, 0.22, log = TRUE), # -log(dur_latent)
          dnorm(-params_est["dur_early_est"], 0, 0.22, log = TRUE), # -log(dur_early)
          dnorm(-params_est["dur_progress_est"], 0, 0.22, log = TRUE), # -log(dur_progress)
          dbeta(expit(params_est["ifr_est"]), 1.5, 200, log = TRUE) + params_est["ifr_est"] - 2 * log(exp(params_est["ifr_est"]) + 1), # logit(ifr)
          dbeta(expit(params_est["rho_death_est"]), 8, 2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
          dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
          dtnorm(exp(params_est["alpha0_est"]), mean = 4, sd = 2, a = 0, log = TRUE) + params_est["alpha0_est"], # log(alpha0)
          dbeta(expit(params_est["alpha1_est"]), 3, 1, log = TRUE) + params_est["alpha1_est"] - 2 * log(exp(params_est["alpha1_est"]) + 1), # logit(alpha1)
          dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"]) # -0.5 * log(kappa)
    }

  priors <- list(logprior = logprior,
                 to_estimation_scale = to_estimation_scale,
                 from_estimation_scale = from_estimation_scale)


  # Build par_initializer ---------------------------------------------------
  true_pars =
    c(R0       = 1.5,    # basic reproduction number
      dur_latent = 1,
      dur_early   = 1,      # infectious period duration = 2 days
      dur_progress = 1,
      ifr = 0.06,
      rho_death = 0.7,
      phi_death = 2.2,
      alpha0   = 4, # beta-binomial intercept
      alpha1   = 0.8,    # beta-binomial slope
      kappa    = 2.2)

  parameters =
    c(beta = true_pars[["R0"]] / popsize / (true_pars[["dur_early"]] + true_pars[["dur_progress"]]), # R0 = beta * P / mu
      gamma = 1 / true_pars[["dur_latent"]],
      nu_early = 1 / true_pars[["dur_early"]],
      mu_rec = (1 - true_pars[["ifr"]]) / true_pars[["dur_progress"]],
      mu_death = true_pars[["ifr"]] / true_pars[["dur_progress"]],
      rho_death = true_pars[["rho_death"]],
      phi_death = true_pars[["phi_death"]],
      alpha0 = true_pars[["alpha0"]],
      alpha1 = true_pars[["alpha1"]],
      kappa = true_pars[["kappa"]])

  n_params <- length(parameters)
  par_initializer = function() {
    priors$from_estimation_scale(priors$to_estimation_scale(parameters) + rnorm(n_params, 0, 0.1))
  }


  # Build Kernel ------------------------------------------------------------
  # specify the kernel
  mcmc_kern <-
    mcmc_kernel(
      parameter_blocks =
        list(parblock(
          pars_nat = c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa"),
          # par_est should have different names from pars_nat
          pars_est = c("R0_est", "dur_latent_est", "dur_early_est", "dur_progress_est", "ifr_est", "rho_death_est", "phi_death_est", "alpha0_est", "alpha1_est", "kappa_est"),
          priors = priors,
          # alg = "mvnss",
          alg = "mvnmh",
          sigma = diag(0.01, n_params),
          initializer = par_initializer,
          control =
            # mvnss_control(stop_adaptation = 1e2))),
            mvnmh_control(stop_adaptation = 150000,
                          scale_cooling = 0.85,
                          scale_constant = 1,
                          step_size = 0.25))),
      lna_ess_control = lna_control(bracket_update_iter = 5e3,
                                    joint_initdist_update = FALSE))


  # Fit Model ---------------------------------------------------------------

  n_chains <- 4
  thinning_interval <- 100
  iterations <- 350000

  stem_fit_list <- foreach(chain = 1:n_chains,
                           .packages = "stemr",
                           .export = ls()) %dorng% {
                             fit_stem(stem_object = stem_object,
                                      method = "ode", # or "lna"
                                      mcmc_kern = mcmc_kern,
                                      iterations = iterations,
                                      thinning_interval = thinning_interval,
                                      print_progress = 1e3)
                           }


  multi_chain_stem_fit <- list()
  multi_chain_stem_fit$stem_fit_list <- stem_fit_list
  multi_chain_stem_fit$n_iterations <- 2000
  multi_chain_stem_fit$n_chains <- n_chains
  multi_chain_stem_fit$thinning_interval <- thinning_interval

  multi_chain_stem_fit
}


fit_and_save_city_model <- function(city_name, folder_name, prev_folder_name, time_interval_in_days = 3) {

  start_date <- str_sub(folder_name, end = 10)
  end_date <- str_sub(folder_name, start = 12)

  init_states_C <- find_initial_conditions(city_name = city_name, folder_name = folder_name, prev_folder_name = prev_folder_name)

  city_data <- create_city_data(city_name = city_name, start_date = start_date, end_date = end_date, time_interval_in_days)
  multichain_stem_fit <- fit_city_model(dat = city_data, init_states = init_states_C$init_states, C = init_states_C$C)

  fs::dir_create(path = "code", model_name, "oc_cities", folder_name)
  write_rds(multichain_stem_fit, here("code", model_name, "oc_cities", folder_name, str_c(folder_name, "_", str_replace_all(city_name, " ", "-"), ".rds")))
}
