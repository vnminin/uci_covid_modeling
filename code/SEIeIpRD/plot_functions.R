library(tidyverse)
library(tidybayes)

# Prep Cummulative Incidence ----------------------------------------------
prep_cumulative_inc <- function(model_objects, oc_post) {
  
  prepped_epi_curves <- prep_epi_curves(oc_post = oc_post, model_objects = model_objects) %>%
    pivot_wider(names_from = g_text, values_from = epi_curves)
    
  posterior <- oc_post %>% 
    spread_draws(frac_carrs) %>%
    mutate(frac_susc =1 - frac_carrs) %>%
    select(.draw, frac_susc) %>%
    left_join(prepped_epi_curves, by = ".draw") %>% 
    mutate("Cumulative Incidence" = (frac_susc - S) * model_objects$popsize) %>% 
    mutate("Cumulative Death" = D * model_objects$popsize) %>%
    mutate(Prevalence = round((E + IE + IP) * model_objects$popsize)) %>% 
    select(.draw, t, `Cumulative Incidence`, `Cumulative Death`, Prevalence) %>% 
    pivot_longer(cols = c(starts_with("Cumulative"), Prevalence)) %>% 
    select(-.draw) %>% 
    mutate(usage = ifelse(t <= model_objects$actual_last_day, "train", "test"))
  
  prep_epi_curves(oc_post = oc_post, model_objects = model_objects) %>%
    pivot_wider(names_from = g_text,
                values_from = epi_curves,
                values_fn = list(epi_curves = list)) %>% 
    unnest(-t) %>% 
    transmute() %>% 
    ungroup()
  
  observed <- model_objects$lumped_ochca_covid_forecast %>% 
    mutate("Cumulative Incidence" = cumsum(new_cases),
           "Cumulative Death" = cumsum(new_deaths)) %>% 
    select(t = end_date, "Cumulative Incidence", "Cumulative Death") %>% 
    pivot_longer(-t, values_to = "observed") %>% 
    mutate(usage = ifelse(t <= model_objects$actual_last_day, "train", "test"))
  
  list(posterior = posterior, observed = observed)
}

# Prep Deaths -------------------------------------------------------------
prep_death <- function(model_objects, oc_post, ci_width = 0.8) {
  
  posterior <- oc_post %>%
    spread_draws(deaths_pp[i]) %>% 
    mutate(t = model_objects$dates_pp[i]) %>% 
    ungroup() %>% 
    select(deaths = deaths_pp, t)
  
  observed <- model_objects$lumped_ochca_covid %>% 
    mutate(observed = new_deaths) %>% 
    select(t = end_date, observed)
  
  list(posterior = posterior, observed = observed)
}

# Prep Percent Positivity -------------------------------------------------------------
prep_pos<- function(model_objects, oc_post, ci_width = 0.8) {
  
  observed <- model_objects$lumped_ochca_covid %>% 
    mutate(observed = new_cases/new_tests) %>% 
    dplyr::select(t = end_date, observed)
  
  tests <- data.frame(model_objects$dates_pp, model_objects$tests)%>%
    rename("t"="model_objects.dates_pp", "new_tests"="model_objects.tests")
  


  posterior <- oc_post %>%
    spread_draws(cases_pp[i]) %>% 
    mutate(t = model_objects$dates_pp[i]) %>% 
    ungroup() %>% 
    dplyr::select(cases = cases_pp, t) %>%
    #filter(t<=max(observed$t)) %>%
    left_join(tests, by="t") %>%
    mutate(percent_positive = cases / new_tests) %>%
    dplyr::select(t, percent_positive)
  
  list(posterior = posterior, observed = observed)
}


# Prep pp Time Series -----------------------------------------------------
prep_pp <- function(model_objects, oc_post) {
  posterior <- oc_post %>%
    spread_draws(cases_pp[i], deaths_pp[i]) %>% 
    mutate(t = model_objects$dates_pp[i]) %>% 
    ungroup() %>% 
    select(-starts_with("."), -i) %>% 
    rename_all(~str_remove(., "_pp")) %>% 
    left_join(tibble(t = model_objects$dates_pp, tests = model_objects$tests)) %>% 
    mutate(prop_pos = cases / tests) %>% 
    pivot_longer(-t) %>% 
    mutate(usage = ifelse(t <= model_objects$actual_last_day, "train", "test"))
  
  observed <- model_objects$lumped_ochca_covid_forecast %>% 
    select(t = end_date, cases = new_cases, deaths = new_deaths, tests = new_tests) %>% 
    mutate(prop_pos = cases / tests) %>% 
    pivot_longer(-t, values_to = "observed") %>% 
    mutate(usage = ifelse(t <= model_objects$actual_last_day, "train", "test"))
  
  list(posterior = posterior, observed = observed)
}

# prepped_pp <- prep_pp(model_objects, oc_post)
# ci_width <- c(0.5, 0.8, 0.95)
# 
# prepped_pp$posterior %>%
#   group_by(t, name) %>%
#   median_qi(.width = ci_width) %>%
#   ggplot(aes(t, value)) +
#   geom_lineribbon(size = 0.8, color = "steelblue4") +
#   geom_point(data = prepped_pp$observed, aes(t, observed), fill = "black") +
#   facet_wrap(. ~ name, ncol = 2, scales = "free_y") +
#   ggtitle("Observed & predicted counts in 3 day periods",
#           subtitle = str_c("Posterior Median &", str_c(percent(ci_width), collapse = ", "), "credible intervals", sep = " ")) +
#   xlab("Date") +
#   # ylab("Reported deaths due to COVID-19") +
#   # theme(legend.position = c(0.15, 0.8), legend.title = element_text(size=8), legend.text = element_text(size=9), legend.background = element_rect(fill="transparent")) +
#   theme(legend.position = "bottom", legend.title = element_text(size=8), legend.text = element_text(size=9), legend.background = element_rect(fill="transparent")) +
#   scale_fill_brewer(name="Credibility level", labels=str_c(percent(rev(ci_width))), guide = guide_legend(title.position = "top", direction = "horizontal")) +
#   scale_x_date(breaks=c("8 day"), date_labels = "%b %d") +
#   theme(legend.position = "none") +
#   ylab("")

## Plot Epi Curves ---------------------------------------------------------
prep_epi_curves <- function(oc_post, model_objects) {
  oc_post %>%
    spread_draws(epi_curves[i]) %>%
    ungroup() %>% 
    dplyr::select(i, epi_curves, .draw) %>%
    mutate(g = (i - 1) %% model_objects$n_curves + 1) %>%
    mutate(t = model_objects$dates_pp[ceiling(i / model_objects$n_curves)]) %>% 
    mutate(g_text = model_objects$epi_curve_names[g]) %>%
    mutate(g_text = fct_reorder(g_text, g))  %>%
    dplyr::select(-i, -g) %>% 
    group_by(t, g_text)
}


# Prevalence --------------------------------------------------------------

prep_prevalence <- function(oc_post, model_objects) {
  prep_epi_curves(oc_post = oc_post, model_objects = model_objects) %>%
    pivot_wider(names_from = g_text,
                values_from = epi_curves,
                values_fn = list(epi_curves = list)) %>% 
    unnest(-t) %>% 
    transmute(Prevalence = round((E + IE + IP) * model_objects$popsize)) %>% 
    ungroup()
}

# Plot Priors and Posteriors -------------------------------------------------------------
prep_priors_posts<-function(oc_post, oc_prior){
posterior <- oc_post %>% 
  tidy_draws() %>% 
  select(-ends_with("__"), -matches("\\[\\d+\\]"), -starts_with(".")) %>% 
  pivot_longer(cols = everything()) %>% 
  rename(variable = name) %>% 
  mutate(type = "posterior")

prior <- oc_prior %>% 
  tidy_draws() %>% 
  select(-ends_with("__"), -matches("\\[\\d+\\]"), -starts_with(".")) %>% 
  pivot_longer(cols = everything()) %>% 
  rename(variable = name) %>% 
  mutate(type = "prior")

# All vars
allvars <- bind_rows(posterior, prior) %>%
  group_by(variable, type)

return(allvars)
}

# prepped_epi_curves <- prep_epi_curves(oc_post, model_object)
# ci_width <- c(0.5, 0.8, 0.95)
# 
# prepped_epi_curves %>%
#   median_qi(.width = ci_width) %>%
#   ggplot(aes(t, epi_curves)) +
#   geom_lineribbon() +
#   facet_wrap(~ g_text, scales = "free_y") +
#   scale_fill_brewer() +
#   scale_x_date(breaks=c("14 day"), date_labels = "%b %d") +
#   theme(legend.position = "none") +
#   ylab(NULL) +
#   xlab(NULL)

# Prep Deaths CA -------------------------------------------------------------
prep_death_ca <- function(model_objects, ca_post, ci_width = 0.8) {
  
  posterior <- ca_post %>%
    spread_draws(deaths_pp[i]) %>% 
    mutate(t = model_objects$dates_pp[i]) %>% 
    ungroup() %>% 
    select(deaths = deaths_pp, t)
  
  observed <- model_objects$lumped_ca_covid %>% 
    mutate(observed = new_deaths) %>% 
    ungroup()%>%
    select(t = end_date, observed)
  
  list(posterior = posterior, observed = observed)
}

# Prep Percent Positivity CA -------------------------------------------------------------
prep_pos_ca<- function(model_objects, ca_post, ci_width = 0.8) {
  
  
  observed <- model_objects$lumped_ca_covid %>% 
    mutate(observed = new_cases/new_tests) %>% 
    ungroup()%>%
    select(t = end_date, observed)
  
  
  tests<-model_objects$lumped_ca_covid%>%
    ungroup()%>%
    select(t=end_date, new_tests)
  
  posterior <- ca_post %>%
    spread_draws(cases_pp[i]) %>% 
    mutate(t = model_objects$dates_pp[i]) %>% 
    ungroup() %>% 
    select(cases = cases_pp, t)%>%
    filter(t<=max(observed$t))%>%
    left_join(tests, by="t")%>%
    mutate(
      percent_positive=cases/new_tests
    )%>%
    select(t, percent_positive)
  
  
  
  list(posterior = posterior, observed = observed)
}

# Prep Cummulative Incidence CA ----------------------------------------------
prep_cumulative_inc_ca <- function(model_objects, ca_post) {
  
  S_D_curves <- ca_post %>% 
    spread_draws(epi_curves[i]) %>% 
    ungroup() %>% 
    select(i, epi_curves, .draw) %>% 
    mutate(g = (i - 1) %% model_objects$n_curves + 1) %>% 
    mutate(g_text = model_objects$epi_curve_names[g]) %>% 
    mutate(g_text = fct_reorder(g_text, g)) %>% 
    mutate(t = model_objects$dates_pp[ceiling(i / model_objects$n_curves)]) %>% 
    select(-i, - g) %>% 
    pivot_wider(names_from = g_text, values_from = epi_curves) %>%
    select(.draw, t, S, D)
  
  posterior <- ca_post %>% 
    spread_draws(frac_susc) %>% 
    select(.draw, frac_susc) %>%
    left_join(S_D_curves, by = ".draw") %>% 
    mutate("Cumulative Incidence" = (frac_susc - S) * model_objects$popsize) %>% 
    mutate("Cumulative Death" = D * model_objects$popsize) %>%
    select(.draw, t, starts_with("Cumulative")) %>% 
    pivot_longer(cols = starts_with("Cumulative")) %>% 
    select(-.draw)
  
  observed <- model_objects$lumped_ca_covid %>% 
    ungroup()%>%
    mutate("Cumulative Incidence" = cumsum(new_cases),
           "Cumulative Death" = cumsum(new_deaths)) %>% 
    select(t = end_date, "Cumulative Incidence", "Cumulative Death") %>% 
    pivot_longer(-t, values_to = "observed")
  
  list(posterior = posterior, observed = observed)
}

multi_spread_draws <- function(...) {
  oc_data <- oc_prior %>%
    spread_draws(...) %>% 
    mutate(.dist = "Prior") %>% 
    bind_rows(oc_post %>%
                spread_draws(...) %>% 
                mutate(.dist = "Posterior")) %>% 
    mutate(.loc = "oc")
  
  ca_data <- ca_prior %>%
    spread_draws(...) %>% 
    mutate(.dist = "Prior") %>% 
    bind_rows(ca_post %>%
                spread_draws(...) %>% 
                mutate(.dist = "Posterior")) %>% 
    mutate(.loc = "ca")
  
  bind_rows(oc_data, ca_data)
}

multi_gather_draws <- function(...) {
  oc_data <- oc_prior %>%
    gather_draws(...) %>% 
    mutate(.dist = "Prior") %>% 
    bind_rows(oc_post %>%
                gather_draws(...) %>% 
                mutate(.dist = "Posterior")) %>% 
    mutate(.loc = "oc")
  
  ca_data <- ca_prior %>%
    gather_draws(...) %>% 
    mutate(.dist = "Prior") %>% 
    bind_rows(ca_post %>%
                gather_draws(...) %>% 
                mutate(.dist = "Posterior")) %>% 
    mutate(.loc = "ca")
  
  bind_rows(oc_data, ca_data)
}