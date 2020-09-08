# uci_covid_modeling

This github contains code for generating situation reports on the dynamics and future trends of COVID-19 in Orange County, CA. The manuscript associated with this work is [here](https://arxiv.org/abs/2009.02654), and the current Orange County situation report is [available here](https://vnminin.github.io/uci_covid_modeling/). 

## Navigation
Files used to generate website pages are stored in the [analysis folder](https://github.com/vnminin/uci_covid_modeling/tree/master/analysis). The most recent report is generated with the [index](https://github.com/vnminin/uci_covid_modeling/blob/master/analysis/index.Rmd) file. Previous reports are generated by files corresponding to their date ranges. 

Functions used in website generation, stan code, and stan output are stored in the [code folder](https://github.com/vnminin/uci_covid_modeling/tree/master/code). 

Aggregated data used as the model input are available in the [data folder](https://github.com/vnminin/uci_covid_modeling/tree/master/data)


## Data
We use data aggregated from data from provided by Orange County, California Health Care Agency. Crucially, we exclude repeat tests given to patients who have tested positive for COVID-19 (which happens when patients are hospitalized). Our data may not match publicly available sources. 

## Methodology
Our analysis relies on a mechanistic six compartment model of the COVID-19 pandemic. We then use Bayesian inference to provide inference on key disease dynamics and make predictions on future observed cases and deaths. Further descriptions of the model are available in the manuscript linked above. 

We use the [Stan](https://mc-stan.org/) platform for our Bayesian inference. You can download [rStan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) to use Stan inside of R. See the [stan code](https://github.com/vnminin/uci_covid_modeling/tree/master/code/SEIeIpRD) for modeling and inference details. 

To replicate our analysis, use the [website code](https://github.com/vnminin/uci_covid_modeling/blob/master/analysis/index.Rmd). 

To conduct your own analysis of COVID-19 trends using our model you will need the following data:

1. Daily number of reported cases
2. Daily number of tests (viral not antibody)
3. Daily number of reported deaths

Due to reporting delays, we recommend limiting the end date of your data set to well before the actual end date of available data, as past counts are often updated retroactively weeks later. In our analysis we then condensed this daily data into three day periods. You can adjust this as needed. You will also need to specify a population size. It is likely you will want to adjust the priors for initial conditions, this can be done in create_model_objects (below).


## Example code
Below is a brief example of how to run our model to generate posterior distributions using functions available in the code folder of this repository.
```
library(rstan)
library(here)

model_name <- "SEIeIpRD"
loc_name <- "oc"

source(here("code", model_name, "helper_functions.R"))

data_set <- data

start_date <- "YYYY-mm-dd"
end_date <- "YYYY-mm-dd"

#specify the beta parameters for the initial condition priors, and popsize
frac_not_susceptible <- c(x, y)
frac_infected_infectious <-c(w,z)
frac_infectious_early <-c(u,v)

new_popsize <- N


#create the list which will be put into stan
model_object <- create_model_objects(data_set, 
                                      first_day=start_date, 
                                      last_day=end_date,
                                      time_interval_in_days=3
                                      popsize=new_popsize, 
                                      frac_carrs_beta=frac_not_susceptible,
                                      frac_carrs_infec_beta=frac_infected_infectious,
                                      frac_infec_early_beta=frac_infectious_early)


#set initial values and controls as needed, and run the code
oc_post <- stan(file = here("code", SEIeIpRD, str_c(SEIeIpRD, ".stan")),
                data = model_objects,
                init = init_func,
                seed = 0,
                chains = 4,
                control = control_list)


```
This posterior can then be used for further analysis. Code for our specific analyses are available in the [website code](https://github.com/vnminin/uci_covid_modeling/blob/master/analysis/index.Rmd). 

## Citation
Fintzi J, Bayer D, Goldstein I, Lumbard K, Ricotta E, Warner S, Busch LM, Strich JR, Chertow DS, Parker DM, Boden-Albala B, Dratch A, Chhuon R, Quick N, Zahn M, Minin VN. Using multiple data streams to estimate and forecast SARS-CoV-2 transmission dynamics, with application to the virus spread in Orange County, California, https://arxiv.org/abs/2009.02654.
