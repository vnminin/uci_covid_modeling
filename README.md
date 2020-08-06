# uci_covid_modeling

This github contains code for generating situation reports on the dynamics and future trends of COVID-19 in Orange County, CA. The manuscript associated with this work is [here](), and the current Orange County situation report is [available here](https://vnminin.github.io/uci_covid_modeling/). 


## Methodology
Our analysis relies on a mechanistic six compartment model of the COVID-19 epidemic. We then use Bayesian inference to provide inference on key disease dynamics and make predictions on future observed cases and deaths. Further descriptions of the model are available in the manuscript linked above. 

We use rStan to conduct Bayesian inference. Here is our [model code](https://github.com/vnminin/uci_covid_modeling/tree/master/code/SEIeIpRD). 

To replicate our analysis, use the [website code](https://github.com/vnminin/uci_covid_modeling/blob/master/analysis/index.Rmd). 

To conduct your own analysis of COVID-19 trends using our model you will need the following data:

1. Daily number of reported cases
2. Daily number of tests (viral not antibody)
3. Daily number of reported deaths

Due to reporting delays, we recommend limiting the end date of your data set to well before the actual end date of available data, as past counts are often updated retroactively weeks later. In our analysis we then condensed this daily data into three day periods. You can adjust this as needed. You will also need to specify a population size. It is likely you will want to adjust the priors for initial conditions, this can be done in create_model_objects (below).


## Example code
Below is a brief example of how to run our model using functions available in the code folder of this repository.
```
data_set <- data

start_date <- "YYYY-mm-dd"
end_date <- "YYYY-mm-dd"

#specify the beta parameters for the prior for the fraction in the population not susceptible
frac_not_susceptible <- c(x, y)

#create the list which will be put into stan
model_object <- create_model_objects(data_set, 
                                      first_day=start_date, 
                                      last_day=end_date,
                                      time_interval_in_days=3
                                      popsize=3175692L, 
                                      frac_carrs_beta=frac_not_susceptible)


#set initial values and controls as needed, and run the code
oc_post <- stan(file = here("code", SEIeIpRD, str_c(SEIeIpRD, ".stan")),
                data = model_objects,
                init = init_func,
                seed = 0,
                chains = 4,
                control = control_list)


```
This posterior can then be used for further analysis. Code for our specific analyses are available in the [website code](https://github.com/vnminin/uci_covid_modeling/blob/master/analysis/index.Rmd). 
