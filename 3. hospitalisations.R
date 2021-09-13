### Poisson regression for hospital data

library(tidyverse)
library(here)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library( tidybayes )
library( loo )
library(furrr)

plan(multisession, workers = 10)

setwd( here() )

# select subset of data for regression analysis
# based on frequency of sampling and start of vaccination
# here we take September 2020 up to and including February 2021
startday <- as.Date("2020-09-01")
lastday <- as.Date("2021-08-15")     # 2021-04-12

#load( "./output/fit_pspline_2021-07-22.rda" )
load( "./output/posteriors_2021-08-25.rda")
load("df_viralload_human_regions.RData")

df_fractions <- df_viralload_human_regions %>% 
  select( municipality, rwzi=RWZI, starts_with( "frac" )) %>% 
  unique()

source( "functions.R")

# Calculate the percentage of vaccinated individuals per municipality
df_vaccins <- calc_vax(readxl::read_excel("Vaccinatiegraad.xlsx"),
                       read.delim("Ziekenhuisopnames.txt", sep = " "))
# Calculate median load per municipality from posterior
#  also sums up the population in municipalities
df_muni <- calc_df_muni(df_posteriors,startday,lastday)

# Save df_muni en df_vaccins
save(df_vaccins,file = paste0(Sys.Date(),"df_vaccins_age.RData"))
save(df_muni,file = paste0(Sys.Date(), "df_muni_age.RData"))

future:::ClusterRegistry("stop")

# run Stan model
fit_hospitalization = stan(
  "hospitalizations.stan",
  model_name = "wastewater_model",
  data = compose_data(
          max_delay = 0,
          ref_load = 19,
          delay_vax = 14,
          df_muni),
  init = initials_hosp,
  chains = 4,
  iter = 100,
  refresh = 10,
  control = list(adapt_delta = 0.7, max_treedepth = 12)
)

# some output and checks
traceplot(fit_hospitalization, pars = c("prevention_vax","hosp_rate[16]", "hosp_rate[249]", "hosp_rate[252]", "hosp_rate[292]"))

df_posteriors_hosp <- fit_hospitalization %>% 
  recover_types( df_muni ) %>% 
  stan_split(10,c("expected_hospitalizations","simulated_hospitalizations")) %>%
  future_map(function(x){spread_draws(x,c(expected_hospitalizations,simulated_hospitalizations)[date,municipality])}) %>% 
  bind_rows() %>%
  left_join( df_muni )


# model selection based on predictive performance 
# out-of-sample prediction for time series a la Vehtari possible?
# waic ok, loo_ic meh - sort this out!
loo_output = loo(fit_hospitalization, cores = 10, is_method = "psis")
loo_output
LL <- extract_log_lik(fit_hospitalization, parameter_name='log_likes_hospital', merge_chains=F)
waic(LL)
r_eff <- relative_eff(exp(LL), cores = 10) # costly, take no more than 1000-2000 samples
loo(LL, r_eff=r_eff)
rm(r_eff)

# save, rmeove, and load fit
save(fit_hospitalization, file = str_c( outdir_out, "fit_hosp_age", Sys.Date(), ".rda"))
save(df_posteriors_hosp, df_fractions, file = str_c( outdir_out, "posteriors_hosp_age", Sys.Date(), ".rda"))

