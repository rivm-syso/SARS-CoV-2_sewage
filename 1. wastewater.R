library(tidyverse)
library(here)
library(rstan)
library(loo)
library( tidybayes )

###
# User Settings
###
startday <- as.Date("2020-09-01")
lastday <- Sys.Date()

num_knots     <- 10
spline_degree <- 3
n_chains <- 10
setwd( here() )
source( "functions.R" )

###
# Model Run
###

#options(buildtools.check = NULL)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load( "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Generated data/df_viralload_human_regions.RData")
lastday <- min( lastday, max(df_viralload_human_regions$Datum))
startday <- max( startday, min(df_viralload_human_regions$Datum))
    
df_fractions <- df_viralload_human_regions %>% 
  select( municipality, rwzi=RWZI, municipality_pop=Inwoneraantal_municipality, starts_with( "frac" )) %>% 
  unique() 

# prepare data for Stan
df_sewage <- read_df_sewage( df_viralload_human_regions )

# make Stan datalist
datalist <- compose_data( 
  df_sewage %>% 
    select( rwzi, concentration, date ),
  spline_degree = spline_degree,
  num_knots = num_knots,
  mode = 0 # determines sampling temperature
)

# run Stan model
fit <- stan(
  "wastewater.stan",
  model_name = "wastewater_model",
  data = datalist,
  init = initials,
  chains = n_chains,
  iter = 80,
  thin = 5,
  control = list(adapt_delta = 0.7, max_treedepth = 12),
  pars = c("k", "x0", "sigma_observations", "RWvar", "load" )
)

traceplot(fit, pars = c("k", "x0", "sigma_observations", "RWvar"))

df_posteriors <- fit %>% 
  recover_types( df_sewage ) %>% 
  spread_draws( load[date,rwzi],
                x0, k ) %>% 
  left_join( df_sewage ) # Get original data back in


save(fit, file = str_c( outdir_out, "fit_pspline_", Sys.Date(), ".rda"))
save(df_posteriors, df_fractions, file = str_c( outdir_out, "posteriors_", Sys.Date(), ".rda"))

# model selection based on expected predictive performance
# notice that loo_ic does not perform well, presumably bc these
# are timeseries data. for the moment, rely on waic
loo_output = loo(fit, cores = 10, is_method = "psis")
loo_output
LL <- extract_log_lik(fit, parameter_name='log_likes_water', merge_chains=F)
waic(LL)
r_eff <- relative_eff(exp(LL), cores = n_chains)
loo(LL, r_eff=r_eff)
rm(r_eff) # big thing, remove whenever possible
