library(tidyverse)
library(here)
library(rstan)
library(loo)
library(tidybayes)
library(furrr)

setwd( here() )
if(!exists("functions_sourced")){
  source( "0. functions.R" )
}
if(!exists("settings_sourced")){
  source( "0. settings.R" )
}

###
# Model Run
###

#options(buildtools.check = NULL)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load( viralload_filename )
lastday <- min( lastday, max(df_viralload_human_regions$Datum))
startday <- max( startday, min(df_viralload_human_regions$Datum))

df_fractions <- df_viralload_human_regions %>% 
  select( municipality, rwzi=RWZI, vr = VRname, municipality_pop=Inwoneraantal_municipality, starts_with( "frac" )) %>% 
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
  here(runname,"model_code_backup","wastewater.stan"),
  model_name = "wastewater_model",
  data = datalist,
  init = initials,
  chains = n_chains,
  warmup = 400,
  iter = 800,
  thin = 4,
  control = list(adapt_delta = 0.7, max_treedepth = 12),
  pars = c("k", "x0", "sigma_observations", "RWvar", "load" )
)

print(traceplot(fit, pars = c("k", "x0", "sigma_observations", "RWvar")))

df_posteriors <- fit %>%
  recover_types( df_sewage ) %>%
  stan_split(10,c("load"),c("a_individual","load_population","log_likes_water")) %>%
  future_map(function(x){spread_draws(x,load[date,rwzi], x0, k )}) %>%
  bind_rows() %>%
  ungroup() %>%
  mutate( date=as.Date(date) )


save(fit, file = here( runname, "output", "model_data", "fit_pspline.RData"))
save(df_posteriors, df_fractions, df_sewage, file = here( runname, "output", "model_data", "posteriors.RData"))

# Removed model selection for now.

# model selection based on expected predictive performance
# notice that loo_ic does not perform well, presumably bc these
# are timeseries data. for the moment, rely on waic
# loo_output = loo(fit, cores = 10, is_method = "psis")
# loo_output
# LL <- extract_log_lik(fit, parameter_name='log_likes_water', merge_chains=F)
# waic(LL)
# r_eff <- relative_eff(exp(LL), cores = n_chains)
# loo(LL, r_eff=r_eff)
# rm(r_eff) # big thing, remove whenever possible
