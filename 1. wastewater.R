library(tidyverse)
library(here)
library(rstan)
library(loo)
library(tidybayes)
library(furrr)
library(knitr) # For writing parameters to html file

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
  "wastewater.stan",
  model_name = "wastewater_model",
  data = datalist,
  init = initials,
  chains = n_chains,
  warmup = 400,
  iter = 800,
  thin = 4,
  control = list(adapt_delta = 0.7, max_treedepth = 12),
  pars = c("k", "x0", "sigma_observations", "RWvar", "load", "a_population" )
)

print(traceplot(fit, pars = c("k", "x0", "sigma_observations", "RWvar")))

summary( fit, pars = c("k", "x0", "sigma_observations", "RWvar", "a_population") ) %>% 
  pluck( "summary" ) %>%
  as_tibble( rownames = "parameter" ) %>% 
  select( parameter, mean, `2.5%`, `97.5%`, n_eff, Rhat ) %>% 
  kable( format="html", digits = 3 ) %>% 
  cat( file = here( runname, "output", "model_data", "parameters.html") )

df_posteriors <- fit %>%
  recover_types( df_sewage ) %>%
  stan_split(10,c("load"),c("load_population"),startday,lastday) %>%
  future_map(function(x){spread_draws(x,load[date,rwzi], x0, k )}) %>%
  bind_rows() %>%
  ungroup() %>%
  mutate( date=as.Date(date) )


save(fit, file = here( runname, "output", "model_data", "fit_pspline.RData"))
save(df_posteriors, df_fractions, df_sewage, file = here( runname, "output", "model_data", "posteriors.RData"))

# Clean up the enviroment by removing objects we no longer need
rm(datalist,fit)
invisible(gc()) # Just to be sure as fit can be rather large
