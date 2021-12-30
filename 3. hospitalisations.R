### Poisson regression for hospital data

library( tidyverse )
library( here )
library( rstan )
library( tidybayes )
library( loo )
library( furrr )

setwd( here() )

if(!exists("functions_sourced")){
  source( "0. functions.R" )
}
if(!exists("settings_sourced")){
  source( "0. settings.R" )
}

#options(buildtools.check = NULL)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Number of cores for furrr
plan(multisession, workers = 10)

#
# Load files
# 
load_if_needed( "df_viralload_human_regions", viralload_filename )
load_if_needed( list("df_posteriors","df_fractions"),
                here( runname, "output", "model_data", "posteriors.Rdata" ) )

df_vaccins <- calc_vax(startday,lastday) # Calculate the percentage of vaccinated individuals per municipality

# Clip first and last day such that all data frames span the same time period
lastday  <- min( lastday,  max(df_posteriors$date), max(df_vaccins$date), max(df_viralload_human_regions$Datum))
startday <- max( startday, min(df_posteriors$date), min(df_vaccins$date), min(df_viralload_human_regions$Datum))

# Calculate median load per municipality from posterior
#  also sums up the population in municipalities
df_muni <- calc_df_muni(df_posteriors, df_vaccins, startday,lastday, 20 )

rm( df_posteriors )

# Save calculated data frames per age group
# save(df_vaccins,file = here( runname, "output", "model_data", "df_vaccins_age.RData"))
save(df_muni,file = here( runname, "output", "model_data", "df_muni_age.RData"))

# TODO: is this really needed?
future:::ClusterRegistry("stop")

# run Stan model
fit_hospitalization = stan(
  "hospitalizations.stan",
  model_name = "wastewater_model",
  data = compose_data(
    max_delay = max_delay,
    ref_load = ref_load,
    delay_vax = delay_vax,
    df_muni %>% mutate( date=as.factor(date))),
  init = initials_hosp,
  chains = 4,
  iter = 200,
  refresh = 10,
  control = list(adapt_delta = .95, max_treedepth = 12)
)

print(traceplot(fit_hospitalization,"mean_hosp_rate"))
print(traceplot(fit_hospitalization,"sigma_hosp_rate"))
print(traceplot(fit_hospitalization,"prevention_vax"))

df_posteriors_hosp <- fit_hospitalization %>%
  recover_types(df_muni) %>%
  spread_draws(hosp_rate[age_group,municipality],
               prevention_vax[age_group])

# model selection based on predictive performance 
# out-of-sample prediction for time series a la Vehtari possible?
# waic ok, loo_ic meh - sort this out!
# loo_output = loo(fit_hospitalization, cores = 10, is_method = "psis")
# loo_output
# LL <- extract_log_lik(fit_hospitalization, parameter_name='log_likes_hospital', merge_chains=F)
# waic(LL)
# r_eff <- relative_eff(exp(LL), cores = 10) # costly, take no more than 1000-2000 samples
# loo(LL, r_eff=r_eff)
# rm(r_eff)

save(fit_hospitalization, file = here( runname, "output", "model_data", "fit_hosp_age.RData"))
save(df_posteriors_hosp, df_fractions, file = here( runname, "output", "model_data", "posteriors_hosp_age.RData"))

