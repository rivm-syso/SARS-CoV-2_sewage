### Poisson regression for hospital data

library( tidyverse )
library( here )
library( rstan )
library( tidybayes )
library( loo )
library( furrr )

setwd( here() )

source( "functions.R")
source( "settings.R" )

#options(buildtools.check = NULL)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

plan(multisession, workers = 10)

load( viralload_filename )
lastday <- min( lastday, max(df_viralload_human_regions$Datum))
startday <- max( startday, min(df_viralload_human_regions$Datum))

# only when not exist
load( here( outdir_out, "model_data", "fit_pspline_2021-10-28.rda" ) )
load( here( outdir_out, "model_data", "posteriors_2021-10-28.rda" ) )

df_fractions <- df_viralload_human_regions %>% 
  select( municipality, rwzi=RWZI, municipality_pop=Inwoneraantal_municipality, starts_with( "frac" )) %>% 
  unique()

# Calculate the percentage of vaccinated individuals per municipality
df_vaccins <- calc_vax(startday,lastday)

# Calculate median load per municipality from posterior
#  also sums up the population in municipalities
df_muni <- calc_df_muni(df_posteriors, df_vaccins, startday,lastday)

# Sometimes the vaccin file is not up to date
df_muni <- df_muni %>% filter( !is.na(population))
  
rm( df_posteriors )

# Save df_muni en df_vaccins
save(df_vaccins,file = here( outdir_out, "model_data", str_c("df_vaccins_age", Sys.Date(),".RData")))
save(df_muni,file = here( outdir_out, "model_data", str_c("df_muni_age", Sys.Date(),".RData")))

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
               prevention_vax[age_group]) #%>%
  # slice_sample(n = 100) %>%
  # right_join(df_muni) %>%
  # # construct the expected and simulated hospitalizations
  # group_by(municipality,age_group) %>%
  # group_split() %>%
  # lapply(function(df){
  #   arrange(df,date) %>%
  #     mutate(percentage_vax = lag(percentage_vax, n = 14, default = 0),
  #            expected_hospitalizations_cf = hosp_rate * 10^(load-19) * population,
  #            expected_hospitalizations = expected_hospitalizations_cf *
  #                                           (1 - prevention_vax*percentage_vax),
  #            simulated_hospitalizations = rpois(nrow(df),expected_hospitalizations),
  #            simulated_hospitalizations_cf = rpois(nrow(df),expected_hospitalizations_cf))
  # }) %>%
  # bind_rows() 

save(df_posteriors_hosp, df_muni,
     file = paste0(outdir_res,Sys.Date(), "df_posteriors.RData"))
save(fit_hospitalization, file = str_c( outdir_res, "fit_hosp", Sys.Date(), ".rda"))

df_posteriors_hosp <- fit_hospitalization %>% 
  recover_types( df_muni ) %>% 
  stan_split(10,c("expected_hospitalizations","simulated_hospitalizations"),
             c("hosp_parameter","sum_load")) %>%
  future_map(function(x){spread_draws(x,c(expected_hospitalizations,simulated_hospitalizations)[date,municipality])}) %>% 
  bind_rows() %>%
  left_join( df_muni, c("date", "municipality") )

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

# save, rmeove, and load fit
save(fit_hospitalization, file = here(outdir_out, "model_data", str_c( "fit_hosp_age", Sys.Date(), ".rda")))
save(df_posteriors_hosp, df_fractions, file = here(outdir_out, "model_data", str_c( "posteriors_hosp_age", Sys.Date(), ".rda")))

