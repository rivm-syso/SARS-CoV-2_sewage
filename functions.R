# directories
outdir_fig = "./figures/"
if( !dir.exists(outdir_fig)) dir.create(outdir_fig)
outdir_out = "./output/"
if( !dir.exists(outdir_out)) dir.create(outdir_out)
outdir_res = "./results/"
if( !dir.exists(outdir_res)) dir.create(outdir_res)

calc_df_muni <-function(df_posteriors,startday,lastday){
  df_posteriors %>% 
    filter( as.Date(date) >= startday, as.Date(date) <= lastday ) %>%
    group_by( municipality,date ) %>% 
    sample_draws(10) %>%
    ungroup() %>% 
    select( .draw, date, load, rwzi, municipality ) %>%
    left_join(df_fractions, by = c("rwzi","municipality") ) %>% 
    mutate( load_muni = frac_municipality2RWZI * load) %>% 
    group_by( date, municipality, .draw ) %>%
    summarize( load = sum( load_muni ),
               .groups="drop_last") %>% 
    group_by( date, municipality ) %>% 
    # Use parallel computing for speed
    group_split() %>%
    future_map(function(df){summarize(df, load = median(load), # Median_qi can also be used
                                      # Then we first have to group df, and then apply 
                                      # median_qi. Does make the code a bit slower.
                                      date = first(date),
                                      municipality = first(municipality))} ) %>%
    bind_rows() %>%
    # Add the vaccinations, age, and hospitalizations
    left_join(df_vaccins, by = c("municipality","date")) %>%
    group_by(date,municipality) %>%
    summarize(date= first(date),
              load = first(load),
              hospitalizations = sum(hospitalizations),
              percentage_vax = sum(percentage_vax * population) / sum(population),
              population = sum(population)) %>%
    ungroup()
}

calc_vax <- function(df_vaccins,df_ziekenhuisopnames,startday,lastday){
  df_vaccins <- df_vaccins %>%
    filter(between(as.Date(Datum_1eprik),startday,lastday)) %>%
    select("date" = "Datum_1eprik",
           "municipality" = "Gemeente",
           "age_group" = "Leeftijdsgroep5",
           "population" = "Populatie",
           "percentage_vax" = "Vaccinatiegraad_coronit_cims") %>%
    # Percentages kunnen nooit meer dan 1/100% zijn.
    mutate(percentage_vax = if_else(percentage_vax > 100,1,percentage_vax/100)) %>%
    filter(age_group != "Niet vermeld")
  
  # We match each municipality and age_group with the same number of population
  # independent of the day, hence we make a help-tibble with the populations
  df_vaccins_pop <- df_vaccins %>%
    group_by(municipality,age_group) %>%
    summarize(population = first(population)) %>%
    ungroup()
  
  # Once we matched the population, we no longer need the population data 
  # in the vaccination data
  df_vaccins <- df_vaccins %>%
    select(-population)
  
  df_ziekenhuisopnames <- df_ziekenhuisopnames %>%
    filter(between(as.Date(Date_of_statistics),startday,lastday)) %>%
    select("date" = "Date_of_statistics",
           "municipality" = "Gemeente",
           "age_group" = "Leeftijdsgroep5",
           "hospitalizations" = "Hospital_admission") %>%
    # Eerst voegen we de populaties toe aan de vaccinatiedata
    left_join(df_vaccins_pop,by = c("municipality","age_group")) %>%
    # Daarna voegen we de vaccinaties toe, en vullen die aan met nullen
    left_join(df_vaccins,by = c("municipality","age_group","date")) %>%
    mutate(percentage_vax = if_else(is.na(percentage_vax),0,percentage_vax)) %>%
    # Namen verbeteren
    mutate(municipality = str_remove(municipality,"^\'"),           # Den Haag & Bosch 
           municipality = str_remove(municipality," \\(O\\.\\)"),   # Hengelo
           municipality = str_replace(municipality,",","."),        # Nuenen etc.
           municipality = str_replace(municipality,"â","a"),        # Fryslan
           municipality = str_replace(municipality,"ú","u"))        # Sudwest Fryslan

  return(df_ziekenhuisopnames)
}

read_df_sewage <- function( df_viralload_human_regions ){
  df_viralload_human_regions %>%
    mutate(rwzi = as.factor( RWZI ),
           date = Datum,
           hospitalizations = Ziekenhuisopnames_NICE,
           concentration = case_when(
             N12_100000_RWZI != 0    ~ log10(N12_100000_RWZI), 
             N12_100000_RWZI == 0    ~ 0,
             is.na(N12_100000_RWZI)  ~ -1),
           rwzi_persons = Inwoneraantal,
           municipality_persons = Inwoneraantal_municipality,
           municipality = as.factor( municipality ),
           vr = as.factor( VRname ),
           .keep="none") %>% 
    filter(date >= startday, date <= lastday) %>% 
    mutate( date=as.factor(date) ) 
}

# initialize variables
initials = function() {
  list(
    k = 6.0,
    x0 = 12.0,
    sigma_observations = 0.35,
    RWvar = 0.35,
    a_population = c(11, rep(0.2, num_knots + spline_degree - 2 )),
    a_individual = matrix(0.0, nrow = length(levels(df_sewage$rwzi)), 
                          ncol = num_knots + spline_degree - 1)
  )
}


# initialize variables
initials_hosp = function() {
  return(list(
    mean_hosp_rate = 3.0,
    sigma_hosprate = 2,
    hosp_rate = rep(2.5, length( unique( df_muni$municipality ))),
    prevention_vax = 0.8
    # hosp_rate_age = rep(.5, length( unique( df_muni$age_group)) - 1),
    # prevention_vax = rep(.8, length( unique( df_muni$age_group)))
  ))
}
# 
# df_muni <- df_posteriors %>% 
#   select( .draw, date, load, rwzi, municipality, hospitalizations ) %>%
#   left_join(df_fractions ) %>% 
#   mutate( load_muni = frac_municipality2RWZI * load ) %>% 
#   group_by( date, municipality, hospitalizations, .draw ) %>%
#   summarize( load = sum( load_muni ), .groups="drop_last") %>% 
#   median_qi( load ) %>%
#   mutate( date=as.Date(as.character(date))) 


# colorblind-friendly scheme
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

# the estimated probability of detection
probdetection = function(x, x0, k) {
  return(1 / (1 + exp((-k) * (x - x0))))
}

# Split the model fit for future map
stan_split <- function(fit_hospitalization,num_groups,parameters,par_ignore){
  # We split the fit of the model, fit_hospitalizations in multiple groups along
  # the second dimension of the parameters of interest.
  # num_groups is the number of groups we create
  # parameters is a vector of strings of the parameters of interest, i.e. the
  # parameters we want to extract from the model
  # par_ignore is a vector of strings of parameters which are also of high dimension.
  # To make sure we can parallize spread_draws, we already throw out most of these results.
  places <- 1:fit_hospitalization@sim[["n_flatnames"]]
  
  # We create to integer vectors where the results from the parameters are. place_rest
  # will be a vector of all other results.
  place_parameters <- NULL
  place_log_like <- NULL
  place_rest <- NULL
  
  # We find the positions of the parameters we want to extract from the model
  for(i in seq_len(length(parameters))){
    place_parameters[[i]] <- places[str_detect(fit_hospitalization@sim[["fnames_oi"]],
                                               parameters[i])]
    place_rest <- c(place_rest,place_parameters[[i]])
  }
  
  # We find the positions of the parameters of which we want to forget most entries
  for(i in seq_len(length(par_ignore))){
    place_log_like[[i]] <- places[str_detect(fit_hospitalization@sim[["fnames_oi"]],
                                               par_ignore[i])]
    place_rest <- c(place_rest,place_log_like[[i]])
  }

  # Find the postion of all other variables
  place_rest <- sort(setdiff(places,place_rest))
  
  # We divide the results in (almost) equally sized lists
  n_chains <- fit_hospitalization@sim[["chains"]]
  n_days <- fit_hospitalization@par_dims[[parameters[1]]][[1]]
  n_muni <- floor(fit_hospitalization@par_dims[[parameters[1]]][[2]]/num_groups)
  muni_extra <- 0
  
  fit_hosp_list <- NULL
  
  for(counter in 1:num_groups){
    # We create a smaller stan-object
    fit_hosp_list_new <- fit_hospitalization
    
    if(counter == num_groups){
      muni_extra <- fit_hospitalization@par_dims[[parameters[1]]][[2]] - num_groups*n_muni
    }
    # We only keep a subset of the expected and simulated hospitalizations
    places_sub <- ((counter-1)*n_muni*n_days + 1):(counter*n_muni*n_days + muni_extra*n_days)
    # In case of the final worker, we have a bit more municipalities
    n_muni <- n_muni + muni_extra
    
    # We find the place of this group only
    place_hosp <- place_rest
    for(i in seq_len(length(parameters))){
      place_hosp <- c(place_hosp,place_parameters[[i]][places_sub])
    }
    for(i in seq_len(length(par_ignore))){
      place_hosp <- c(place_hosp,place_log_like[[i]][1])
    }
    
    # We find the places of the results we keep and those we discard
    place_hosp <- sort(place_hosp)
    place_discard <- sort(setdiff(places,place_hosp))

    for(i in 1:n_chains){
      fit_hosp_list_new@sim[["samples"]][[i]][place_discard] <- NULL
    }
    
    for(i in seq_len(length(parameters))){
      fit_hosp_list_new@sim[["dims_oi"]][[parameters[i]]][[2]] <- n_muni
    }
    
    for(i in seq_len(length(par_ignore))){
      fit_hosp_list_new@sim[["dims_oi"]][[par_ignore[i]]] <- 
        rep(1,length(fit_hosp_list_new@sim[["dims_oi"]][[par_ignore[1]]]))
    }
    
    
    fit_hosp_list_new@sim[["fnames_oi"]] <- fit_hosp_list_new@sim[["fnames_oi"]][place_hosp]
    fit_hosp_list_new@sim[["n_flatnames"]] <- length(fit_hosp_list_new@sim[["fnames_oi"]])
    
    fit_hosp_list[[counter]] <- fit_hosp_list_new
  }
  
  return(fit_hosp_list)
}



