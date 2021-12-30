load_if_needed <- function( object, filename ){
# Object becomes a list in case we need to load multiple file from the same file
    if( object %>% 
      lapply(function(object){
        length(ls( pattern=object, envir = .GlobalEnv)) == 0}) %>%
      c(recursive = T) %>% any())
    load( filename, envir = .GlobalEnv )
}


calc_df_muni <-function(df_posteriors, df_vaccins, startday,lastday,age = 5){
  # We create the data frame with waste water data and the viral load.
  # The optional input age can either be a single multiple of 5 which will 
  # determine equally sized age groups, or a vector with explicit age groups
  # of the form "5n - 5m-1".
  
  # Fractions of municipalities and VR's in RWZI's
  df_fractions <- df_viralload_human_regions %>%
    select( municipality, rwzi=RWZI, municipality_pop=Inwoneraantal_municipality, starts_with( "frac" )) %>%
    unique()
  
  
  df_muni <- df_posteriors %>% 
    filter( between( date, startday, lastday ) ) %>%
    group_by( municipality,date ) %>% 
    sample_draws(10) %>%
    ungroup() %>% 
    select( .draw, date, load, rwzi, municipality ) %>%
    left_join(df_fractions, by = c("rwzi","municipality") ) %>% 
    mutate( load_muni = frac_municipality2RWZI * load) %>% 
    group_by( date, municipality, .draw ) %>%
    summarize( load = sum( load_muni ),
               .groups="drop_last") %>% 
    # Use parallel computing for speed
    group_split() %>%
    future_map(function(df){summarize(df, load_sd = sd(load), 
                                      load = mean(load),# Median_qi can also be used
                                      # Then we first have to group df, and then apply 
                                      # median_qi. Does make the code a bit slower.
                                      date = first(date),
                                      municipality = first(municipality))} ) %>%
    bind_rows() %>%
    # Add the vaccinations, age, and hospitalizations
    left_join(df_vaccins, by = c("municipality","date")) %>%
    # Extract the lowest age from each age group so that we can reshape the groups
    mutate(age_group = str_extract(age_group,"^[0-9]+") %>% as.numeric())
    
    # Make the new age groups so that we can later group them together
    # Both approaches probably have better ways to achieve their goal
    if(is.numeric(age)){
      df_muni <- df_muni %>%
        mutate(age_group = paste0(floor(age_group/age)*age," - ",
                                  age*floor(age_group/age)+age-1)) 
    } else {
      for(i in seq_len(length(age))){
        age_low = str_extract(age[i],"^[0-9]+") %>% as.numeric()
        age_high = str_extract(age[i],"[0-9]+$") %>% as.numeric()
        df_muni <- df_muni %>%
          mutate(age_group = if_else(between(age_group,age_low,age_high),-i*1.,age_group))
      }
      df_muni <- df_muni %>% mutate(age_group = age[-age_group])
    }
    
  df_muni <- df_muni %>%
    group_by(date,municipality,age_group) %>%
    summarize(load = first(load),
              load_sd = first(load_sd),
              percentage_vax = sum(percentage_vax*population)/sum(population),
              hospitalizations = sum(hospitalizations), 
              population = sum(population),
              .groups="drop") %>%
    mutate(age_group = as.factor(age_group))
    
    return(df_muni)
}


calc_vax <- function( startday,lastday){
  df_vaccins <- read.csv(vaccin_filename) %>%
    filter(between(as.Date(Datum_1eprik),startday,lastday)) %>%
    select("date" = "Datum_1eprik",
           "municipality" = "Gemeente",
           "age_group" = "Leeftijdsgroep",
           "population" = "Populatie",
           "percentage_vax" = "Vaccinatiegraad_coronit_cims") %>%
    # Percentages kunnen nooit meer dan 1/100% zijn.
    mutate(percentage_vax = if_else(percentage_vax > 100,1,percentage_vax/100)) %>%
    filter(age_group != "Niet vermeld")
  
  # Update the final day to match with the vaccin data
  lastday <- as.Date(max(df_vaccins$date))
  
  # We match each municipality and age_group with the same number of population
  # independent of the day, hence we make a help-tibble with the populations
  df_vaccins_pop <- df_vaccins %>%
    group_by(municipality,age_group) %>%
    summarize(population = first(population), .groups="drop")
  
  # Once we matched the population, we no longer need the population data 
  # in the vaccination data
  df_vaccins <- df_vaccins %>%
    select(-population)
  
  df_ziekenhuisopnames <- read.csv(hosp_filename) %>%
    filter(between(as.Date(AdmissionDate_Pid),startday,lastday)) %>%
    select("date" = "AdmissionDate_Pid",
           "municipality" = "Gemeente",
           "age_group" = "Leeftijdsgroep5",
           "hospitalizations" = "Hospital_admissions") %>%
    filter(!is.na(municipality) & (age_group != "Niet vermeld")) %>%
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
           municipality = str_replace(municipality,"ú","u")) %>%         # Sudwest Fryslan
    mutate( date=as.Date(date))
  
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
    # mean_hosp_rate = 3.0,
    # sigma_hosprate = 2,
    #hosp_rate = rep(2.5, length( unique( df_muni$municipality ))),
    hosp_rate_age = rep(1,length(levels(df_muni$age_group))),
    prevention_vax = rep(.8, length( levels( df_muni$age_group)))
  ))
}

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

functions_sourced <- T