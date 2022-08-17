load_if_needed <- function( object, filename ){
# Object becomes a list in case we need to load multiple file from the same file
    if( object %>% 
      lapply(function(object){
        length(ls( pattern=object, envir = .GlobalEnv)) == 0}) %>%
      c(recursive = T) %>% any())
    load( filename, envir = .GlobalEnv )
}


download_hospitalization <- function(startday, lastday, nationwide = T){
  
  download.file("https://data.rivm.nl/covid-19/COVID-19_ziekenhuisopnames.csv",
                destfile = here(runname, "output", "model_data", "hospitalizations.csv"),
                mode = "wb")
  
  df_hospital <- read.csv2(here(runname, "output", "model_data", "hospitalizations.csv")) %>%
    select("date" = "Date_of_statistics", "hospitalization" = "Hospital_admission",
           "municipality" = "Municipality_name") %>%
    # There are a few difference in municipality names between the open data 
    # and our own data
    mutate(municipality = str_remove(municipality,"^\'"),           # Den Haag & Bosch 
           municipality = str_remove(municipality," \\(O\\.\\)"),   # Hengelo
           municipality = str_replace(municipality,",","."),        # Nuenen etc.
           municipality = str_replace(municipality,"â","a"),        # Fryslan
           municipality = str_replace(municipality,"ú","u")) %>%         # Sudwest Fryslan
    mutate(date = as.Date(date)) %>%
    filter(between(date,startday,lastday))
  
  # Aggregate the results if we want the nationwide number of hospitalizations
  if(nationwide){
    df_hospital <- df_hospital %>% 
      group_by(date) %>% 
      summarize(hospitalization = log10(sum(hospitalization)), .groups = "drop")
  } else {
    df_hospital <- df_hospital %>%
      group_by(municipality) %>%
      arrange(date) %>%
      mutate(hospitalization = log10(mean_run(hospitalization,7,-3)))
  }
  
  return(df_hospital)
  
}

download_testresults <- function(startday, lastday){
  
  download.file("https://data.rivm.nl/covid-19/COVID-19_aantallen_gemeente_per_dag.csv",
                destfile = here(runname, "output", "model_data", "tests_results.csv"),
                mode = "wb")
  
  df_hospital <- read.csv2(here(runname, "output", "model_data", "tests_results.csv")) %>%
    select("date" = "Date_of_publication", "test" = "Total_reported",
           "municipality" = "Municipality_name") %>%
    # There are a few difference in municipality names between the open data 
    # and our own data
    mutate(municipality = str_remove(municipality,"^\'"),           # Den Haag & Bosch 
           municipality = str_remove(municipality," \\(O\\.\\)"),   # Hengelo
           municipality = str_replace(municipality,",","."),        # Nuenen etc.
           municipality = str_replace(municipality,"â","a"),        # Fryslan
           municipality = str_replace(municipality,"ú","u")) %>%         # Sudwest Fryslan
    mutate(date = as.Date(date)) %>%
    filter(between(date,startday,lastday))

  return(df_hospital)
  
}

calc_df_load_municipality <- function(df_posteriors,df_fractions){
  # We redistribute the loads per STP over the different municipalities
  df_posteriors %>%
    group_by(date) %>%
    group_split() %>%
    future_map(function(df){
      left_join(df,df_fractions, by = "rwzi") %>%
        mutate(load = frac_municipality2RWZI * 10^load) %>%
        group_by(municipality,date,.draw) %>%
        summarize(load = log10(sum(load)), .groups = "drop_last") %>%
        median_qi(load)
    }) %>%
    bind_rows()
}

read_df_sewage <- function( df_viralload_human_regions ){
  df_viralload_human_regions %>%
    mutate(rwzi = as.factor( RWZI ),
           date = Datum,
           concentration = case_when(
             RNA_100000_RWZI != 0    ~ log10(RNA_100000_RWZI), 
             RNA_100000_RWZI == 0    ~ 0,
             is.na(RNA_100000_RWZI)  ~ -1),
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
    a_individual = matrix(0.0, ncol = length(levels(df_sewage$rwzi)), 
                          nrow = num_knots + spline_degree - 1)
  )
}

# colorblind-friendly scheme
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

# the estimated probability of detection
probdetection = function(x, x0, k) {
  return(1 / (1 + exp((-k) * (x - x0))))
}

# Split the model fit for future map
stan_split <- function(fit_hospitalization,num_groups,parameters,par_ignore,
                            startday,lastday){
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
                                               paste0("^",parameters[i],"(\\[|$)"))]
    place_rest <- c(place_rest,place_parameters[[i]])
  }
  
  # We find the positions of the parameters of which we want to forget most entries
  for(i in seq_len(length(par_ignore))){
    place_log_like[[i]] <- places[str_detect(fit_hospitalization@sim[["fnames_oi"]],
                                               paste0("^",par_ignore[i],"(\\[|$)"))]
    place_rest <- c(place_rest,place_log_like[[i]])
  }

  # Find the postion of all other variables
  place_rest <- sort(setdiff(places,place_rest))
  
  # We divide the results in (almost) equally sized lists
  n_chains <- fit_hospitalization@sim[["chains"]]
  n_days <- as.numeric(lastday - startday) + 1
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
