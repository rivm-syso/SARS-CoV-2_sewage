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
    mutate( date = as.factor(date),
            age_group = as.factor(age_group))
}

# Nog even goed kijken naar startdatum en einddatum hierin verwerken
calc_vax <- function(df_vaccins,df_ziekenhuisopnames){
  df_vaccins <- df_vaccins %>%
    filter(between(as.Date(Week_1eprik),startday,lastday)) %>%
    mutate(Week_1eprik = format(Week_1eprik,"%G-%V")) %>%
    select("week" = "Week_1eprik",
           "municipality" = "Gemeente",
           "age_group" = "Leeftijdsgroep5",
           "age_population" = "Populatie",
           "percentage_vax" = "Vaccinatiegraad_coronit_cims") %>%
    # Percentages kunnen nooit meer dan 1 zijn.
    mutate(percentage_vax = if_else(percentage_vax > 100,1,percentage_vax/100)) %>%
    filter(age_group != "Niet vermeld")
  
  df_ziekenhuisopnames <- df_ziekenhuisopnames %>%
    filter(between(as.Date(Date_of_statistics),startday,lastday)) %>%
    select("date" = "Date_of_statistics",
           "municipality" = "Gemeente",
           "age_group" = "Leeftijdsgroep5",
           "hospitalizations" = "Hospital_admission") %>%
    mutate(week = format(as.Date(date),"%G-%V")) %>%
    # Eerst voegen we de populaties toe aan de vaccinatiedata
    left_join(df_vaccins %>% select(municipality,age_group,age_population) %>% unique(),
                                by = c("municipality","age_group")) %>%
    # Daarna voegen we de vaccinaties toe, en vullen die aan met nullen
    left_join(df_vaccins %>% select(municipality,age_group,week,percentage_vax) %>% unique(),
              by = c("municipality","age_group","week")) %>%
    mutate(percentage_vax = if_else(is.na(percentage_vax),0,percentage_vax)) %>%
    select( - week) %>%
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
    hosp_rate_age = rep(.5, length( unique( df_muni$age_group)) - 1),
    prevention_vax = rep(.8, length( unique( df_muni$age_group)))
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
stan_split <- function(fit_hospitalization,num_groups,parameters){
  
  places <- 1:fit_hospitalization@sim[["n_flatnames"]]
  
  place_parameters <- NULL
  place_rest <- "log_likes_hospital"
  for(i in seq_len(length(parameters))){
    # Find the positions of the expected hospitalizations
    place_parameters[[i]] <- places[str_detect(fit_hospitalization@sim[["fnames_oi"]],
                                               parameters[i])]
    place_rest <- paste0(place_rest,"|",parameters[i])
  }
  
  # Find the positions of the log_likes_hospital to get rid of most entries
  place_log_like <- places[str_detect(fit_hospitalization@sim[["fnames_oi"]],
                                      "log_likes_hospital")]
  
  # Find the postion of all other variables
  place_rest <- places[str_detect(fit_hospitalization@sim[["fnames_oi"]],
                                  place_rest,negate = T)]
  
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
    
    place_hosp <- c(place_log_like[1],place_rest)
    for(i in seq_len(length(parameters))){
      place_hosp <- c(place_hosp,place_parameters[[i]][places_sub])
    }
    # We find the places of the results we keep and those we discard
    place_hosp <- sort(place_hosp)
    place_discard <- sort(setdiff(places,place_hosp))
    
    # The dimension of the expected/simulated hospitalizations is n_days and n_muni
    # fit_hosp_list_new@par_dims[["expected_hospitalizations"]][[2]] <- n_muni
    # fit_hosp_list_new@par_dims[["simulated_hospitalizations"]][[2]] <- n_muni
    # fit_hosp_list_new@par_dims[["log_likes_hospital"]] <- c(1,1)
    
    for(i in 1:n_chains){
      fit_hosp_list_new@sim[["samples"]][[i]][place_discard] <- NULL
    }
    
    for(i in seq_len(length(parameters))){
      fit_hosp_list_new@sim[["dims_oi"]][[parameters[i]]][[2]] <- n_muni
    }
    fit_hosp_list_new@sim[["dims_oi"]][["log_likes_hospital"]] <- 
      rep(1,length(fit_hosp_list_new@sim[["dims_oi"]][["log_likes_hospital"]]))
    
    
    fit_hosp_list_new@sim[["fnames_oi"]] <- fit_hosp_list_new@sim[["fnames_oi"]][place_hosp]
    fit_hosp_list_new@sim[["n_flatnames"]] <- length(fit_hosp_list_new@sim[["fnames_oi"]])
    
    fit_hosp_list[[counter]] <- fit_hosp_list_new
  }
  
  return(fit_hosp_list)
}



