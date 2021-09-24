# directories
outdir_fig = "./figures/"
if( !dir.exists(outdir_fig)) dir.create(outdir_fig)
outdir_out = "./output/"
if( !dir.exists(outdir_out)) dir.create(outdir_out)
outdir_res = "./results/"
if( !dir.exists(outdir_res)) dir.create(outdir_res)

calc_df_muni <-function(df_posteriors){
  df_posteriors %>% 
    select( .draw, date, load, rwzi, municipality, hospitalizations ) %>%
    left_join(df_fractions, by = c("rwzi", "municipality") ) %>% 
    mutate( load_muni = frac_municipality2RWZI * load ) %>% 
    group_by( date, municipality, hospitalizations, .draw ) %>%
    summarize( load = sum( load_muni ), 
               municipality_pop = sum(municipality_pop),
               .groups="drop_last") %>% 
    group_by( date, municipality, hospitalizations, municipality_pop ) %>% 
    median_qi( load ) %>%
    mutate( date=as.Date(as.character(date)))
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
    hosp_rate = rep(2.5, length( unique( df_muni$municipality)))
  ))
}

# colorblind-friendly scheme
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

# the estimated probability of detection
probdetection = function(x, x0, k) {
  return(1 / (1 + exp((-k) * (x - x0))))
}