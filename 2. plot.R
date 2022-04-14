library( tidyverse )
library( tidybayes )
library( patchwork )
library( here )
library( furrr )

if(!exists("functions_sourced")){
  source( "0. functions.R" )
}
if(!exists("settings_sourced")){
  source( "0. settings.R" )
}

load_if_needed( "df_viralload_human_regions", viralload_filename )

lastday <- min( lastday, max(df_viralload_human_regions$Datum))
startday <- max( startday, min(df_viralload_human_regions$Datum))
# rm( df_viralload_human_regions ) # Needed to construct df_muni later on

load_if_needed( list("df_posteriors","df_fractions","df_sewage"),
                here( runname, "output", "model_data", "posteriors.Rdata" ) )


###
# Plot prob of detection
###
df_posteriors %>%
  filter(date == startday & rwzi == first(rwzi)) %>%
  expand_grid( conc=seq(11, 13, by = 0.1) ) %>%
  mutate( p_det = probdetection(conc, x0,k) ) %>% 
  ggplot( ) +
  geom_line( aes( x = conc, y = p_det, group= .draw ),
             alpha = .01, color = cbPalette[5] ) +
  scale_x_continuous( expression(paste("Log(load) (10"^"-5"," persons"^"-1",")")),
                      limits = c(11, 13),
                      expand = c(0, 0),  breaks = c(11, 11.5, 12, 12.5, 13)  ) +
  scale_y_continuous( "Prob(detection)",  limits = c(0, 1),
                      expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1)  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave( here(runname, "figures", "model_data", "prob_detection.png"), width = 6.5, height = 4.5, units = "in")
ggsave( here(runname, "figures", "model_data", "prob_detection.tiff"), width = 6.5, height = 4.5, units = "in")
ggsave( here(runname, "figures", "model_data", "prob_detection.svg"), width = 6.5, height = 4.5, units = "in")


###
# Plot the Netherlands and write CSV
###
df_posteriors %>%         # group_by(..) %>% group_split is more robust than group_split(..)
  group_by( date ) %>%    # if the tibble is already grouped then group_by overwrite those groups
  group_split() %>%       # whereas group_split only add new groups and keeps the original groups
  future_map(function(df){
    left_join(df,df_sewage %>% mutate(date = as.Date(date)), by = c("date","rwzi")) %>%
      group_by(date, .draw) %>%
      summarize( weighted_concentration = log10(weighted.mean(10^load,rwzi_persons )), .groups="drop" )
  }) %>%
  bind_rows() %>%
  group_by( date ) %>% 
  median_qi( weighted_concentration ) %T>%
  write.csv(here(runname,"output","Netherlands","posterior.csv")) %>%
  ggplot( aes( x = date, y = weighted_concentration,
               ymin = .lower, ymax = .upper ) ) +
  geom_ribbon(alpha = 0.25) +
  geom_line( color = cbPalette[5] ) +
#  labs( title = "Viral load in the Netherlands",
#       subtitle = "1 August 2020 - 8 February  2022") +
  coord_cartesian(ylim = c(11, 15)) +
  scale_x_date( "Date", breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") + 
  scale_y_continuous( expression(paste("Log(load) (10"^"-5"," persons"^"-1",")"))) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(color = cbPalette[5] ),
    plot.subtitle = element_text(color = cbPalette[5] ),
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave( here( runname,"Figures", "Netherlands", "posterior.png" ), width = 6.5, height = 4.5, units = "in")
ggsave( here( runname,"Figures", "Netherlands", "posterior.tiff" ), width = 6.5, height = 4.5, units = "in")
ggsave( here( runname,"Figures", "Netherlands", "posterior.svg" ), width = 6.5, height = 4.5, units = "in")

###
# Plot each RWZI, write CSV
###
df_posteriors %>% 
  left_join(df_sewage %>%
              select(date,rwzi,concentration) %>%
              unique() %>%
              mutate(date = as.Date(date)), by = c("date","rwzi")) %>%
  select( date, load, concentration, rwzi ) %>%
  group_by( rwzi) %>% 
  group_split() %>%
  future_map(function(x){
    group_by(x,date, rwzi,measurement = concentration) %>%
      median_qi(load) %>%
      mutate( zeromeasurement = ifelse(measurement == 0, load, -2),
              date=as.Date(as.character(date)),
              rwzi=as.character(rwzi))}, .options = furrr_options(seed = T)) %>%
  bind_rows() %>%
  group_by( rwzi ) %>% 
  group_split() %>%
  future_walk( function(x, y = x$rwzi[1]){
    p <- ggplot(x, mapping = aes(x = date,y = load,ymin = .lower, ymax=.upper) ) +
      geom_point(aes(y = measurement), color = cbPalette[5], size = 2.5, alpha = .5) +
      geom_point(aes(y = zeromeasurement), color = cbPalette[5], size = 2, shape = 1, alpha = .5) +
      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[5]) +
      coord_cartesian(ylim = c(11, 15)) +
      scale_x_date("Date",  breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") + #"%b"
      scale_y_continuous(expression(paste("Log(load) (10"^"-5"," persons"^"-1",")"))) +
      ggtitle(y) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(color = cbPalette[5]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none" )
    ggsave( here( runname, "figures", "RWZI", str_c( y, ".png")),
            plot = p, width = 6.5, height = 4.5, units = "in")
    
    write_csv(x, here( runname, "output", "RWZI", str_c( y, ".csv")))})

##
# make figure for manuscript
# largest installations are 16, 249, 252, 292
##
df_posteriors %>% 
  left_join(df_sewage %>%
              select(date,rwzi,rwzi_persons,concentration) %>%
              unique() %>%
              mutate(date = as.Date(date)), by = c("date","rwzi")) %>%
  # Filter dynamical on the largest STP instead of fixed positions
  filter(rwzi_persons %in%
               (df_sewage$rwzi_persons %>% 
               unique() %>% 
               sort(decreasing = T) %>%
               .[1:9])) %>%
  arrange(desc(rwzi_persons)) %>%
  select( date, load, concentration, rwzi ) %>%
  group_by(rwzi) %>%
  group_split() %>%
  future_map(function(x){
    group_by(x, date, rwzi, measurement=concentration ) %>% 
    median_qi( load ) %>%
    mutate( zeromeasurement = ifelse(measurement == 0, load, -2),
            date=as.Date(as.character(date)),
            rwzi=as.character(rwzi))})  %>%
  bind_rows() %T>%
  write.csv(here(runname,"output","manuscript","rwzi.csv")) %>%
  ggplot(mapping = aes(x = date,y = load,ymin = .lower, ymax=.upper) ) +
    geom_point(aes(y = measurement), color = cbPalette[5], size = 1, alpha = .5) +
    geom_point(aes(y = zeromeasurement), color = cbPalette[5], size = 1, shape = .5) +
    geom_ribbon(alpha = 0.25) +
    geom_line(color = cbPalette[5]) +
    coord_cartesian(ylim = c(11, 15)) +
    scale_x_date("Date",  breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") + #"%b"
    scale_y_continuous(expression(paste("Log(load) (10"^"-5"," persons"^"-1",")"))) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(color = cbPalette[5]),
      panel.grid.major = element_line(size = 0.7),
      panel.grid.minor = element_blank(),
      legend.position = "none" ) +
    facet_wrap(vars(rwzi), ncol = 3) +
    theme(strip.background = element_blank(),
          strip.text = element_text(color = cbPalette[5],size = 30),
          axis.title = element_text(size = 30))

ggsave(  here( runname, "figures", "manuscript", "figure2_9plants.png"),width = 16, height = 11, units = "in")
ggsave(  here( runname, "figures", "manuscript", "figure2_9plants.tiff"),width = 16, height = 11, units = "in")
ggsave(  here( runname, "figures", "manuscript", "figure2_9plants.svg"),width = 16, height = 11, units = "in")


###
# Compare 0,7,14 days before last date
###
df_posteriors %>%
  filter( date %in% (max(date) - c(0, 7, 14)  ) ) %>% 
  group_by( rwzi ) %>%
  group_split() %>%
  future_map(function(x){
    group_by(x, rwzi, date ) %>% 
    median_qi( load ) }) %>%
  bind_rows() %>%
  select(-.width, -.point, -.interval ) %>% 
  pivot_wider( names_from=date, values_from=c(load, .lower, .upper ) ) %>% 
  write_csv2(  here( runname, "output", "output_compare.csv"))

###
# Municipality level
###

# First create a tibble with the load on the municipality level
df_posteriors_municipality <- calc_df_load_municipality(df_posteriors,df_fractions)

df_posteriors_municipality %>%
  group_by(municipality) %>%
  group_split() %>%
  future_walk(function(x, y = x$municipality[1]){

    p <- ggplot(x, mapping = aes(x = date)) +
      geom_line(aes(y = load), color = cbPalette[5]) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper ), alpha = 0.25) +
      ggtitle(y) +     
      coord_cartesian(ylim = c(11, 15)) +
      scale_x_date( "Date",  breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") + 
      scale_y_continuous(name = expression(paste("Log(load) (10"^"-5"," persons"^"-1",")"))) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(color = cbPalette[5]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none")
    
    ggsave( here( runname, "figures", "municipality", str_c(y, ".png" )),
            plot = p, width = 6.5, height = 4.5, units = "in" )
    
    write_csv(x, here( runname, "output", "municipality", str_c( y, ".csv")))
  })

###
# Compare 0,7,14 days before last date, by municipality
###
df_posteriors_municipality %>% 
  filter( date %in% (max(date) - c(0, 7, 14)  ) ) %>% 
  select(-.width, -.point, -.interval ) %>% 
  pivot_wider( names_from=date, values_from=c(load, .lower, .upper ) ) %>% 
  write_csv2( here( runname, "output", "output_compare_muni.csv"))

###
# make figure for manuscript
# largest municipalities are 16, 249, 252, 292
### 

df_posteriors_municipality %>%
  left_join(df_fractions %>%
              group_by(municipality) %>%
              summarize(population = first(municipality_pop) / first(frac_municipality2RWZI)),
            by = "municipality") %>%
  filter(population %in%
           (df_fractions %>%
              group_by(municipality) %>%
              summarize(population = first(municipality_pop) / first(frac_municipality2RWZI)) %>%
              .$population %>%
              sort(decreasing = T) %>%
              .[1:9])) %T>%
  write.csv(here(runname,"output","manuscript","municipalities.csv")) %>%
  ggplot(mapping = aes(x = date, y = load, ymin = .lower, ymax = .upper)) +
  geom_line(color = cbPalette[5]) +
  geom_ribbon(alpha = 0.25) +
  coord_cartesian(ylim = c(11, 15)) +
  scale_x_date( "Date",  breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") + 
  scale_y_continuous(name = expression(paste("Log(load) (10"^"-5"," persons"^"-1",")"))) +
  theme_bw(base_size = 20) +
  theme(
        plot.title = element_text(color = cbPalette[5]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  facet_wrap(vars(municipality), ncol = 3) +
  theme(strip.background = element_blank(),
        strip.text = element_text(color = cbPalette[5],size = 30),
        axis.title = element_text(size = 30))

ggsave(here( runname, "figures", "manuscript", "figure3_9municipalities.png" ),width = 16, height = 11,units = "in")
ggsave(here( runname, "figures", "manuscript", "figure3_9municipalities.tiff" ),width = 16, height = 11,units = "in")
ggsave(here( runname, "figures", "manuscript", "figure3_9municipalities.svg" ),width = 16, height = 11,units = "in")

###
# Plots by safety region
###
df_posteriors %>%
  group_by(date) %>%
  group_split() %>%
  future_map(function(df){
    left_join(df,df_fractions, by = "rwzi") %>%
      mutate(load = frac_VR2RWZI * 10^load) %>%
      group_by(vr,date,.draw) %>%
      summarize(load = log10(sum(load)), .groups = "drop_last") %>%
      median_qi(load)
  }) %>%
  bind_rows %>%
  group_by( vr ) %>% 
  group_split() %>%
  future_walk( function(x,y = x$vr[1]){  
    p <-ggplot(x) +
      geom_ribbon(aes( x=date, ymin=.lower, ymax=.upper), alpha = 0.25) +
      geom_line(aes( x=date, y=load), color = cbPalette[5]) +
      ggtitle(y) +
      coord_cartesian(ylim = c(11, 14.5)) +
      scale_x_date("Date",  breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") + 
      scale_y_continuous(expression(paste("Log(load) (10"^"-5"," persons"^"-1",")"))) +
      theme_bw(base_size = 10) + 
      theme(
        plot.title = element_text(color = cbPalette[5]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none"
      )
    ggsave( here( runname, "figures", "safetyregion", str_c(y, ".png" )),
            plot = p, width = 6.5, height = 4.5, units = "in" )
    write_csv(x, here(runname, "output", "safetyregion", str_c(y,".csv")))})

# Clean up the enviroment by removing objects we no longer need
rm(df_sewage,df_posteriors_municipality)
invisible(gc()) # Just to be sure