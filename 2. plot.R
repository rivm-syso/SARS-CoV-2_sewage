library( tidyverse )
library( tidybayes )
library( patchwork )

source( "functions.R")
source( "settings.R")

load( viralload_filename )
lastday <- min( lastday, max(df_viralload_human_regions$Datum))
startday <- max( startday, min(df_viralload_human_regions$Datum))
rm( df_viralload_human_regions )

#Assuming file was created today
load( here( outdir_out, "model_data", str_c("posteriors_", Sys.Date(), ".rda")))

###
# Plot prob of detection
###
df_posteriors %>%
  group_by(rwzi) %>%
  slice_sample( n=100 ) %>%
  ungroup() %>% 
  mutate( sample=1:n() ) %>% 
  expand_grid( conc=seq(11, 13, by = 0.1) ) %>%
  mutate( p_det = probdetection(conc, x0,k) ) %>% 
  ggplot( ) +
  geom_line( aes( x = conc, y = p_det, group=sample ),
             alpha = .01, color = cbPalette[1] ) +
  scale_x_continuous( "Log(load)", limits = c(11, 13),
                      expand = c(0, 0),  breaks = c(11, 11.5, 12, 12.5, 13)  ) +
  scale_y_continuous( "Prob(detection)",  limits = c(0, 1),
                      expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1)  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(str_c( outdir_fig,  "prob_detection.png"), width = 6.5, height = 4.5, units = "in")

###
# Plot the Netherlands
###
df_posteriors %>%
  group_by( date, .draw ) %>% 
  summarize( weighted_concentration = sum(load * rwzi_persons) / sum( rwzi_persons ) ) %>% 
  group_by( date ) %>% 
  median_qi( weighted_concentration ) %>% 
  mutate( date = as.Date( as.character( date )) ) %>% 
  ggplot( aes( x = date, y = weighted_concentration,
               ymin = .lower, ymax = .upper ) ) +
  geom_ribbon(alpha = 0.25) +
  geom_line( color = cbPalette[1] ) +
  ggtitle("THE NETHERLANDS") +
  coord_cartesian(ylim = c(11, 14.5)) +
  scale_x_date( "Date", date_breaks = "2 month", date_labels = "%m/%y") + 
  scale_y_continuous( "Log(load)") +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(color = cbPalette[1] ),
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave( here( outdir_fig, "Netherlands", "posterior.png" ), width = 6.5, height = 4.5, units = "in")

###
# Plot each RWZI, write CSV
###
df_posteriors %>% 
  select( date, load, concentration, rwzi ) %>%
  group_by( date, rwzi, measurement=concentration ) %>% 
  median_qi( load ) %>%
  mutate( zeromeasurement = ifelse(measurement == 0, load, -2),
          date=as.Date(as.character(date)),
          rwzi=as.character(rwzi)) %>% 
  group_by( rwzi ) %>% 
  group_walk( function(x,y){
    p <- ggplot(x, mapping = aes(x = date,y = load,ymin = .lower, ymax=.upper) ) +
      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[6]) +
      geom_point(aes(y = measurement), color = cbPalette[6], size = 2.5) +
      geom_point(aes(y = zeromeasurement), color = cbPalette[6], size = 2, shape = 1) +
      coord_cartesian(ylim = c(11, 15)) +
      scale_x_date("Date", date_breaks = "2 month", date_labels = "%m/%y") + #"%b"
      scale_y_continuous("Log(load)") +
      ggtitle(y[1,1]) +
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[6]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none" )
    ggsave( here( outdir_fig, "RWZI", str_c( y[1,1], ".png")),
            plot = p, width = 6.5, height = 4.5, units = "in")
    
    write_csv(x, here( outdir_out, "RWZI", str_c( y[1,1], ".csv")))})

##
# make figure for manuscript
# largest installations are 16, 249, 252, 292
##
df_posteriors %>% 
  select( date, load, concentration, rwzi ) %>%
  filter( as.numeric(rwzi) %in% c(123, 16, 82, 69, 294, 267, 154, 187, 24)) %>% 
  group_by( date, rwzi, measurement=concentration ) %>% 
  median_qi( load ) %>%
  mutate( zeromeasurement = ifelse(measurement == 0, load, -2),
          date=as.Date(as.character(date)),
          rwzi=as.character(rwzi)) %>% 
  group_by( rwzi ) %>% 
  group_map( function(x,y){
    p <- ggplot(x, mapping = aes(x = date,y = load,ymin = .lower, ymax=.upper) ) +
      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[7]) +
      geom_point(aes(y = measurement), color = cbPalette[7], size = 2.5) +
      geom_point(aes(y = zeromeasurement), color = cbPalette[7], size = 2, shape = 1) +
      coord_cartesian(ylim = c(11, 15)) +
      scale_x_date("Date", date_breaks = "2 month", date_labels = "%m/%y") + #"%b"
      scale_y_continuous("Log(load)") +
      ggtitle(y[1,1]) +
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[7]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none" )}) %>% 
  reduce( `+` ) + # Patchwork composer of plots
  plot_layout(ncol = 3)
ggsave(  here( outdir_fig, "manuscript", "figure2_9plants.png"),width = 16, height = 11, units = "in")


###
# Compare 0,7,14 days before last date
###
df_posteriors %>%
  ungroup() %>% 
  mutate( date = as.Date( as.character( date))) %>% 
  filter( date %in% (max(date) - c(0, 7, 14)  ) ) %>% 
  group_by( rwzi, date ) %>% 
  median_qi( load ) %>% 
  select(-.width, -.point, -.interval ) %>% 
  pivot_wider( names_from=date, values_from=c(load, .lower, .upper ) ) %>% 
  write_csv2( str_c( outdir_out, "output_compare.csv"))

###
# Municipality level
###

df_vaccins <- calc_vax(startday,lastday)
df_muni <- calc_df_muni( df_posteriors, df_vaccins, startday, lastday )

df_muni %>%
  mutate( date=as.Date(date)) %>% #TODO dates as dates!
  group_by( municipality ) %>% 
  group_walk( function(x,y){
    max_h <- max(x$hospitalizations, na.rm=TRUE )
    p<- ggplot(x, mapping = aes(x = date)) +
      geom_line(aes(y = load), color = cbPalette[7]) +
      # TODO: put confidence bounds back in.
      #geom_ribbon(aes(ymin = .lower, ymax = .upper ), alpha = 0.25) +
      geom_point( aes(x = date, y = 10 + (hospitalizations*5)/max_h), color = cbPalette[7]) +
      ggtitle(y[1,1]) +
      scale_y_continuous(name = "Log(load)", sec.axis = sec_axis(~(-2*max_h +(.* max_h/5)), name = "Hospitalizations")) +
      coord_cartesian(ylim = c(11, 15)) +
      scale_x_date( "Date", date_breaks = "2 month", date_labels = "%m/%y") + #"%b"
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[7]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none")
    ggsave( here(outdir_fig, "municipality", str_c(y[1,1], ".png" )),
            plot = p, width = 6.5, height = 4.5, units = "in" )
    write_csv(x, here( outdir_out, "municipality", str_c( y[1,1], ".csv"))) } )

###
# Compare 0,7,14 days before last date, by municipality
###
df_posteriors %>% 
  select( .draw, date, load, rwzi, municipality, hospitalizations ) %>%
  left_join( df_fractions ) %>% 
  mutate( date = as.Date( as.character( date))) %>% 
  filter( date %in% (max(date) - c(0, 7, 14)  ) ) %>% 
  mutate( load_muni = frac_RWZI2municipality * 10^load ) %>% 
  group_by( date, municipality, hospitalizations, .draw ) %>%
  summarize( load = log10(sum( load_muni )), .groups="drop_last") %>% 
  group_by( municipality, date ) %>% 
  median_qi( load ) %>% 
  select(-.width, -.point, -.interval ) %>% 
  pivot_wider( names_from=date, values_from=c(load, .lower, .upper ) ) %>% 
  write_csv2( str_c( outdir_out, "output_compare_muni.csv"))

###
# make figure for manuscript
# largest municipalities are 16, 249, 252, 292
### 

df_muni %>% 
  mutate( date=as.Date(date) ) %>% 
  filter( municipality %in% c("Amsterdam","Rotterdam","s-Gravenhage",
                              "Utrecht","Eindhoven","Groningen","Tilburg","Almere","Breda" )) %>% 
  group_by( municipality ) %>% 
  group_map( function(x,y){
    max_h <- max(x$hospitalizations)
    p<- ggplot(x, mapping = aes(x = date)) +
      geom_line(aes(y = load), color = cbPalette[7]) +
      #geom_ribbon(aes(ymin = .lower, ymax = .upper ), alpha = 0.25) +
      geom_point( aes(x = date, y = 10 + (hospitalizations*5)/max_h), color = cbPalette[7]) +
      ggtitle(y[1,1]) +
      scale_y_continuous(name = "Log(load)", sec.axis = sec_axis(~(-2*max_h +(.* max_h/5)), name = "Hospitalizations")) +
      coord_cartesian(ylim = c(11, 15)) +
      scale_x_date( "Date", date_breaks = "2 month", date_labels = "%m/%y") + #"%b"
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[7]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none") }) %>% 
  reduce( `+` ) +
  plot_layout( ncol=3)
ggsave(here( outdir_fig, "manuscript",  "figure3_9municipalities.png" ),width = 16, height = 11,units = "in")

###
# Plots by safety region
# TODO: the linking to VR fractions should be precalculated!
###
df_posteriors %>%
  group_by( date, vr, .draw ) %>%
  mutate( vr_persons = sum(rwzi_persons) ) %>% 
  summarize( load = log10( sum( rwzi_persons * 10^load ) / vr_persons ), .groups="drop_last" ) %>% 
  median_qi( load ) %>%
  mutate( vr=as.character(vr),
          date=as.Date(as.character(date))) %>% 
  group_by( vr ) %>% 
  group_walk( function(x,y){  
    p <-ggplot(x) +
      geom_ribbon(aes( x=date, ymin=.lower, ymax=.upper), alpha = 0.25) +
      geom_line(aes( x=date, y=load), color = cbPalette[5]) +
      ggtitle(y[1,1]) +
      coord_cartesian(ylim = c(11, 14.5)) +
      scale_x_date("Date",date_breaks = "2 month", date_labels = "%m/%y") + 
      scale_y_continuous("Log(load)") +
      theme_bw(base_size = 20) + 
      theme(
        plot.title = element_text(color = cbPalette[5]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none"
      )
    ggsave( here( outdir_fig, "safetyregion", str_c(y[1,1], ".png" )),
            plot = p, width = 6.5, height = 4.5, units = "in" )
    write_csv(x, here(outdir_out, "safetyregion", str_c(y[1,1],".csv")))})
