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

load_if_needed( "df_posteriors_hosp", here( runname, "output", "model_data", "posteriors_hosp.RData") )
load_if_needed( "^df_muni", here( runname, "output", "model_data", "df_muni.RData")) # Note ^ in regular expression to avoid finding "calc_df_muni"

df_plot_hosp <- df_muni %>% 
  group_by(municipality) %>% 
  group_split() %>%
  future_map(function(x){
    left_join(x,df_posteriors_hosp, by = c("municipality", "age_group")) %>%
      mutate(expected_hospitalizations_cf = hosp_rate * 10^(load - ref_load) * population,
             expected_hospitalizations = expected_hospitalizations_cf *
               (1 - prevention_vax*percentage_vax),
             simulated_hospitalizations = rpois(nrow(.),expected_hospitalizations),
             simulated_hospitalizations_cf = rpois(nrow(.),expected_hospitalizations_cf)) %>%
      bind_rows(group_by(.,date,municipality,.chain,.iteration,.draw) %>%
                  summarize(age_group = "Totaal",
                            across(contains("hospitalizations"),sum),
                            .groups = "drop")) %>%
      group_by(date,municipality,hospitalizations,age_group) %>%
      median_qi(expected_hospitalizations,simulated_hospitalizations,
                expected_hospitalizations_cf,simulated_hospitalizations_cf) %>%
      mutate(date = as.Date(as.character(date)))
  },.options = furrr_options(seed = T)) %>%
  bind_rows() %>%
  group_by(municipality) %>%
  group_split()

save(df_plot_hosp, file = here( runname, "output" , "model_data", "df_plot.RData"))

##### We plot the fitted hospitalizations for each municipality #####
df_plot_hosp %>% future_map(function(x){
    p <- list()
    for(i in sort(unique(x$age_group))){
      y <- filter(x, age_group == i)
      p[[i]] <- ggplot(y, mapping = aes(x = date, y = expected_hospitalizations,
                                        ymin = simulated_hospitalizations.lower, 
                                        ymax = simulated_hospitalizations.upper)) +
        geom_ribbon(alpha = 0.25) +
        geom_line(color = cbPalette[6]) + 
        geom_point(aes(y = hospitalizations), color = cbPalette[6], size = 2.5) +
        scale_x_date("Date", date_breaks = "1 month", date_labels = "%m/%y") + 
        scale_y_continuous("Hospitalizations") + 
        coord_cartesian(ylim = c(0, max(y$simulated_hospitalizations.upper,
                                        y$hospitalizations,
                                        y$expected_hospitalizations)+1)) +
        ggtitle(i) +
        theme_bw(base_size = 20) +
        theme(
          plot.title = element_text(color = cbPalette[6]),
          panel.grid.major = element_line(size = 0.7),
          panel.grid.minor = element_blank(),
          legend.position = "none"
        )
    }
    
    ggsave( here(runname,"figures", "Leeftijd", str_c("hosp_", x$municipality[1], ".png")),
            plot = wrap_plots(p, ncol = 2),
            width = 18, height = 4*ceiling(length(p)/2), units = "in")
    
    write_csv(x,  here(runname, "output", "Leeftijd", str_c("hosp_", x$municipality[1], ".csv")))
})

#### We include the counterfactuals in the case we did not have vaccinations ####

df_plot_hosp %>% future_map(function(x){
  p <- list()
  for(i in sort(unique(x$age_group))){
    y <- filter(x, age_group == i)
    p[[i]] <- ggplot(y, mapping = aes(x = date, y = expected_hospitalizations,
                                      ymin = simulated_hospitalizations.lower, 
                                      ymax = simulated_hospitalizations.upper)) +

      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[6]) + 
      geom_ribbon(alpha = 0.25, mapping = aes(ymin = simulated_hospitalizations_cf.lower,
                                              ymax = simulated_hospitalizations_cf.upper)) +
      geom_line(color = cbPalette[7], mapping = aes(y = expected_hospitalizations_cf)) + 
      geom_point(aes(y = hospitalizations), color = cbPalette[6], size = 2.5) +
      scale_x_date("Date", date_breaks = "1 month", date_labels = "%m/%y") + 
      scale_y_continuous("Hospitalizations") + 

      coord_cartesian(ylim = c(0, max(x$simulated_hospitalizations.upper,
                                      x$hospitalizations,
                                      x$expected_hospitalizations)+1)) +
      ggtitle(y[[1]]) +
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[6]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none"
      )
  }
  
  ggsave( here(runname, "figures", "Leeftijd", str_c("cf_hosp_", x$municipality[1], ".png")),
          plot = wrap_plots(p, ncol = 2),
          width = 18, height = 4*ceiling(length(p)/2), units = "in")
  
  write_csv(x, here( runname, "output", "Leeftijd", str_c("cf_hosp_",x$municipality[1], ".csv")))
})
