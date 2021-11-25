library( tidyverse )
library( tidybayes )
library( patchwork )
library( here )
library( furrr )

source( "functions.R")
source( "settings.R")

load( here( outdir_out, "model_data", str_c( "posteriors_hosp_age", Sys.Date(), ".rda")))
load( here( outdir_out, "model_data", str_c("df_muni_age", Sys.Date(),".RData")))

df_plot_hosp <- df_muni %>% 
  group_by(municipality) %>% 
  group_split() %>%
  future_map(function(x){
    left_join(x,df_posteriors_hosp, by = c("municipality", "age_group")) %>%
      group_by(age_group) %>%
      group_split() %>%
      lapply(function(df){
        arrange(df,date) %>%
          mutate(percentage_vax = lag(percentage_vax, n = 14, default = 0),
                 expected_hospitalizations_cf = hosp_rate * 10^(load-19) * population,
                 expected_hospitalizations = expected_hospitalizations_cf *
                   (1 - prevention_vax*percentage_vax),
                 simulated_hospitalizations = rpois(nrow(.),expected_hospitalizations),
                 simulated_hospitalizations_cf = rpois(nrow(.),expected_hospitalizations_cf))}) %>%
      bind_rows() %>%
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

save(df_plot_hosp, file = here( outdir_out, "model_data", str_c("df_plot", Sys.Date(),".RData")))

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
    
    ggsave( paste0( outdir_fig,"Leeftijd/hosp_", x$municipality[1], ".png"),
            plot = wrap_plots(p, ncol = 2),
            width = 18, height = 4*ceiling(length(p)/2), units = "in")
    
    write_csv(x, paste0( outdir_out, "Leeftijd/hosp_",x$municipality[1], ".csv"))
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
  
  ggsave( paste0( outdir_fig,"Leeftijd/cf_hosp_", x$municipality[1], ".png"),
          plot = wrap_plots(p, ncol = 2),
          width = 18, height = 4*ceiling(length(p)/2), units = "in")
  
  write_csv(x, paste0( outdir_out, "Leeftijd/cf_hosp_",x$municipality[1], ".csv"))
})