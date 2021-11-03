library( tidyverse )
library( tidybayes )
library( patchwork )
library( here )
library( furrr )

source( "functions.R")
source( "settings.R")

load( here( outdir_out, "model_data", str_c( "posteriors_hosp_age", Sys.Date(), ".rda")))


df_posteriors_hosp %>%
  group_by(date,municipality,hospitalizations,load) %>%
  median_qi(expected_hospitalizations,simulated_hospitalizations) %>%
  mutate(date = as.Date(as.character(date))) %>%
  group_by(municipality) %>%
  group_walk( function( x,y ){
    p <- ggplot(x, mapping = aes(x = date, y = expected_hospitalizations,
                                            ymin = simulated_hospitalizations.lower, 
                                            ymax = simulated_hospitalizations.upper)) +
      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[6]) + 
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
    ggsave( here(outdir_fig, "municipality_hosp", str_c(y[[1]], ".png")),
            plot = p, width = 6.5, height = 4.5, units = "in")
    
    write_csv(x, here( outdir_out, "municipality_hosp", str_c(y[[1]], ".csv")))})

stop("Tot hier")

# Plaatje voor Michiel met counterfactual met "verouderde" vaccinitatie data

df_plot <- df_posteriors_hosp %>% 
  group_by(municipality) %>%
  group_split() %>%
  lapply(function(x){
    y <- x %>% filter(date == "2020-09-01") %>%
      mutate(hosp_rate = expected_hospitalizations/(10^(load-19)*municipality_pop)) %>% 
      select(municipality,hosp_rate,.draw)
    x <- left_join(x,y, by = c("municipality",".draw")) %>%
      mutate(expected_hospitalizations_novax = hosp_rate*10^(load-19)*municipality_pop,
        simulated_hospitalizations_novax = rpois(nrow(.),expected_hospitalizations_novax))
  }) %>%
  bind_rows() %>%
  group_by(date,municipality,hospitalizations,load) %>%
  #slice_sample( n=100 ) %>%
  median_qi(expected_hospitalizations,simulated_hospitalizations,
            expected_hospitalizations_novax,simulated_hospitalizations_novax) %>%
  mutate(date = as.Date(as.character(date))) 

df_plot %>%
  group_by(municipality) %>%
  group_split()%>%
  lapply(function(x){
    p <- ggplot(x, mapping = aes(x = date, y = expected_hospitalizations,
                                 ymin = simulated_hospitalizations.lower, 
                                 ymax = simulated_hospitalizations.upper)) +
      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[7]) + 
      geom_point(aes(y = hospitalizations), color = cbPalette[7], size = 2.5) +
      geom_line(aes( y = expected_hospitalizations_novax), color = cbPalette[6]) +
      geom_ribbon(aes(ymin = simulated_hospitalizations_novax.lower, 
                      ymax = simulated_hospitalizations_novax.upper), alpha = 0.25) +
      scale_x_date("Date", date_breaks = "1 month", date_labels = "%m/%y") + 
      scale_y_continuous("Hospitalizations") + #,
      # sec.axis = sec_axis(trans =  ~.* conversion, 
      #                     name = "Virusvracht")) + 
      coord_cartesian(ylim = c(0, max(x$simulated_hospitalizations_novax.upper,
                                      x$hospitalizations,
                                      x$expected_hospitalizations_novax)+1)) +
      ggtitle(x$municipality[1]) +
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[6]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    ggsave( paste0( outdir_fig,"3hosp_", x$municipality[1], ".png"),
            plot = p, width = 6.5, height = 4.5, units = "in")})

df_plot %>% ungroup %>%
  filter(municipality %in% (df_plot %>% ungroup() %>%
                              select(municipality,hospitalizations) %>% 
                              group_by(municipality) %>% 
                              summarize(hospitalizations = max(hospitalizations)) %>%
                              arrange(desc(hospitalizations)) %>%
                              .$municipality %>%
                              .[1:4])) %>%
  group_by(municipality) %>%
  group_split() %>%
  lapply(function(x){
    ggplot(x, mapping = aes(x = date, y = expected_hospitalizations,
                                 ymin = simulated_hospitalizations.lower, 
                                 ymax = simulated_hospitalizations.upper)) +
      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[7]) + 
      geom_point(aes(y = hospitalizations), color = cbPalette[7], size = 2.5) +
      geom_line(aes( y = expected_hospitalizations_novax), color = cbPalette[6]) +
      geom_ribbon(aes(ymin = simulated_hospitalizations_novax.lower, 
                      ymax = simulated_hospitalizations_novax.upper), alpha = 0.25) +
      scale_x_date("Date", date_breaks = "1 month", date_labels = "%m/%y") + 
      scale_y_continuous("Hospitalizations") + #,
      # sec.axis = sec_axis(trans =  ~.* conversion, 
      #                     name = "Virusvracht")) + 
      coord_cartesian(ylim = c(0, max(x$simulated_hospitalizations_novax.upper,
                                      x$hospitalizations,
                                      x$expected_hospitalizations_novax)+1)) +
      ggtitle(x$municipality[1]) +
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[6]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )}) %>%
  cowplot::plot_grid(plotlist = .,nrow = 2) %>%
  ggsave(file = "4_municipality.jpg",width = 12, height = 7.3)
  



df_posteriors_hosp %>% 
  group_by(municipality,date) %>% 
  summarize(hosp_per_million = first(hospitalizations)/first(municipality_pop)*10^6,
            load = first(load), 
            percentage_vax = first(percentage_vax),
            .groups = "drop_last") %>%
  mutate(vaccinatiegraad_2_weken_terug = lag(percentage_vax,14,default = 0),
         date = as.Date(date),
         delta_dominant = if_else(date > as.Date("2021-06-28"),"Y","N")) %>%
  group_split() %>%
  lapply(function(x){
    p <- ggplot(x) + geom_point(aes(x = load, y = hosp_per_million, 
                                    color = vaccinatiegraad_2_weken_terug,
                                    shape = delta_dominant)) +
      scale_colour_gradientn(colours=rainbow(6)) +
      ylab("Ziekenhuisopnames per 1 miljoen inwoners") +
      xlab("log load (modeluitkomst)") +
      ggtitle(x$municipality[1])
    
    ggsave( paste0( outdir_fig,"Virusvracht_ziekenhuis_vaccin/", x$municipality[1], ".png"),
            plot = p, width = 6.5, height = 4.5, units = "in")
    
  })


# Include virusvracht
df_posteriors_hosp %>%
  group_by(date,municipality,hospitalizations,load) %>%
  #slice_sample( n=100 ) %>%
  median_qi(expected_hospitalizations,simulated_hospitalizations) %>%
  mutate(date = as.Date(as.character(date))) %>%
  group_by(municipality) %>%
  group_split()%>%
  lapply(function(x,conversion){
    conversion <- filter(conversion, municipality == x$municipality[1]) %>% .$hosp_rate
    p <- ggplot(x, mapping = aes(x = date, y = expected_hospitalizations,
                                 ymin = simulated_hospitalizations.lower, 
                                 ymax = simulated_hospitalizations.upper)) +
      geom_ribbon(alpha = 0.25) +
      geom_line(color = cbPalette[7]) + 
      geom_point(aes(y = hospitalizations), color = cbPalette[7], size = 2.5) +
      geom_line(aes( y = 10^(load-19)*conversion), color = cbPalette[6]) +
      scale_x_date("Date", date_breaks = "1 month", date_labels = "%m/%y") + 
      scale_y_continuous("Hospitalizations") + #,
                        # sec.axis = sec_axis(trans =  ~.* conversion, 
                        #                     name = "Virusvracht")) + 
      coord_cartesian(ylim = c(0, max(x$simulated_hospitalizations.upper,
                                      x$hospitalizations,
                                      x$expected_hospitalizations)+1)) +
      ggtitle(x$municipality[1]) +
      theme_bw(base_size = 20) +
      theme(
        plot.title = element_text(color = cbPalette[6]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    ggsave( paste0( outdir_fig,"1hosp_", x$municipality[1], ".png"),
            plot = p, width = 6.5, height = 4.5, units = "in")},C)





>>>>>>> ee0b8dcec2ee0061a9ea5142f77231901534bef4

# make figure for manuscript
color = cbPalette[7]
plots <- list()
mun_number <- c(16, 249, 252, 292, 91, 111, 283, 9, 49)
for (m in 1:9) {
  municipality <- mun_number[[m]]
  conc_med <- array(NA, numberofdays)
  conc_ub <- array(NA, numberofdays)
  conc_lb <- array(NA, numberofdays)
  cases <- array(NA, numberofdays)
  for (j in 1:numberofdays) {
    conc_med[j] <-mean(params_hospitalization$expected_hospitalizations[, j, municipality])
    # conc_med[j] <-median(params_hospitalization$simulated_hospitalizations[, j, municipality])
    conc_lb[j] <-
      quantile(params_hospitalization$simulated_hospitalizations[, j, municipality], probs = 0.025)
    conc_ub[j] <-
      quantile(params_hospitalization$simulated_hospitalizations[, j, municipality], probs = 0.975)
    cases[j] <- hospitaldata_matrix[municipality, j]
  }
  df_conc <- data.frame(
    list(
      date = as.Date(1:numberofdays, origin = "2020-09-01"), #startday
      median_conc = conc_med,
      lower_conc = conc_lb,
      upper_conc = conc_ub
    )
  )
  df_hosp <- data.frame(
    list(
      date = as.Date(1:numberofdays, origin = "2020-09-01"),
      cases = cases
    )
  )
  max_y <- max(df_hosp$cases)
  plots[[m]] <- 
    ggplot(
      data = df_conc,
      mapping = aes(x = date)) +
    geom_point(data = df_hosp, aes(x = date, y = cases, color = cbPalette[7])) +
    geom_line(aes(y = median_conc), color = color) +
    geom_ribbon(aes(ymin = lower_conc, ymax = upper_conc), alpha = 0.25) +
    labs(x = "Date" , y = "Hospitalizations") +
    ggtitle(colnames(df_RWZItoMUN)[[1+municipality]]) +
    scale_y_continuous(name = "Hospitalizations") + # libary(scales) , breaks= pretty_breaks()
    coord_cartesian(ylim = c(0, max_y+2)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%m/%y") + #"%b"
    theme_bw(base_size = 20) +
    scale_color_manual(values = color) + # meh, improve
    theme(
      plot.title = element_text(color = color),
      panel.grid.major = element_line(size = 0.7),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

plots[1:6] <- plots[1:6] %>% 
  map(~ .x + theme(axis.title.x = element_blank())) 
plots[1:2] <- plots[1:2] %>% 
  map(~ .x + theme(axis.title.y.right = element_blank())) 
plots[4:5] <- plots[4:5] %>% 
  map(~ .x + theme(axis.title.y.right = element_blank())) 
plots[7:8] <- plots[7:8] %>% 
  map(~ .x + theme(axis.title.y.right = element_blank())) 
plots[2:3] <- plots[2:3] %>% 
  map(~ .x + theme(axis.title.y.left = element_blank())) 
plots[5:6] <- plots[5:6] %>% 
  map(~ .x + theme(axis.title.y.left = element_blank())) 
plots[8:9] <- plots[8:9] %>% 
  map(~ .x + theme(axis.title.y.left = element_blank())) 

figure4 <- 
  grid.arrange(
    plots[[1]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[2]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[3]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[4]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[5]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[6]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[7]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[8]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    plots[[9]] + theme(plot.margin=unit(c(5.5,5,5.5,5.5), "pt")),
    nrow = 3
  )
ggsave(
  paste0(
    outdir_fig, 
    "figure4_9municipalities.png"
  ),
  plot = figure4, 
  width = 16, 
  height = 11, 
  units = "in"
)


# plot the estimated (hyper) distribution of hospitalizations at log(load)=13
# for calculatation of sewage risk levels that correspond to existing risk levels for hospitalizations
pow = 0.0 # e.g., -0.5 log10 step lower than 13, ie 12.5
posterior_thresholds = tibble(sample = 1:n_samples,
                              mean = 10^pow * params_hospitalization$mean_hosp_rate,
                              rate = params_hospitalization$sigma_hosp_rate,
                              p1 = pgamma(4/7, shape = mean*rate, rate = rate),
                              p2 = pgamma(16/7, shape = mean*rate, rate = rate),
                              p3 = pgamma(27/7, shape = mean*rate, rate = rate))

mp1 <- median(posterior_thresholds$p1)
mp2 <- median(posterior_thresholds$p2)
mp3 <- median(posterior_thresholds$p3)
mp1
mp2-mp1
mp3-mp2
1-mp3

posterior_parameters = tibble(sample = 1:n_samples,
                              mean = 10^pow * params_hospitalization$mean_hosp_rate, 
                              rate = params_hospitalization$sigma_hosp_rate)

plot_data =
  pmap_df(posterior_parameters,
          function(sample, mean, rate) {
            tibble(
              sample = sample,
              hosps = seq(0,10, by = 0.1),
              prob_hosps = dgamma(hosps, shape = mean* rate, rate = rate)
            )
          })

plot_hospitalization = ggplot(data = plot_data) +
  geom_rect(aes(xmin=0, xmax=4/7, ymin=0, ymax=Inf), fill = "green", alpha = 0.005) +
  geom_rect(aes(xmin=4/7, xmax=16/7, ymin=0, ymax=Inf), fill = "yellow", alpha = 0.005) +
  geom_rect(aes(xmin=16/7, xmax=27/7, ymin=0, ymax=Inf), fill = "orangered1", alpha = 0.005) +
  geom_rect(aes(xmin=27/7, xmax=8, ymin=0, ymax=Inf), fill = "brown", alpha = 0.005) +
  geom_line(
    aes(group = sample, x = hosps, y = prob_hosps),
    alpha = .02,
    color = cbPalette[2]
  ) +
  labs(x = "Hospitalizations", y = "Density") +
  ggtitle("Expected daily hospitalizations\nper million at log(load) = 13\nper 100,000 per day") +
  #geom_vline(xintercept = 4/7, linetype = "dotted", alpha = 0.5) +
  #geom_vline(xintercept = 16/7) +
  #geom_vline(xintercept = 27/7) +
  scale_x_continuous(
    limits = c(0, 8),
    expand = c(0, 0),
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8)
  ) +
  scale_y_continuous(
    limits = c(0, 0.70),
    expand = c(0, 0),
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

plot_hospitalization

ggsave(
  paste0(
    "./Michiel/corona-wastewater/figures/png/load132.png"
  ),
  plot = plot_hospitalization, 
  width = 6.5, 
  height = 4.5, 
  units = "in"
)

# make figure for manuscript
hist(params_hospitalization$hosp_rate[,mun_number[1]])

posterior_parameters = tibble(sample = 1:n_samples,
                              mean = 10^pow * params_hospitalization$mean_hosp_rate, 
                              rate = params_hospitalization$sigma_hosp_rate)

plot_data =
  pmap_df(posterior_parameters,
          function(sample, mean, rate) {
            tibble(
              sample = sample,
              hosps = seq(0,10, by = 0.1),
              prob_hosps = dgamma(hosps, shape = mean* rate, rate = rate)
            )
          })

posterior_rates_municipalities = tibble(sample = 1:n_samples,
                                        Amsterdam = params_hospitalization$hosp_rate[,mun_number[1]],
                                        Rotterdam = params_hospitalization$hosp_rate[,mun_number[2]],
                                        "The Hague" = params_hospitalization$hosp_rate[,mun_number[3]], 
                                        Utrecht = params_hospitalization$hosp_rate[,mun_number[4]])
posterior_rates_municipalities_long <- posterior_rates_municipalities %>%
  gather(sample) %>%
  rename(City = sample) %>%
  rename(rate = value)

figure5 = ggplot(data = plot_data) +
  geom_line(
    aes(group = sample, x = hosps, y = prob_hosps),
    alpha = .02,
    color = color
  ) +
  geom_histogram(
    data = posterior_rates_municipalities_long,
    aes(x = rate, color = City, fill = City, after_stat(density*0.1)),
    binwidth = 0.1,
    alpha = 0.1
  ) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") + 
  labs(x = expression(paste("Hospitalization rate (",10^-6,"",  day^-1, ")")), y = "Density/Probability") +
  #ggtitle("Expected daily hospitalizations\nper million at log(load) = 13\nper 100,000 per day") +
  scale_x_continuous(
    limits = c(0, 8),
    expand = c(0, 0),
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8)
  ) +
  scale_y_continuous(
    limits = c(0, 0.5),
    expand = c(0, 0),
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

figure5

ggsave(
  paste0(
    outdir_fig, 
    "figure5_9municipalities.png"
  ),
  plot = figure5, 
  width = 8, 
  height = 5, 
  units = "in"
)
