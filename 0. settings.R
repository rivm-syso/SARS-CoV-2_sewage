# Needed for finding the most recent datafiles from EPI
library(tidyverse)
library(here)

# Under this subfolder, all relevant files will be written
runname <- "runs/test1"
dir.create(here(runname), showWarnings = FALSE, recursive = TRUE ) # allows all runs to be in a single subdirectory if runname consists of stacked folders

startday <- as.Date("2020-09-01")
lastday  <- as.Date("2021-11-30")

viralload_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Generated data/df_viralload_human_regions.RData"
vaccin_mapname <- "/rivm/r/COVID-19/Toegang_externen/ZO_Rioolwater/EPI_Tabellen/"

# We need the most recent, i.e. the highest date, file. Which means that alphabetically 
# it is the last file 
vaccin_filename <-  list.files(vaccin_mapname,pattern = "vaccinatiegraad_leeftijdsgroep_gemeente.+csv$",
                               full.names = T) %>% sort(decreasing = T) %>% .[1]
hosp_filename <- list.files(vaccin_mapname,pattern = "ziekenhuisopnames_gemeenten_leeftijdsgroepen.+csv$",
                            full.names = T) %>% sort(decreasing = T) %>% .[1]

dir.create.if.needed <- function(x){
  if( !dir.exists(x)) dir.create(x)
}

# directories
dir.create.if.needed(here(runname, "model_code_backup"))
file.copy( from = list.files( pattern = "^[0-9].*R|^.*stan"), 
           to   = here(runname, "model_code_backup"))

walk( c("figures", "output"),
      function(x){
        outdir <- here(runname, x)
        dir.create.if.needed(outdir)
        dir.create.if.needed(here(outdir, "model_data" ))
        dir.create.if.needed(here(outdir, "Netherlands" ))
        dir.create.if.needed(here(outdir, "RWZI" ))
        dir.create.if.needed(here(outdir, "municipality" ))
        dir.create.if.needed(here(outdir, "municipality_hosp" ))
        dir.create.if.needed(here(outdir, "safetyregion" ))
        dir.create.if.needed(here(outdir, "manuscript" ))
        dir.create.if.needed(here(outdir, "Leeftijd" ))})

num_knots     <- 10
spline_degree <- 3
n_chains      <- 10

settings_sourced <- T