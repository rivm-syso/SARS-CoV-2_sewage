# Needed for finding the most recent datafiles from EPI
library(tidyverse)
library(here)
library(furrr)

# Under this subfolder, all relevant files will be written
runname <- "results/paper_augustus_2020_february_2022_15_knots"
dir.create(here(runname), showWarnings = FALSE, recursive = TRUE ) # allows all runs to be in a single subdirectory if runname consists of stacked folders

startday <- as.Date("2020-08-01")
lastday  <- as.Date("2022-02-08")

# We need the most recent, i.e. the highest date, file. Which means that alphabetically 
# it is the last file 
vaccin_mapname <- "/rivm/r/COVID-19/Toegang_externen/ZO_Rioolwater/EPI_Tabellen/"

vaccin_filename <-  list.files(vaccin_mapname,pattern = "vaccinatiegraad.+_gemeente.+\\.csv$",
                               full.names = T) %>% sort(decreasing = T) %>% .[1]
hosp_filename <- list.files(vaccin_mapname,pattern = "ziekenhuisopnames_gemeenten.+\\.csv$",
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

# Finally, we copy df_viralload to the newly created map if it doesn't exists yet
viralload_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Generated data/df_viralload_human_regions.RData"
if(!file.exists(here(runname,"output", "model_data","df_viralload_human_regions.RData"))){
  file.copy(from = viralload_filename,
            to = here(runname,"output", "model_data","df_viralload_human_regions.RData"))
}
viralload_filename <- here(runname,"output", "model_data","df_viralload_human_regions.RData")

# Settings for wastewater-model
num_knots     <- 15
spline_degree <- 3
n_chains      <- 10

# Settings for hospitalizations-model
max_delay     <- 0
ref_load      <- 19
delay_vax     <- 14

# Number of cores for parallelization
plan(multisession, workers = 10)

settings_sourced <- T