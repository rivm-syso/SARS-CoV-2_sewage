# Needed for finding the most recent datafiles from EPI
library(tidyverse)
library(here)
library(furrr)

# Under this subfolder, all relevant files will be written
runname <- "results/paper"
dir.create(here(runname), showWarnings = FALSE, recursive = TRUE ) # allows all runs to be in a single subdirectory if runname consists of stacked folders

startday <- as.Date("2020-08-01")
lastday  <- as.Date("2022-02-08")

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
        dir.create.if.needed(here(outdir, "safetyregion" ))
        dir.create.if.needed(here(outdir, "manuscript" ))})

# Finally, we copy df_viralload to the newly created map if it doesn't exists yet
viralload_filename <- "df_viralload_human_regions.RData"
if(!file.exists(here(runname,"output", "model_data","df_viralload_human_regions.RData"))){
  file.copy(from = viralload_filename,
            to = here(runname,"output", "model_data","df_viralload_human_regions.RData"))
}
viralload_filename <- here(runname,"output", "model_data","df_viralload_human_regions.RData")

# Settings for wastewater-model
num_knots     <- 15
spline_degree <- 3
n_chains      <- 10

# Number of cores for parallelization
plan(multisession, workers = 10)

settings_sourced <- T