# Under this subfolder, all relevant files will be written
runname <- "test1"
dir.create(here(runname), showWarnings = FALSE )

startday <- as.Date("2020-10-20")
lastday  <- as.Date("2021-10-20")

viralload_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Generated data/df_viralload_human_regions.RData"
vaccin_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Wouter/Data model/Vaccinatiegraad.csv"
hosp_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Wouter/Data model/Ziekenhuisopnames.csv"

dir.create.if.needed <- function(x){
  if( !dir.exists(x)) dir.create(x)
}

# directories
dir.create.if.needed(here(runname, "model_code_backup"))
file.copy( from = list.files( pattern = "^.*R|^.*stan"), 
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