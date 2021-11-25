startday <- as.Date("2020-10-20")
lastday  <- as.Date("2021-10-20")

viralload_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Generated data/df_viralload_human_regions.RData"
vaccin_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Wouter/Data model/Vaccinatiegraad.csv"
hosp_filename <- "/rivm/r/D114007 COVID-19/Surveillance/Modelering/Rioolwater/Wouter/Data model/Ziekenhuisopnames.csv"

# directories
outdir_fig = "./figures/"
if( !dir.exists(outdir_fig)) dir.create(outdir_fig)
if( !dir.exists(here(outdir_fig, "Netherlands" ))) dir.create(here(outdir_fig,"Netherlands"))
if( !dir.exists(here(outdir_fig, "RWZI"))) {
  dir.create(here(outdir_fig, "RWZI"))}
if( !dir.exists(here(outdir_fig, "municipality"))) {
  dir.create(here(outdir_fig, "municipality"))}
if( !dir.exists(here(outdir_fig, "municipality_hosp"))) {
  dir.create(here(outdir_fig, "municipality_hosp"))}
if( !dir.exists(here(outdir_fig, "safetyregion"))) {
  dir.create(here(outdir_fig, "safetyregion"))}

if( !dir.exists(here(outdir_fig, "manuscript"))) {
  dir.create(here(outdir_fig, "manuscript"))}
if( !dir.exists(here(outdir_fig, "Leeftijd"))) {
  dir.create(here(outdir_fig, "Leeftijd"))}


outdir_out = "./output/"
if( !dir.exists(outdir_out)) dir.create(outdir_out)
if( !dir.exists(here(outdir_out, "RWZI"))) {
  dir.create(here(outdir_out, "RWZI"))}
if( !dir.exists(here(outdir_out, "municipality"))) {
  dir.create(here(outdir_out, "municipality"))}
if( !dir.exists(here(outdir_out, "municipality_hosp"))) {
  dir.create(here(outdir_out, "municipality_hosp"))}
if( !dir.exists(here(outdir_out, "safetyregion"))) {
  dir.create(here(outdir_out, "safetyregion"))}
if( !dir.exists(here(outdir_out, "Leeftijd"))) {
  dir.create(here(outdir_out, "Leeftijd"))}

num_knots     <- 10
spline_degree <- 3
n_chains      <- 10