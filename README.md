# Modelling patterns of SARS-CoV-2 circulation in the Netherlands, August 2020-February 2022, revealed by a nationwide sewage surveillance program

* [Michiel van Boven](mailto://michiel.van.boven@rivm.nl), Wouter A Hetebrij, Arno N Swart, Erwin Nagelkerke, Rudolf FHJ van der Beek, Sjors Stouten, Rudolf T Hoogeveen, Fuminari Miura, Astrid Kloosterman, Anne-Merel van der Drift, Anne Welling, Willemijn J Lodder, Ana M de Roda-Husman*

Under review for Epidemics.

---

## Background

Surveillance of SARS-CoV-2 in wastewater offers an unbiased and near real-time tool to track circulation of SARS-CoV-2 at a local scale, next to other epidemic indicators such as hospital admissions and test data. However, individual measurements of SARS-CoV-2 in sewage are noisy, inherently variable, and can be left-censored.

## Aim

We aimed to infer latent virus loads in a comprehensive sewage surveillance program that includes all sewage treatment plants (STPs) in the Netherlands and covers 99.6% of the Dutch population.

## Methods

A multilevel Bayesian penalized spline model was developed and applied to estimate time- and STP-specific virus loads based on water flow adjusted SARS-CoV-2 qRTPCR data from 1‑4 sewage samples per week for each of the >300 STPs.

## Results

The model provided an adequate fit to the data and captured the epidemic upsurges and downturns in the Netherlands, despite substantial day-to-day measurement variation. Estimated STP virus loads varied by more than two orders of magnitude, from approximately 10 12 (virus particles per 100,000 persons per day) in the epidemic trough in August 2020 to almost 10 15 in many STPs in January 2022. Epidemics at the local levels were slightly shifted between STPs and municipalities, which resulted in less pronounced peaks and troughs at the national level.

## Conclusion

Although substantial day-to-day variation is observed in virus load measurements, wastewater-based surveillance of SARS-CoV-2 can track long-term epidemic progression at a local scale in near real-time, especially at high sampling frequency.

---

## Installing and running the code

Clone the repository into a local directory, open the project file `wastewater.Rproj` and run the R files in the specified order: 
1. `source("0. functions.R")`, defines needed functions
2. `source("0. settings.R")`, global settings
3. `source("1. wastewater.R")`, runs the model and stores the output in `./results/`
4. `source("2. plot.R")`, creates plots and stores them in `./figures/`

Make sure to have all needed libraries installed. In particular this code relies on [stan](https://mc-stan.org/)
