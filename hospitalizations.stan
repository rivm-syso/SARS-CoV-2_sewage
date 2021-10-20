/* Stan program for inference using of time- and location-varying        */
/* load of SARS-CoV-2 RNA in Dutch wastewater plants.                    */ 
/* Global load levels can be replaced by global and regional             */
/* loads.                                                                */                       
/* (c) Michiel van Boven, licensed with BSD3 clause (use but mention)    */
/* Created: October 2020                                                 */
/* Updated: January 2021 to integrate with safety region level           */
/* analysis and hospitalization data by admission date and safety region */

/* TODO: p-spline at the country level; hospitalization data by admission date,       */
/* replace round function with interpolation function, mumps-like variance components,*/
/* better: try shifted spline instead of interpolation function                       */
/* update 2/2021: two-step approach: first fit sewage data, then combine with hosps   */
/* update 3/2021: cleanup and choice for shifted spline - slow but hopefully reliable */
/* update 6/3/2021: greatly simplified and cleanup: weighted past as predictor        */
/* using truncated geometric distributions                                            */

/* TODO: max_delay just throws away the first few days instead of modeling a shift 
 between hospitalizations and viral loads. */


data {
  int<lower=1> n;
  int<lower = 1> n_municipality;  // number of regions (province, safety region, municipality)
  int<lower = 1> n_date;          // number of observation days 
  int<lower = 1> n_age_group;     // number of different age-groups                             

  /* Classifications of the data */
  int<lower=1, upper=n_municipality> municipality[n];
  int<lower=1, upper=n_age_group> age_group[n];
  int<lower=1, upper=n_date> date[n];

  /* For vector-based calculations, we work exclusively with vectors */
  vector<lower=1>[n] population;
  vector<lower=0>[n] load;
  vector<lower=0,upper=1>[n] percentage_vax;
  
  int hospitalizations[n];        // now included at the level of sewage plants

  int ref_load;                   // reference sewage load for hospitalization rates
  int delay_vax;                  // Delay between administration of vaccin and effectiveness
                                  // Perhaps make this another parameter
}

/* The vaccins have a delay between administration and effectiveness, hence we 
shift the vaccination coverage in the data */
transformed data {
  real<lower=0, upper = 1> percentage_vax_mat[n_date,n_municipality,n_age_group];
  vector<lower = 0, upper = 1>[n] percentage_vax_delay;

  for( i in 1:n ){
      percentage_vax_mat[date[i],municipality[i],age_group[i]] = percentage_vax[i];
  }
  
  // Shift the vaccination data by delay_vax
  for(i in 1:n){
    if(date[i] - delay_vax <= 0){
      percentage_vax_delay[i] = 0;
    } else {
      percentage_vax_delay[i] = percentage_vax_mat[date[i]-delay_vax,municipality[i],age_group[i]];
    }
  }
  
}

parameters {
  /* hospitalization rates and hospitalization delay */
  vector<lower = 0>[n_age_group] mean_hosp_rate;                 // hospitalization rates with random effect 
  vector<lower = 0>[n_age_group] sigma_hosp_rate;
  matrix<lower = 0>[n_age_group,n_municipality] hosp_rate;       // Can be slightly improved by making array of rowvectors
  
  vector<lower=0,upper=1>[n_age_group] prevention_vax;
}


model {
  vector[n] sum_load = exp( log(10) *(load - ref_load));
  vector[n] hosp_parameter;
  
  /* For a slight increase in speed, we want to use vector-valued calculations */
  vector[n] hosp_rate_vector;
  vector[n] prevention_vax_vector;

  /* Prior on the effectiveness of vaccination */
  prevention_vax ~ beta(27,3);

  for(i in 1:n){
    hosp_rate_vector[i] = hosp_rate[age_group[i],municipality[i]];
    prevention_vax_vector[i] = prevention_vax[age_group[i]];
  }

  /* The hospitalizations ratios for each municipality is governed by the overall
      hosp_rate for the concerning municipality */
 for(i in 1:n_age_group){
    hosp_rate[i,] ~ gamma(sigma_hosp_rate[i] * mean_hosp_rate[i], sigma_hosp_rate[i]);
 }

  hosp_parameter = hosp_rate_vector .*
                    (1 - prevention_vax_vector .* percentage_vax_delay) .*
                    sum_load .*
                    population;
  
  hospitalizations ~ poisson(hosp_parameter);

}

// generated quantities {
//   vector[n] sum_load = exp( log(10) *(load - ref_load));
//   real expected_hospitalizations[n_date,n_municipality,n_age_group];
//   real expected_hospitalizations_cf[n_date,n_municipality,n_age_group];
// 
//   int simulated_hospitalizations[n_date,n_municipality,n_age_group];
//   int simulated_hospitalizations_cf[n_date,n_municipality,n_age_group];
// 
//   for ( i in 1 : n){
//       expected_hospitalizations[date[i],municipality[i],age_group[i]] =
//                         hosp_rate[age_group[i],municipality[i]] *
//                         (1 - prevention_vax[age_group[i]]*percentage_vax_delay[i]) *
//                         sum_load[i] *
//                         population[i];
//                         
//       simulated_hospitalizations[date[i],municipality[i],age_group[i]] =
//                         poisson_rng(expected_hospitalizations[date[i],municipality[i],age_group[i]]);
// 
//       expected_hospitalizations_cf[date[i],municipality[i],age_group[i]] =
//                         hosp_rate[age_group[i],municipality[i]] *
//                         sum_load[i] *
//                         population[i];
//                         
//       simulated_hospitalizations_cf[date[i],municipality[i],age_group[i]] =
//                         poisson_rng(expected_hospitalizations_cf[date[i],municipality[i],age_group[i]]);
//   }
// 
// }
