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
  int<lower = 1> n_municipality;                                      // number of regions (province, safety region, municipality)
  int<lower = 1> n_date;                                              // number of observation days 
  int<lower = 1> n_age_group;
  
  int<lower=1, upper=n_municipality> municipality[n];
  int date[n];
  int age_group[n];
  
  vector<lower=0, upper = 100>[n] age;
  vector<lower=1>[n] population;
  vector<lower=0>[n] load;
  vector<lower=0, upper = 1>[n] percentage_vax;
  
  int hospitalizations[n];                                            // now included at the level of sewage plants

  int delay_vax;                                                      // delay between vaccination and effectiveness
  int max_delay;                                                      // maximum shift between sewage data and hospitalizations
  int ref_load;                                                       // reference sewage load for hospitalization rates
}

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
  real<lower = 0> mean_hosp_rate;                                      // hospitalization rates with random effect
  real<lower = 0> sigma_hosp_rate;                                     // hospitalization rates with random effect
  vector<lower = 0>[n_municipality] hosp_rate;                 // hospitalization rates with random effect 
  simplex[2] hosp_rate_age;
  vector<lower = 0, upper = 1>[n_age_group] prevention_vax;       // effect of vaccination on hospitalizations

}

model {
  vector[n] sum_load = exp( log(10) *(load - ref_load));
  vector[n] hosp_parameter;
  /* (hyper)parameter for hospitalization rates */
  hosp_rate ~ gamma(sigma_hosp_rate * mean_hosp_rate, sigma_hosp_rate); 
  
  /* hyper parameter effect of vaccination */
  prevention_vax ~ uniform(0,1); // beta-distribution would be better
  /* log-likelihood contributions */
  
   for ( i in 1 : n){
     hosp_parameter[i] = hosp_rate[municipality[i]] * 
                        (hosp_rate_age[1]*log(age[i]) + hosp_rate_age[2]) * 
                        sum_load[i] *
                        (1 - prevention_vax[age_group[i]]*percentage_vax_delay[i]) *
                        population[i];
  }
  
  hospitalizations ~ poisson(hosp_parameter);

}

generated quantities {
  vector[n] sum_load = exp( log(10) *(load - ref_load));
  vector[n] hosp_parameter;

  matrix[n_date, n_municipality] expected_hospitalizations = rep_matrix(0,n_date, n_municipality); // for calculation of credible intervals
  int simulated_hospitalizations[n_date, n_municipality] = rep_array(0,n_date, n_municipality); // for calculation of prediction intervals

  for ( i in 1 : n){
     hosp_parameter[i] = hosp_rate[municipality[i]] *
                        (hosp_rate_age[1]*log(age[i]) + hosp_rate_age[2]) *
                        sum_load[i] *
                        (1 - prevention_vax[age_group[i]]*percentage_vax_delay[i]) *
                        population[i];

    expected_hospitalizations[date[i],municipality[i]] = expected_hospitalizations[date[i],municipality[i]] + hosp_parameter[i];
    simulated_hospitalizations[date[i],municipality[i]] = simulated_hospitalizations[date[i],municipality[i]] + poisson_rng(hosp_parameter[i]);
  }

}

// generated quantities {
//   vector[n] sum_load = exp( log(10) *(load - ref_load));
//   vector[n] hosp_parameter;
//   vector[n] hosp_parameter_novax;
// 
//   matrix[n_date, n_municipality] expected_hospitalizations = rep_matrix(0,n_date, n_municipality); // for calculation of credible intervals
//   int simulated_hospitalizations[n_date, n_municipality] = rep_array(0,n_date, n_municipality); // for calculation of prediction intervals
//   
//   matrix[n_date, n_municipality] expected_hospitalizations_novax = rep_matrix(0,n_date, n_municipality); // for calculation of credible intervals
//   int simulated_hospitalizations_novax[n_date, n_municipality] = rep_array(0,n_date, n_municipality); // for calculation of prediction intervals
// 
//   for ( i in 1 : n){
//      hosp_parameter[i] = hosp_rate[municipality[i]] * 
//                         (hosp_rate_age[1]*log(age[i]) + hosp_rate_age[2]) * 
//                         sum_load[i] *
//                         (1 - prevention_vax*percentage_vax_delay[i]) *
//                         population[i];
//                         
//     expected_hospitalizations[date[i],municipality[i]] = expected_hospitalizations[date[i],municipality[i]] + hosp_parameter[i];
//     simulated_hospitalizations[date[i],municipality[i]] = simulated_hospitalizations[date[i],municipality[i]] + poisson_rng(hosp_parameter[i]);
//   
//     hosp_parameter_novax[i] = hosp_rate[municipality[i]] * 
//                         (hosp_rate_age[1]*log(age[i]) + hosp_rate_age[2]) * 
//                         sum_load[i] *
//                         population[i];
//                         
//     expected_hospitalizations_novax[date[i],municipality[i]] = expected_hospitalizations_novax[date[i],municipality[i]] + hosp_parameter_novax[i];
//     simulated_hospitalizations_novax[date[i],municipality[i]] = simulated_hospitalizations_novax[date[i],municipality[i]] + poisson_rng(hosp_parameter_novax[i]);
//   }
// 
// }
