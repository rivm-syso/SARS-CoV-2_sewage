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
  int<lower = 0> n_age_group;                                         // number of different age groups
  int<lower=1, upper=n_municipality> municipality[n];
  int<lower=1, upper = n_age_group> age_group[n];
  int<lower=1> age_population[n];
  int date[n];
  real<lower=0> load[n];
  real<lower=0, upper = 1> percentage_vax[n];
  int hospitalizations[n];                                            // now included at the level of sewage plants

  int delay_vax;                                                      // delay between vaccination and effectiveness
  int max_delay;                                                      // maximum shift between sewage data and hospitalizations
  int ref_load;                                                       // reference sewage load for hospitalization rates
}

transformed data {
  int hospitalizations_mat[n_date, n_municipality,n_age_group];   
  int population_mat[n_date,n_municipality,n_age_group];
  // matrix [n_date, n_municipality] load_municipality;           // loads by municipality
  real load_mat[n_date,n_municipality];
  real<lower=0, upper = 1> percentage_vax_mat[n_date,n_municipality,n_age_group];
  int daysforanalysis = n_date - max_delay;
  
  for( i in 1:n ){
      hospitalizations_mat[date[i],municipality[i],age_group[i]] = hospitalizations[i];
      population_mat[date[i],municipality[i],age_group[i]] = age_population[i];
      percentage_vax_mat[date[i],municipality[i],age_group[i]] = percentage_vax[i];
      load_mat[date[i], municipality[i]] = load[i];
  }
  
  for( n_obs in 1:n_date){
    // This assumes that the first data is before vaccination started, and we
    // implement the delay already. We start from the last day
    if(n_date - n_obs + 1 - delay_vax <= 0){
      for(i in 1:n_municipality){
        for(j in 1:n_age_group){
          percentage_vax_mat[n_date - n_obs + 1,i,j] = 0;
        }
      }
    } else {
        percentage_vax_mat[n_date - n_obs + 1,,] = percentage_vax_mat[n_date - n_obs + 1 - delay_vax,,];
    }
  }
  
}

parameters {
  /* hospitalization rates and hospitalization delay */
  real<lower = 0> mean_hosp_rate;                                      // hospitalization rates with random effect
  real<lower = 0> sigma_hosp_rate;                                     // hospitalization rates with random effect
  vector<lower = 0>[n_municipality] hosp_rate;                 // hospitalization rates with random effect 
  simplex[n_age_group] hosp_rate_age;
  real<lower = 0, upper = 1> prevention_vax[n_age_group];       // effect of vaccination on hospitalizations
 }

transformed parameters { 
  real log_likes_hospital[daysforanalysis, n_municipality,n_age_group] ;   // log-likelihood contributions of hospitalization data
  real sum_load;
  

  /* notice that this assumes sufficient padding of data to the left */
  for ( i in 1 : n_municipality ) { 
    for(j in 1 : n_age_group){
      for ( n_obs in 1 : daysforanalysis ) { 
	      int ttrue = max_delay + n_obs; 
	      sum_load = 10^(load_mat[ttrue,i] - ref_load);
	      log_likes_hospital[n_obs, i,j] = poisson_lpmf(hospitalizations_mat[ttrue,i,j] | hosp_rate[i] * hosp_rate_age[j] * sum_load * (1 - prevention_vax[j]*percentage_vax_mat[ttrue,i,j]) * population_mat[1,i,j]);
      }
	  } 
  }
}

model {   
  /* (hyper)parameter for hospitalization rates */
  hosp_rate ~ gamma(sigma_hosp_rate * mean_hosp_rate, sigma_hosp_rate); 
 
  /* hyper parameter effect of vaccination */
  prevention_vax ~ uniform(0,1); // beta-distribution would be better
  /* log-likelihood contributions */
  for ( i in 1 : n_municipality ) { 
      for(j in 1 : n_age_group){
        for ( n_obs in 1 : daysforanalysis ) { 
          target += log_likes_hospital[n_obs, i,j]; 
        }
      }
  }

}

generated quantities {
  // real log_lik = sum(log_likes_hospital);  // posterior log-likelihood; mode=1 to calculate WBIC
  real expected_hospitalizations[n_date, n_municipality]; // for calculation of credible intervals
  int simulated_hospitalizations[n_date, n_municipality]; // for calculation of prediction intervals
  real s_load;
  real x;
  
  real expected_hospitalizations_novax[n_date, n_municipality]; // for calculation of credible intervals
  int simulated_hospitalizations_novax[n_date, n_municipality]; // for calculation of prediction intervals


  // posterior quantities
  for ( i in 1 : n_municipality ) {
      for ( t in 1 : n_date ) {	//CHECK - OK
        if ( t > max_delay) {	//simplify for ms/github
		      s_load = 10^(load_mat[t,i] - ref_load);
		      expected_hospitalizations[t, i] = 0;
          simulated_hospitalizations[t, i] = 0;
	        for( j in 1 : n_age_group) {
            x = hosp_rate[i] * hosp_rate_age[j] * s_load * (1 - prevention_vax[j]*percentage_vax_mat[t,i,j]) * population_mat[1,i,j];
            expected_hospitalizations[t, i] = expected_hospitalizations[t, i] + x;
            simulated_hospitalizations[t, i] =simulated_hospitalizations[t, i] + poisson_rng(x);
	        }
        } else {
	      expected_hospitalizations[t, i] = 0;
        simulated_hospitalizations[t, i] = 0;
	  }
   }
  }
  
    for ( i in 1 : n_municipality ) {
      for ( t in 1 : n_date ) {	//CHECK - OK
        if ( t > max_delay) {	//simplify for ms/github
		      s_load = 10^(load_mat[t,i] - ref_load);
		      expected_hospitalizations_novax[t, i] = 0;
          simulated_hospitalizations_novax[t, i] = 0;
	        for( j in 1 : n_age_group) {
            x = hosp_rate[i] * hosp_rate_age[j] * s_load  * population_mat[1,i,j];
            expected_hospitalizations_novax[t, i] = expected_hospitalizations_novax[t, i] + x;
            simulated_hospitalizations_novax[t, i] =simulated_hospitalizations_novax[t, i] + poisson_rng(x);
	        }
        } else {
	      expected_hospitalizations_novax[t, i] = 0;
        simulated_hospitalizations_novax[t, i] = 0;
	  }
   }
  }

}
