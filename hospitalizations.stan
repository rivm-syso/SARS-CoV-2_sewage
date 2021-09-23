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


data {
  int<lower=1> n;
  int<lower = 1> n_municipality;                                      // number of regions (province, safety region, municipality)
  int<lower = 1> n_date;                                              // number of observation days 
  int<lower=1, upper=n_municipality> municipality[n];
  int<lower=1> municipality_pop[n];
  int date[n];
  real<lower=0> load[n];
  int hospitalizations[n];                                            // now included at the level of sewage plants

  int max_delay;                                                      // maximum shift between sewage data and hospitalizations
  int ref_load;                                                       // reference sewage load for hospitalization rates
}

transformed data {
  int hospitalizations_mat[n_date, n_municipality];         
  matrix [n_date, n_municipality] load_municipality;           // loads by municipality
  int municipality_pop_vec[n_municipality];
  real load_mat[n_date,n_municipality];
  int daysforanalysis = n_date - max_delay;
  
  for( i in 1:n_municipality ){
    municipality_pop_vec[i] = 0;
  }
  
  for( i in 1:n ){
      hospitalizations_mat[date[i],municipality[i]] = hospitalizations[i];
      load_mat[date[i], municipality[i]] = load[i];
      municipality_pop_vec[municipality[i]] = municipality_pop_vec[municipality[i]] + municipality_pop[i];
  }
}

parameters {
  /* hospitalization rates and hospitalization delay */
  real<lower = 0> mean_hosp_rate;                                      // hospitalization rates with random effect
  real<lower = 0> sigma_hosp_rate;                                     // hospitalization rates with random effect
  vector<lower = 0>[n_municipality] hosp_rate;                 // hospitalization rates with random effect 
 }

transformed parameters { 
  matrix [daysforanalysis, n_municipality] log_likes_hospital;   // log-likelihood contributions of hospitalization data
  real sum_load;

  /* notice that this assumes sufficient padding of data to the left */
  for ( i in 1 : n_municipality ) { 
    for ( n_obs in 1 : daysforanalysis ) { 
	    int ttrue = max_delay + n_obs; 
	    sum_load = 10^(load_mat[ttrue,i] - ref_load);
	    log_likes_hospital[n_obs, i] = poisson_lpmf(hospitalizations_mat[ttrue,i] | hosp_rate[i] * sum_load * municipality_pop_vec[i]);
	  } 
  }
}

model {   
  /* (hyper)parameter for hospitalization rates */
  hosp_rate ~ gamma(sigma_hosp_rate * mean_hosp_rate, sigma_hosp_rate); 
 
  /* log-likelihood contributions */
  target += sum(log_likes_hospital);  
}

generated quantities {
  real log_lik = sum(log_likes_hospital);  // posterior log-likelihood; mode=1 to calculate WBIC
 /* real expected_hospitalizations[n_date, n_municipality]; // for calculation of credible intervals
  int simulated_hospitalizations[n_date, n_municipality]; // for calculation of prediction intervals
  real s_load;
  real x;

  // posterior quantities
  for ( i in 1 : n_municipality ) {
    for ( t in 1 : n_date ) {	//CHECK - OK
      if ( t > max_delay) {	//simplify for ms/github
		    s_load = 10^(load_municipality[t,i] - ref_load);
        x = hosp_rate[i] * s_load * municipality_pop_vec[i];
	      expected_hospitalizations[t, i] = x;
        simulated_hospitalizations[t, i] = poisson_rng(x);
    } else {
	      expected_hospitalizations[t, i] = 0;
        simulated_hospitalizations[t, i] = 0;		 
	  }
   } 
  }*/ 
} 



