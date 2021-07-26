/* Stan program for inference using of time- and location-varying        */
/* load of SARS-CoV-2 RNA in Dutch wastewater plants.                    */ 
/* Global load levels can be replaced by global and regional             */
/* loads.                                                                */                       
/* (c) Michiel van Boven, licensed with BSD3 clause (use but mention)    */
/* Created: October 2020                                                 */
/* Updated: January 2021 to integrate with safety region level           */
/* analysis and hospitalization data by admission date and safety region */
/* January 2021: this branch only estimates loads at safety regions      */


/* Changes by Arno:
--- Moved various calculations into the transformed data block
--- Take as much variables as possible from a single data frame,
    to make it easy to work with tidybayes, and to keep the R environment clean
*/

functions {
  /* construct b-spline basis functions; taken from Kharratzadeh's example */
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order); // need the forward declaration  
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t))
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  } 
}

data {
  /* spline data */
  int num_knots;                                                      // number of knots
  int spline_degree;                                                  // spline degree (order - 1)
 
  /* number of installations, regions, and observation times */ 
  int<lower = 1> n;
  int<lower = 1> n_rwzi;                                              // number of wastewater installations 
  int<lower = 1> n_date;                                             // number of observation days 
  int<lower = 1, upper=n_date> date[n];
  real                         concentration[n];
  int<lower = 1, upper=n_rwzi> rwzi[n];
 
  /* sampling temperature */
  int<lower = 0, upper = 1> mode;                                     // 0 = estimation; 1 = WBIC 

}

transformed data {
  
  /* spline parameters */ 
  // Use this when available in v2.27:
  //vector[n_date] knots = linspaced_vector(num_knots, 0, n_date );
  vector[num_knots] knots;

  int num_basis = num_knots + spline_degree - 1;                      // total number of B-splines
  matrix[num_basis, n_date] B;                                        // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2 * spline_degree + num_knots] ext_knots;                    // extended knots
  
  /* sampling temperature */  
  real<lower = 0, upper = 1> watanabe_beta;                           // 0 = normal; 1 = WBIC
  real ts[n_date];
  
  /* wastewater data */
  matrix[n_date, n_rwzi] wastewater;             // -1 = no measurement; 0 = left-censored
  wastewater           = rep_matrix(-1, n_date, n_rwzi );
  for( i in 1:n ){
      wastewater[date[i],rwzi[i]] = concentration[i];
      ts[date[i]] = date[i];
  }
  
  /* construction of splines */	    
  knots[1]=1;
  for( t in 2:num_knots){
    knots[t] = knots[t-1] + (n_date - 1) * 1.0 / (num_knots - 1); //1.0 avoids rouding
  }
  
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1 : num_basis){
    B[ind,:] = to_row_vector(build_b_spline(ts, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  B[num_basis, n_date] = 1;                                     // not needed if ts is sorted 
	
  /* sampling mode: 0 = normal sampling; 1 = WBIC sampling */
  if ( mode == 0 ) {
    watanabe_beta = 1.0;
  }
  else { // mode == 1
    watanabe_beta = 1.0/log(n_rwzi * n_date);    // check carefully; remove -1s
  }
}

parameters {
  /* splines for wastewater intensities */
  matrix [num_basis, n_rwzi] a_individual;  // a_individual[,i] is weights vector for i-th installation 
  vector [num_basis]         a_population;  // a_population is row weights vector
    
  /* parameters of detection function */
  real <lower = 0> x0;                      // x0 of logistic function determining p(positive) 
  real <lower = 2> k;                       // k of logistic function determining p(positive) 
  
  /* sd/var of observations */ 
  real<lower = 0> sigma_observations;       // one sigma/var for all sewage plants

  /* variance of RW1/2 parameter of population load */
  real<lower = 0> RWvar;                    // see Lang & Brezger (2004)
}

transformed parameters { 
  matrix<lower=0> [n_date, n_rwzi] load;    // estimated loads per sewage plant
  matrix [n_date, n_rwzi] log_likes_water;  // log-likelihood contributions of sewage data
  vector [num_basis] weights;               // regression coefficients for the population spline 
      
  
  // P-spline for population level load
  weights[1] = a_population[1];
  for (i in 2 : num_basis) {				        // RW1 smoothing prior (Lang & Brezger, 2004) 
    weights[i] = weights[i-1] + a_population[i];  
  }

  /* estimated loads at the plant level 
      = population load + deviation     */
  load = B' * (rep_matrix( weights, n_rwzi ) + a_individual);
  
  log_likes_water = rep_matrix(0, n_date, n_rwzi );
  
  /* likelihood contributions */
  for ( i in 1 : n_rwzi ) { 
    for ( t in 1 : n_date ) {
	    if (wastewater[t,i] > -0.5) { // if -1 then no measurement
        real alpha = k * (load[t, i] - x0);
        if (wastewater[t, i] < 1.0) // if 0 then no RNA detected
	        log_likes_water[t,i] = bernoulli_logit_lpmf(0 | alpha ); 
	      else 
		      log_likes_water[t,i] = bernoulli_logit_lpmf(1 | alpha ) + normal_lpdf(wastewater[t,i] | load[t, i], sigma_observations);
      }	
    } 
  }
}

model {   
  // prior on RW1 variance
  RWvar ~ inv_gamma(1, 0.0005);                                       // Lang & Brezger (2004) 

  /* population/regional priors */
  a_population[1] ~ normal(12, 1);
  for (s in 2 : num_basis) {
	 a_population[s] ~ normal(0, sqrt(RWvar));                          // tune with LOOIC or p-spline           
  }
  
  /* sewage plant specific priors */
  for (s in 1:num_basis) {
	 a_individual[s] ~ normal(0, 1);
  }

  /* priors for logistic detection function */
  x0 ~ normal(12, 0.5);                                               // a priori 50% detection at log-load 10-14
 
  /* log-likelihood contributions */
  target += watanabe_beta * sum(log_likes_water);  
}

// generated quantities {
//  // real simulated_measurements[n_date, n_rwzi];   // prediction intervals 
//   //real weighted_concentration[n_date];           // posterior weighted population concentration 
//   real log_lik;                                  // posterior log-likelihood; WBIC sampling when mode = 1
//   
//   /* yields estimate of WBIC when mode = 1 */
//   log_lik = sum(log_likes_water);
// 
//  } 
