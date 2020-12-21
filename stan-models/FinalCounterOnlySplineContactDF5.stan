data {
  int <lower=1> M; 
  int <lower=1> N0; 
  int<lower=1> N[M]; 
  int<lower=1> N2; 
  int cases[N2,M]; 
  int deaths[N2, M]; 
  int deathsh[N2, M];
  int deathsc[N2, M];
  int hosp[N2, M]; 
  int hospin[N2, M];
  int hospout[N2, M];
  int icu[N2, M]; 
  int icuin[N2, M];
  int icuout[N2, M];
  matrix[N2,M] contacts;
  matrix[N2, M] f1; // h * s
  matrix[N2, M] f2; // h * s
  matrix[N2, M] f3; // h * s
  matrix[N2, M] f4;
  matrix[N2, M] f5;
  matrix[N2, M] f6;
  matrix[N2, M] f7;
  matrix[N2,M] spline1;
  matrix[N2,M] spline2;
  matrix[N2,M] spline3;
  matrix[N2,M] spline4;
  matrix[N2,M] spline5;
  int EpidemicStart[M];
  real pop[M];
  real SI[N2]; 
}

parameters {
real<lower=0> mu[M]; // intercept for Rt
real<lower=0> y[M];
real beta_1[M];
real beta_2[M];
real beta_3[M];
real beta_4[M];
real beta_5[M];
real beta_6[M];
real<lower=0> phi_dc;
real<lower=0> phi_dh;
real<lower=0> phi_c;
real<lower=0> phi_h;
real<lower=0> phi_hi;
real<lower=0> phi_ho;
real<lower=0> phi_i;
real<lower=0> phi_ii;
real<lower=0> phi_io;
real<lower=0> tau;
real <lower=0> ifr_noise_1[M];
real <lower=0> ifr_noise_2[M];
real <lower=0> ifr_noise_3[M];
real <lower=0> ifr_noise_5[M];
real <lower=0> ifr_noise_7[M];
real<lower=0> var_phi_dc;
real<lower=0> var_phi_dh;
real<lower=0> var_phi_c;
real<lower=0> var_phi_h;
real<lower=0> var_phi_hi;
real<lower=0> var_phi_ho;
real<lower=0> var_phi_i;
real<lower=0> var_phi_ii;
real<lower=0> var_phi_io;
real<lower=0> guma;
}

transformed parameters {
matrix[N2, M] prediction = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso  = rep_matrix(0,N2,M);
matrix[N2, M] Rt = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj = Rt;
matrix[N2, M] lp2 = rep_matrix(0,N2,M);
{
  matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
  
  for (m in 1:M){
   for (i in 2:N0){
  cumm_sum[i,m] = cumm_sum[i-1,m] + y[m];
  }
  prediction[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                  lp2[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+spline4[i,m]*beta_4[m]+spline5[i,m]*beta_5[m]+contacts[i,m]*beta_6[m];
           Rt[i,m] = mu[m]*exp(lp2[i,m]);
        
                   }
              
              
  Rt_adj[1:N0,m] = Rt[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution=0;
  for(j in 1:(i-1)) {
  convolution += prediction[j, m] * SI[i-j];
  }
  cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i-1,m];
  Rt_adj[i,m] = ((pop[m]-cumm_sum[i,m]) / pop[m]) * Rt[i,m];
  prediction[i, m] = Rt_adj[i,m] * convolution;
  }
  
  
  
  E_cases[1, m]= 1e-15 * prediction[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases[i,m] += prediction[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi[1, m]= 1e-15 * prediction[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi[i,m] += prediction[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso[1, m]= 1e-15 * E_hospitalsi[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso[i,m] += E_hospitalsi[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals[1, m] =  E_hospitalsi[1,m]-E_hospitalso[1,m];
  for (i in 2:N2){
  E_hospitals[i,m]+=  E_hospitals[i-1,m]+  E_hospitalsi[i,m]-E_hospitalso[i,m];
  }
  
  E_icusi[1, m]= 1e-15 * E_hospitalsi[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi[i,m] += E_hospitalsi[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso[1, m]= 1e-15 * E_icusi[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso[i,m] += E_icusi[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus[1, m] = E_icusi[1,m]-E_icuso[1,m];
  for (i in 2:N2){
  E_icus[i,m]+=  E_icus[i-1,m]+  E_icusi[i,m]-E_icuso[i,m];
  }
  
  E_deathsc[1, m]=  1e-15 * prediction[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc[i,m] +=  prediction[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh[1, m]=  1e-15*E_hospitalsi[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh[i,m] +=  E_hospitalsi[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths[i,m] += E_deathsh[i,m]+E_deathsc[i,m];
  }
  
  }
}
}
model {
tau ~ exponential(0.03);
for (m in 1:M){
y[m] ~ exponential(1/tau);
}

var_phi_dc ~ normal(8,2);
var_phi_dh ~ normal(8,2);
var_phi_c ~ normal(8,2);//normal(8,2)
var_phi_h ~ normal(8,2);
var_phi_hi ~ normal(8,2);
var_phi_ho ~ normal(8,2);
var_phi_i ~ normal(8,2);
var_phi_ii ~ normal(8,2);
var_phi_io ~ normal(8,2);

phi_dc ~ normal(0,var_phi_dc); // normal(0,5)
phi_dh ~ normal(0,var_phi_dh); // normal(0,5)
phi_c ~ normal(0,var_phi_c);// normal(0,8.5)
phi_h ~ normal(0,var_phi_h);// normal(0,6)
phi_hi ~ normal(0,var_phi_hi);// normal(0,6)
phi_ho ~ normal(0,var_phi_ho);// normal(0,6)
phi_i ~ normal(0,var_phi_i);// normal(0,6)
phi_ii ~ normal(0,var_phi_ii);// normal(0,6)
phi_io ~ normal(0,var_phi_io);// normal(0,6)
for (m in 1:M){
mu[m] ~ normal(3.28, 0.25); 
beta_1[m] ~ normal(0,guma);
beta_2[m] ~ normal(0,guma);
beta_3[m] ~ normal(0,guma);
beta_4[m] ~ normal(0,guma);
beta_5[m] ~ normal(0,guma);
beta_6[m] ~ normal(0,0.5);
}
guma ~ normal(0,0.5);
for (m in 1:M){
ifr_noise_1[m] ~ normal(1,.5);//.1
ifr_noise_2[m] ~ normal(1,.5);
ifr_noise_3[m] ~ normal(1,.5);
ifr_noise_5[m] ~ normal(1,.5);
ifr_noise_7[m] ~ normal(1,.5);
}

for(m in 1:M){
deathsh[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_deathsh[EpidemicStart[m]:N[m], m], phi_dh);
}
for(m in 1:M){
deathsc[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_deathsc[EpidemicStart[m]:N[m], m], phi_dc);
}
for(m in 1:M){
cases[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_cases[EpidemicStart[m]:N[m], m], phi_c);
}

for(m in 1:M){
hosp[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_hospitals[EpidemicStart[m]:N[m], m], phi_h);
}

for(m in 1:M){
hospin[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_hospitalsi[EpidemicStart[m]:N[m], m], phi_hi);
}

for(m in 1:M){
hospout[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_hospitalso[EpidemicStart[m]:N[m], m], phi_ho);
}
for(m in 1:M){
icu[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_icus[EpidemicStart[m]:N[m], m], phi_i);
}

for(m in 1:M){
icuin[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_icusi[EpidemicStart[m]:N[m], m], phi_ii);
}

for(m in 1:M){
icuout[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_icuso[EpidemicStart[m]:N[m], m], phi_io);
}
}

generated quantities {
matrix[N2, M] prediction0 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths0  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh0  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc0  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases0  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals0  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi0  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso0  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus0  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi0  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso0  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_0 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_0 = Rt_0;
matrix[N2, M] lp2_0 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum0 = rep_matrix(0,N2,M);
{
  
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum0[i,m] = cumm_sum0[i-1,m] + y[m];
  }
  prediction0[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
              lp2_0[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+spline4[i,m]*beta_4[m]+spline5[i,m]*beta_5[m]-0.5*beta_6[m];
           Rt_0[i,m] = mu[m]*exp(lp2_0[i,m]);
                }
              
              
  Rt_adj_0[1:N0,m] = Rt_0[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution0=0;
  for(j in 1:(i-1)) {
  convolution0 += prediction0[j, m] * SI[i-j];
  }
  cumm_sum0[i,m] = cumm_sum0[i-1,m] + prediction0[i-1,m];
  Rt_adj_0[i,m] = ((pop[m]-cumm_sum0[i,m]) / pop[m]) * Rt_0[i,m];
  prediction0[i, m] = Rt_adj_0[i,m] * convolution0;
  }
  
  E_cases0[1, m]= 1e-15 * prediction0[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases0[i,m] += prediction0[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi0[1, m]= 1e-15 * prediction0[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi0[i,m] += prediction0[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso0[1, m]= 1e-15 * E_hospitalsi0[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso0[i,m] += E_hospitalsi0[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals0[1, m] =  E_hospitalsi0[1,m]-E_hospitalso0[1,m];
  for (i in 2:N2){
  E_hospitals0[i,m]+=  E_hospitals0[i-1,m]+  E_hospitalsi0[i,m]-E_hospitalso0[i,m];
  }
  
  E_icusi0[1, m]= 1e-15 * E_hospitalsi0[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi0[i,m] += E_hospitalsi0[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso0[1, m]= 1e-15 * E_icusi0[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso0[i,m] += E_icusi0[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus0[1, m] = E_icusi0[1,m]-E_icuso0[1,m];
  for (i in 2:N2){
  E_icus0[i,m]+=  E_icus0[i-1,m]+  E_icusi0[i,m]-E_icuso0[i,m];
  }
  
  E_deathsc0[1, m]=  1e-15 * prediction0[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc0[i,m] +=  prediction0[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh0[1, m]=  1e-15*E_hospitalsi0[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh0[i,m] +=  E_hospitalsi0[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths0[i,m] += E_deathsh0[i,m]+E_deathsc0[i,m];
  }
  
  



  
  }
  

  
}
}
