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
  matrix[N2,M] grocery;
  matrix[N2,M] grocery_1;
  matrix[N2,M] grocery_2;
  matrix[N2,M] grocery_3;
  matrix[N2,M] grocery_4;
  matrix[N2,M] grocery_5;
  matrix[N2,M] grocery_6;
  matrix[N2,M] work;
  matrix[N2,M] work_1;
  matrix[N2,M] work_2;
  matrix[N2,M] work_3;
  matrix[N2,M] work_4;
  matrix[N2,M] work_5;
  matrix[N2,M] work_6;
  matrix[N2,M] home;
  matrix[N2,M] home_1;
  matrix[N2,M] home_2;
  matrix[N2,M] home_3;
  matrix[N2,M] home_4;
  matrix[N2,M] home_5;
  matrix[N2,M] home_6;
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
real alpha_1;
real alpha_2;
real alpha_3;
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
matrix[N2, M] lp = rep_matrix(0,N2,M);
matrix[N2, M] lp2 = rep_matrix(0,N2,M);
{
  matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
  
  for (m in 1:M){
   for (i in 2:N0){
  cumm_sum[i,m] = cumm_sum[i-1,m] + y[m];
  }
  prediction[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp[i,m]=grocery[i,m]*(-(alpha_1))+work[i,m]*(-(alpha_2))+home[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
           Rt[i,m] = mu[m]*2*inv_logit(lp[i,m])*exp(lp2[i,m]);
        
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
alpha_1 ~ normal(0,0.5);
alpha_2 ~ normal(0,0.5);
alpha_3 ~ normal(0,0.5);
for (m in 1:M){
mu[m] ~ normal(3.28, 0.25); 
beta_1[m] ~ normal(0,guma);
beta_2[m] ~ normal(0,guma);
beta_3[m] ~ normal(0,guma);
beta_4[m] ~ normal(0,0.5);
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
matrix[N2, M] prediction1 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths1  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh1  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc1  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases1  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals1  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi1  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso1  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus1  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi1  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso1  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_1 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_1 = Rt_1;
matrix[N2, M] lp_1 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_1 = rep_matrix(0,N2,M);
matrix[N2, M] prediction2 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths2  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh2  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc2  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases2  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals2  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi2  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso2  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus2  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi2  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso2  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_2 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_2 = Rt_2;
matrix[N2, M] lp_2 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_2 = rep_matrix(0,N2,M);
matrix[N2, M] prediction3 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths3  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh3  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc3  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases3  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals3  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi3  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso3  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus3  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi3  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso3  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_3 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_3 = Rt_3;
matrix[N2, M] lp_3 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_3 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum0 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum1 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum2 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum3 = rep_matrix(0,N2,M);
matrix[N2, M] prediction4 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths4  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh4  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc4  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases4  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals4  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi4  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso4  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus4  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi4  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso4  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_4 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_4 = Rt_4;
matrix[N2, M] lp_4 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_4 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum4 = rep_matrix(0,N2,M);
matrix[N2, M] prediction5 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths5  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh5  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc5  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases5  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals5  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi5  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso5  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus5  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi5  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso5  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_5 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_5 = Rt_5;
matrix[N2, M] lp2_5 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum5 = rep_matrix(0,N2,M);
matrix[N2, M] prediction6 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths6  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh6  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc6  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases6  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals6  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi6  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso6  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus6  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi6  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso6  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_6 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_6 = Rt_6;
matrix[N2, M] lp_6 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_6 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum6 = rep_matrix(0,N2,M);
matrix[N2, M] prediction7 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths7  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh7  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc7  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases7  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals7  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi7  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso7  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus7  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi7  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso7  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_7 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_7 = Rt_7;
matrix[N2, M] lp_7 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_7 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum7 = rep_matrix(0,N2,M);

matrix[N2, M] prediction8 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths8  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh8  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc8  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases8  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals8  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi8  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso8  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus8  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi8  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso8  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_8 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_8 = Rt_8;
matrix[N2, M] lp_8 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_8 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum8 = rep_matrix(0,N2,M);

matrix[N2, M] prediction9 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths9  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh9  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc9  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases9  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals9  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi9  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso9  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus9  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi9  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso9  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_9 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_9 = Rt_9;
matrix[N2, M] lp_9 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_9 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum9 = rep_matrix(0,N2,M);

matrix[N2, M] prediction10 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths10  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh10  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc10  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases10  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals10  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi10  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso10  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus10  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi10  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso10  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_10 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_10 = Rt_10;
matrix[N2, M] lp_10 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_10 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum10 = rep_matrix(0,N2,M);


matrix[N2, M] prediction11 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths11  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh11  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc11  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases11  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals11  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi11  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso11  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus11  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi11  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso11  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_11 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_11 = Rt_11;
matrix[N2, M] lp_11 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_11 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum11 = rep_matrix(0,N2,M);


matrix[N2, M] prediction12 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths12  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh12  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc12  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases12  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals12  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi12  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso12  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus12  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi12  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso12  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_12 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_12 = Rt_12;
matrix[N2, M] lp_12 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_12 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum12 = rep_matrix(0,N2,M);


matrix[N2, M] prediction13 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths13  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh13  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc13  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases13  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals13  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi13  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso13  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus13  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi13  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso13  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_13 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_13 = Rt_13;
matrix[N2, M] lp_13 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_13 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum13 = rep_matrix(0,N2,M);


matrix[N2, M] prediction14 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths14  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh14  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc14  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases14  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals14  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi14  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso14  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus14  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi14  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso14  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_14 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_14 = Rt_14;
matrix[N2, M] lp_14 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_14 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum14 = rep_matrix(0,N2,M);


matrix[N2, M] prediction15 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths15  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh15  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc15  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases15  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals15  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi15  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso15  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus15  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi15  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso15  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_15 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_15 = Rt_15;
matrix[N2, M] lp_15 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_15 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum15 = rep_matrix(0,N2,M);


matrix[N2, M] prediction16 = rep_matrix(0,N2,M);
matrix[N2, M] E_deaths16  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsh16  = rep_matrix(0,N2,M);
matrix[N2, M] E_deathsc16  = rep_matrix(0,N2,M);
matrix[N2, M] E_cases16  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitals16  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalsi16  = rep_matrix(0,N2,M);
matrix[N2, M] E_hospitalso16  = rep_matrix(0,N2,M);
matrix[N2, M] E_icus16  = rep_matrix(0,N2,M);
matrix[N2, M] E_icusi16  = rep_matrix(0,N2,M);
matrix[N2, M] E_icuso16  = rep_matrix(0,N2,M);
matrix[N2, M] Rt_16 = rep_matrix(0,N2,M);
matrix[N2, M] Rt_adj_16 = Rt_16;
matrix[N2, M] lp_16 = rep_matrix(0,N2,M);
matrix[N2, M] lp2_16 = rep_matrix(0,N2,M);
matrix[N2, M] cumm_sum16 = rep_matrix(0,N2,M);

{
  
 
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum7[i,m] = cumm_sum7[i-1,m] + y[m];
  }
  prediction7[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_7[i,m]=grocery[i,m]*(-(alpha_1))+work[i,m]*(-(alpha_2))+home[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_7[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_7[i,m] = mu[m]*2*inv_logit(lp_7[i,m])*exp(lp2_7[i,m]);
                 }
              
              
  Rt_adj_7[1:N0,m] = Rt_7[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution7=0;
  for(j in 1:(i-1)) {
  convolution7 += prediction7[j, m] * SI[i-j];
  }
  cumm_sum7[i,m] = cumm_sum7[i-1,m] + prediction7[i-1,m];
  Rt_adj_7[i,m] = ((pop[m]-cumm_sum7[i,m]) / pop[m]) * Rt_7[i,m];
  prediction7[i, m] = Rt_adj_7[i,m] * convolution7;
  }
  
  
  
  E_cases7[1, m]= 1e-15 * prediction7[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases7[i,m] += prediction7[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi7[1, m]= 1e-15 * prediction7[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi7[i,m] += prediction7[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso7[1, m]= 1e-15 * E_hospitalsi7[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso7[i,m] += E_hospitalsi7[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals7[1, m] =  E_hospitalsi7[1,m]-E_hospitalso7[1,m];
  for (i in 2:N2){
  E_hospitals7[i,m]+=  E_hospitals7[i-1,m]+  E_hospitalsi7[i,m]-E_hospitalso7[i,m];
  }
  
  E_icusi7[1, m]= 1e-15 * E_hospitalsi7[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi7[i,m] += E_hospitalsi7[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso7[1, m]= 1e-15 * E_icusi7[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso7[i,m] += E_icusi7[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus7[1, m] = E_icusi7[1,m]-E_icuso7[1,m];
  for (i in 2:N2){
  E_icus7[i,m]+=  E_icus7[i-1,m]+  E_icusi7[i,m]-E_icuso7[i,m];
  }
  
  E_deathsc7[1, m]=  1e-15 * prediction7[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc7[i,m] +=  prediction7[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh7[1, m]=  1e-15*E_hospitalsi7[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh7[i,m] +=  E_hospitalsi7[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths7[i,m] += E_deathsh7[i,m]+E_deathsc7[i,m];
  }
  
  }
  
  
  
  
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum0[i,m] = cumm_sum0[i-1,m] + y[m];
  }
  prediction0[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
              lp2_0[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
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
  
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum5[i,m] = cumm_sum5[i-1,m] + y[m];
  }
  prediction5[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
              lp2_5[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_5[i,m] = mu[m]*exp(lp2_5[i,m]);
                }
              
              
  Rt_adj_5[1:N0,m] = Rt_5[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution5=0;
  for(j in 1:(i-1)) {
  convolution5 += prediction5[j, m] * SI[i-j];
  }
  cumm_sum5[i,m] = cumm_sum5[i-1,m] + prediction5[i-1,m];
  Rt_adj_5[i,m] = ((pop[m]-cumm_sum5[i,m]) / pop[m]) * Rt_5[i,m];
  prediction5[i, m] = Rt_adj_5[i,m] * convolution5;
  }
  
  E_cases5[1, m]= 1e-15 * prediction5[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases5[i,m] += prediction5[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi5[1, m]= 1e-15 * prediction5[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi5[i,m] += prediction5[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso5[1, m]= 1e-15 * E_hospitalsi5[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso5[i,m] += E_hospitalsi5[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals5[1, m] =  E_hospitalsi5[1,m]-E_hospitalso5[1,m];
  for (i in 2:N2){
  E_hospitals5[i,m]+=  E_hospitals5[i-1,m]+  E_hospitalsi5[i,m]-E_hospitalso5[i,m];
  }
  
  E_icusi5[1, m]= 1e-15 * E_hospitalsi5[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi5[i,m] += E_hospitalsi5[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso5[1, m]= 1e-15 * E_icusi5[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso5[i,m] += E_icusi5[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus5[1, m] = E_icusi5[1,m]-E_icuso5[1,m];
  for (i in 2:N2){
  E_icus5[i,m]+=  E_icus5[i-1,m]+  E_icusi5[i,m]-E_icuso5[i,m];
  }
  
  E_deathsc5[1, m]=  1e-15 * prediction5[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc5[i,m] +=  prediction5[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh5[1, m]=  1e-15*E_hospitalsi5[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh5[i,m] +=  E_hospitalsi5[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths5[i,m] += E_deathsh5[i,m]+E_deathsc5[i,m];
  }
  
  



  
  }
  
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum1[i,m] = cumm_sum1[i-1,m] + y[m];
  }
  prediction1[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_1[i,m]=(-1)*(-(alpha_1))+(-1)*(-(alpha_2))+1*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_1[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
            Rt_1[i,m] = mu[m]*2*inv_logit(lp_1[i,m])*exp(lp2_1[i,m]);
                  }
              
              
  Rt_adj_1[1:N0,m] = Rt_1[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution1=0;
  for(j in 1:(i-1)) {
  convolution1 += prediction1[j, m] * SI[i-j];
  }
  cumm_sum1[i,m] = cumm_sum1[i-1,m] + prediction1[i-1,m];
  Rt_adj_1[i,m] = ((pop[m]-cumm_sum1[i,m]) / pop[m]) * Rt_1[i,m];
  prediction1[i, m] = Rt_adj_1[i,m] * convolution1;
  }
  
  
  
  E_cases1[1, m]= 1e-15 * prediction1[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases1[i,m] += prediction1[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi1[1, m]= 1e-15 * prediction1[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi1[i,m] += prediction1[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso1[1, m]= 1e-15 * E_hospitalsi1[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso1[i,m] += E_hospitalsi1[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals1[1, m] =  E_hospitalsi1[1,m]-E_hospitalso1[1,m];
  for (i in 2:N2){
  E_hospitals1[i,m]+=  E_hospitals1[i-1,m]+  E_hospitalsi1[i,m]-E_hospitalso1[i,m];
  }
  
  E_icusi1[1, m]= 1e-15 * E_hospitalsi1[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi1[i,m] += E_hospitalsi1[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso1[1, m]= 1e-15 * E_icusi1[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso1[i,m] += E_icusi1[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus1[1, m] = E_icusi1[1,m]-E_icuso1[1,m];
  for (i in 2:N2){
  E_icus1[i,m]+=  E_icus1[i-1,m]+  E_icusi1[i,m]-E_icuso1[i,m];
  }
  
  E_deathsc1[1, m]=  1e-15 * prediction1[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc1[i,m] +=  prediction1[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh1[1, m]=  1e-15*E_hospitalsi1[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh1[i,m] +=  E_hospitalsi1[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths1[i,m] += E_deathsh1[i,m]+E_deathsc1[i,m];
  }
  
  }
  
  
  
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum6[i,m] = cumm_sum6[i-1,m] + y[m];
  }
  prediction6[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_6[i,m]=(-1)*(-(alpha_1))+(-1)*(-(alpha_2))+1*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_6[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
            Rt_6[i,m] = mu[m]*2*inv_logit(lp_6[i,m])*exp(lp2_6[i,m]);
                  }
              
              
  Rt_adj_6[1:N0,m] = Rt_6[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution6=0;
  for(j in 1:(i-1)) {
  convolution6 += prediction6[j, m] * SI[i-j];
  }
  cumm_sum6[i,m] = cumm_sum6[i-1,m] + prediction6[i-1,m];
  Rt_adj_6[i,m] = ((pop[m]-cumm_sum6[i,m]) / pop[m]) * Rt_6[i,m];
  prediction6[i, m] = Rt_adj_6[i,m] * convolution6;
  }
  
  
  
  E_cases6[1, m]= 1e-15 * prediction6[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases6[i,m] += prediction6[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi6[1, m]= 1e-15 * prediction6[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi6[i,m] += prediction6[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso6[1, m]= 1e-15 * E_hospitalsi6[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso6[i,m] += E_hospitalsi6[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals6[1, m] =  E_hospitalsi6[1,m]-E_hospitalso6[1,m];
  for (i in 2:N2){
  E_hospitals6[i,m]+=  E_hospitals6[i-1,m]+  E_hospitalsi6[i,m]-E_hospitalso6[i,m];
  }
  
  E_icusi6[1, m]= 1e-15 * E_hospitalsi6[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi6[i,m] += E_hospitalsi6[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso6[1, m]= 1e-15 * E_icusi6[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso6[i,m] += E_icusi6[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus6[1, m] = E_icusi6[1,m]-E_icuso6[1,m];
  for (i in 2:N2){
  E_icus6[i,m]+=  E_icus6[i-1,m]+  E_icusi6[i,m]-E_icuso6[i,m];
  }
  
  E_deathsc6[1, m]=  1e-15 * prediction6[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc6[i,m] +=  prediction6[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh6[1, m]=  1e-15*E_hospitalsi6[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh6[i,m] +=  E_hospitalsi6[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths6[i,m] += E_deathsh6[i,m]+E_deathsc6[i,m];
  }
  
  }
  
  
  
  
  
 
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum2[i,m] = cumm_sum2[i-1,m] + y[m];
  }
  prediction2[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_2[i,m]=grocery_1[i,m]*(-(alpha_1))+work_1[i,m]*(-(alpha_2))+home_1[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_2[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
           Rt_2[i,m] = mu[m]*2*inv_logit(lp_2[i,m])*exp(lp2_2[i,m]);
                 }
              
              
  Rt_adj_2[1:N0,m] = Rt_2[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution2=0;
  for(j in 1:(i-1)) {
  convolution2 += prediction2[j, m] * SI[i-j];
  }
  cumm_sum2[i,m] = cumm_sum2[i-1,m] + prediction2[i-1,m];
  Rt_adj_2[i,m] = ((pop[m]-cumm_sum2[i,m]) / pop[m]) * Rt_2[i,m];
  prediction2[i, m] = Rt_adj_2[i,m] * convolution2;
  }
  
  
  
  E_cases2[1, m]= 1e-15 * prediction2[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases2[i,m] += prediction2[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi2[1, m]= 1e-15 * prediction2[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi2[i,m] += prediction2[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso2[1, m]= 1e-15 * E_hospitalsi2[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso2[i,m] += E_hospitalsi2[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals2[1, m] =  E_hospitalsi2[1,m]-E_hospitalso2[1,m];
  for (i in 2:N2){
  E_hospitals2[i,m]+=  E_hospitals2[i-1,m]+  E_hospitalsi2[i,m]-E_hospitalso2[i,m];
  }
  
  E_icusi2[1, m]= 1e-15 * E_hospitalsi2[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi2[i,m] += E_hospitalsi2[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso2[1, m]= 1e-15 * E_icusi2[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso2[i,m] += E_icusi2[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus2[1, m] = E_icusi2[1,m]-E_icuso2[1,m];
  for (i in 2:N2){
  E_icus2[i,m]+=  E_icus2[i-1,m]+  E_icusi2[i,m]-E_icuso2[i,m];
  }
  
  E_deathsc2[1, m]=  1e-15 * prediction2[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc2[i,m] +=  prediction2[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh2[1, m]=  1e-15*E_hospitalsi2[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh2[i,m] +=  E_hospitalsi2[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths2[i,m] += E_deathsh2[i,m]+E_deathsc2[i,m];
  }
  
  }
  
     
     
     
     
     
     
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum8[i,m] = cumm_sum8[i-1,m] + y[m];
  }
  prediction8[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_8[i,m]=grocery_1[i,m]*(-(alpha_1))+work_1[i,m]*(-(alpha_2))+home_1[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_8[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_8[i,m] = mu[m]*2*inv_logit(lp_8[i,m])*exp(lp2_8[i,m]);
                 }
              
              
  Rt_adj_8[1:N0,m] = Rt_8[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution8=0;
  for(j in 1:(i-1)) {
  convolution8 += prediction8[j, m] * SI[i-j];
  }
  cumm_sum8[i,m] = cumm_sum8[i-1,m] + prediction8[i-1,m];
  Rt_adj_8[i,m] = ((pop[m]-cumm_sum8[i,m]) / pop[m]) * Rt_8[i,m];
  prediction8[i, m] = Rt_adj_8[i,m] * convolution8;
  }
  
  
  
  E_cases8[1, m]= 1e-15 * prediction8[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases8[i,m] += prediction8[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi8[1, m]= 1e-15 * prediction8[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi8[i,m] += prediction8[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso8[1, m]= 1e-15 * E_hospitalsi8[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso8[i,m] += E_hospitalsi8[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals8[1, m] =  E_hospitalsi8[1,m]-E_hospitalso8[1,m];
  for (i in 2:N2){
  E_hospitals8[i,m]+=  E_hospitals8[i-1,m]+  E_hospitalsi8[i,m]-E_hospitalso8[i,m];
  }
  
  E_icusi8[1, m]= 1e-15 * E_hospitalsi8[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi8[i,m] += E_hospitalsi8[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso8[1, m]= 1e-15 * E_icusi8[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso8[i,m] += E_icusi8[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus8[1, m] = E_icusi8[1,m]-E_icuso8[1,m];
  for (i in 2:N2){
  E_icus8[i,m]+=  E_icus8[i-1,m]+  E_icusi8[i,m]-E_icuso8[i,m];
  }
  
  E_deathsc8[1, m]=  1e-15 * prediction8[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc8[i,m] +=  prediction8[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh8[1, m]=  1e-15*E_hospitalsi8[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh8[i,m] +=  E_hospitalsi8[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths8[i,m] += E_deathsh8[i,m]+E_deathsc8[i,m];
  }
  
  }
  
  
  
  
  
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum3[i,m] = cumm_sum3[i-1,m] + y[m];
  }
  prediction3[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_3[i,m]=grocery_2[i,m]*(-(alpha_1))+work_2[i,m]*(-(alpha_2))+home_2[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_3[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
           Rt_3[i,m] = mu[m]*2*inv_logit(lp_3[i,m])*exp(lp2_3[i,m]);
                 }
              
              
  Rt_adj_3[1:N0,m] = Rt_3[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution3=0;
  for(j in 1:(i-1)) {
  convolution3 += prediction3[j, m] * SI[i-j];
  }
  cumm_sum3[i,m] = cumm_sum3[i-1,m] + prediction3[i-1,m];
  Rt_adj_3[i,m] = ((pop[m]-cumm_sum3[i,m]) / pop[m]) * Rt_3[i,m];
  prediction3[i, m] = Rt_adj_3[i,m] * convolution3;
  }
  
  
  
  E_cases3[1, m]= 1e-15 * prediction3[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases3[i,m] += prediction3[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi3[1, m]= 1e-15 * prediction3[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi3[i,m] += prediction3[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso3[1, m]= 1e-15 * E_hospitalsi3[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso3[i,m] += E_hospitalsi3[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals3[1, m] =  E_hospitalsi3[1,m]-E_hospitalso3[1,m];
  for (i in 2:N2){
  E_hospitals3[i,m]+=  E_hospitals3[i-1,m]+  E_hospitalsi3[i,m]-E_hospitalso3[i,m];
  }
  
  E_icusi3[1, m]= 1e-15 * E_hospitalsi3[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi3[i,m] += E_hospitalsi3[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso3[1, m]= 1e-15 * E_icusi3[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso3[i,m] += E_icusi3[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus3[1, m] = E_icusi3[1,m]-E_icuso3[1,m];
  for (i in 2:N2){
  E_icus3[i,m]+=  E_icus3[i-1,m]+  E_icusi3[i,m]-E_icuso3[i,m];
  }
  
  E_deathsc3[1, m]=  1e-15 * prediction3[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc3[i,m] +=  prediction3[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh3[1, m]=  1e-15*E_hospitalsi3[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh3[i,m] +=  E_hospitalsi3[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths3[i,m] += E_deathsh3[i,m]+E_deathsc3[i,m];
  }
  
  }
  
  
  
  
  
  for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum9[i,m] = cumm_sum9[i-1,m] + y[m];
  }
  prediction9[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_9[i,m]=grocery_2[i,m]*(-(alpha_1))+work_2[i,m]*(-(alpha_2))+home_2[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_9[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_9[i,m] = mu[m]*2*inv_logit(lp_9[i,m])*exp(lp2_9[i,m]);
                 }
              
              
  Rt_adj_9[1:N0,m] = Rt_9[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution9=0;
  for(j in 1:(i-1)) {
  convolution9 += prediction9[j, m] * SI[i-j];
  }
  cumm_sum9[i,m] = cumm_sum9[i-1,m] + prediction9[i-1,m];
  Rt_adj_9[i,m] = ((pop[m]-cumm_sum9[i,m]) / pop[m]) * Rt_9[i,m];
  prediction9[i, m] = Rt_adj_9[i,m] * convolution9;
  }
  
  
  
  E_cases9[1, m]= 1e-15 * prediction9[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases9[i,m] += prediction9[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi9[1, m]= 1e-15 * prediction9[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi9[i,m] += prediction9[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso9[1, m]= 1e-15 * E_hospitalsi9[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso9[i,m] += E_hospitalsi9[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals9[1, m] =  E_hospitalsi9[1,m]-E_hospitalso9[1,m];
  for (i in 2:N2){
  E_hospitals9[i,m]+=  E_hospitals9[i-1,m]+  E_hospitalsi9[i,m]-E_hospitalso9[i,m];
  }
  
  E_icusi9[1, m]= 1e-15 * E_hospitalsi9[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi9[i,m] += E_hospitalsi9[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso9[1, m]= 1e-15 * E_icusi9[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso9[i,m] += E_icusi9[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus9[1, m] = E_icusi9[1,m]-E_icuso9[1,m];
  for (i in 2:N2){
  E_icus9[i,m]+=  E_icus9[i-1,m]+  E_icusi9[i,m]-E_icuso9[i,m];
  }
  
  E_deathsc9[1, m]=  1e-15 * prediction9[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc9[i,m] +=  prediction9[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh9[1, m]=  1e-15*E_hospitalsi9[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh9[i,m] +=  E_hospitalsi9[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths9[i,m] += E_deathsh9[i,m]+E_deathsc9[i,m];
  }
  
  }
  
  
  
  
  
   for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum4[i,m] = cumm_sum4[i-1,m] + y[m];
  }
  prediction4[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_4[i,m]=grocery_3[i,m]*(-(alpha_1))+work_3[i,m]*(-(alpha_2))+home_3[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_4[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
           Rt_4[i,m] = mu[m]*2*inv_logit(lp_4[i,m])*exp(lp2_4[i,m]);
                 }
              
              
  Rt_adj_4[1:N0,m] = Rt_4[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution4=0;
  for(j in 1:(i-1)) {
  convolution4 += prediction4[j, m] * SI[i-j];
  }
  cumm_sum4[i,m] = cumm_sum4[i-1,m] + prediction4[i-1,m];
  Rt_adj_4[i,m] = ((pop[m]-cumm_sum4[i,m]) / pop[m]) * Rt_4[i,m];
  prediction4[i, m] = Rt_adj_4[i,m] * convolution4;
  }
  
  
  
  E_cases4[1, m]= 1e-15 * prediction4[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases4[i,m] += prediction4[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi4[1, m]= 1e-15 * prediction4[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi4[i,m] += prediction4[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso4[1, m]= 1e-15 * E_hospitalsi4[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso4[i,m] += E_hospitalsi4[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals4[1, m] =  E_hospitalsi4[1,m]-E_hospitalso4[1,m];
  for (i in 2:N2){
  E_hospitals4[i,m]+=  E_hospitals4[i-1,m]+  E_hospitalsi4[i,m]-E_hospitalso4[i,m];
  }
  
  E_icusi4[1, m]= 1e-15 * E_hospitalsi4[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi4[i,m] += E_hospitalsi4[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso4[1, m]= 1e-15 * E_icusi4[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso4[i,m] += E_icusi4[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus4[1, m] = E_icusi4[1,m]-E_icuso4[1,m];
  for (i in 2:N2){
  E_icus4[i,m]+=  E_icus4[i-1,m]+  E_icusi4[i,m]-E_icuso4[i,m];
  }
  
  E_deathsc4[1, m]=  1e-15 * prediction4[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc4[i,m] +=  prediction4[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh4[1, m]=  1e-15*E_hospitalsi4[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh4[i,m] +=  E_hospitalsi4[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths4[i,m] += E_deathsh4[i,m]+E_deathsc4[i,m];
  }
  
  }
  
  
  
  
  
  
   for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum10[i,m] = cumm_sum10[i-1,m] + y[m];
  }
  prediction10[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_10[i,m]=grocery_3[i,m]*(-(alpha_1))+work_3[i,m]*(-(alpha_2))+home_3[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_10[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_10[i,m] = mu[m]*2*inv_logit(lp_10[i,m])*exp(lp2_10[i,m]);
                 }
              
              
  Rt_adj_10[1:N0,m] = Rt_10[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution10=0;
  for(j in 1:(i-1)) {
  convolution10 += prediction10[j, m] * SI[i-j];
  }
  cumm_sum10[i,m] = cumm_sum10[i-1,m] + prediction10[i-1,m];
  Rt_adj_10[i,m] = ((pop[m]-cumm_sum10[i,m]) / pop[m]) * Rt_10[i,m];
  prediction10[i, m] = Rt_adj_10[i,m] * convolution10;
  }
  
  
  
  E_cases10[1, m]= 1e-15 * prediction10[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases10[i,m] += prediction10[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi10[1, m]= 1e-15 * prediction10[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi10[i,m] += prediction10[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso10[1, m]= 1e-15 * E_hospitalsi10[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso10[i,m] += E_hospitalsi10[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals10[1, m] =  E_hospitalsi10[1,m]-E_hospitalso10[1,m];
  for (i in 2:N2){
  E_hospitals10[i,m]+=  E_hospitals10[i-1,m]+  E_hospitalsi10[i,m]-E_hospitalso10[i,m];
  }
  
  E_icusi10[1, m]= 1e-15 * E_hospitalsi10[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi10[i,m] += E_hospitalsi10[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso10[1, m]= 1e-15 * E_icusi10[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso10[i,m] += E_icusi10[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus10[1, m] = E_icusi10[1,m]-E_icuso10[1,m];
  for (i in 2:N2){
  E_icus10[i,m]+=  E_icus10[i-1,m]+  E_icusi10[i,m]-E_icuso10[i,m];
  }
  
  E_deathsc10[1, m]=  1e-15 * prediction10[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc10[i,m] +=  prediction10[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh10[1, m]=  1e-15*E_hospitalsi10[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh10[i,m] +=  E_hospitalsi10[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths10[i,m] += E_deathsh10[i,m]+E_deathsc10[i,m];
  }
  
  }




  
   for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum11[i,m] = cumm_sum11[i-1,m] + y[m];
  }
  prediction11[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_11[i,m]=grocery_4[i,m]*(-(alpha_1))+work_4[i,m]*(-(alpha_2))+home_4[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_11[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
           Rt_11[i,m] = mu[m]*2*inv_logit(lp_11[i,m])*exp(lp2_11[i,m]);
                 }
              
              
  Rt_adj_11[1:N0,m] = Rt_11[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution11=0;
  for(j in 1:(i-1)) {
  convolution11 += prediction11[j, m] * SI[i-j];
  }
  cumm_sum11[i,m] = cumm_sum11[i-1,m] + prediction11[i-1,m];
  Rt_adj_11[i,m] = ((pop[m]-cumm_sum11[i,m]) / pop[m]) * Rt_11[i,m];
  prediction11[i, m] = Rt_adj_11[i,m] * convolution11;
  }
  
  
  
  E_cases11[1, m]= 1e-15 * prediction11[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases11[i,m] += prediction11[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi11[1, m]= 1e-15 * prediction11[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi11[i,m] += prediction11[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso11[1, m]= 1e-15 * E_hospitalsi11[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso11[i,m] += E_hospitalsi11[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals11[1, m] =  E_hospitalsi11[1,m]-E_hospitalso11[1,m];
  for (i in 2:N2){
  E_hospitals11[i,m]+=  E_hospitals11[i-1,m]+  E_hospitalsi11[i,m]-E_hospitalso11[i,m];
  }
  
  E_icusi11[1, m]= 1e-15 * E_hospitalsi11[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi11[i,m] += E_hospitalsi11[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso11[1, m]= 1e-15 * E_icusi11[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso11[i,m] += E_icusi11[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus11[1, m] = E_icusi11[1,m]-E_icuso11[1,m];
  for (i in 2:N2){
  E_icus11[i,m]+=  E_icus11[i-1,m]+  E_icusi11[i,m]-E_icuso11[i,m];
  }
  
  E_deathsc11[1, m]=  1e-15 * prediction11[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc11[i,m] +=  prediction11[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh11[1, m]=  1e-15*E_hospitalsi11[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh11[i,m] +=  E_hospitalsi11[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths11[i,m] += E_deathsh11[i,m]+E_deathsc11[i,m];
  }
  
  }





  
   for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum12[i,m] = cumm_sum12[i-1,m] + y[m];
  }
  prediction12[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_12[i,m]=grocery_4[i,m]*(-(alpha_1))+work_4[i,m]*(-(alpha_2))+home_4[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_12[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_12[i,m] = mu[m]*2*inv_logit(lp_12[i,m])*exp(lp2_12[i,m]);
                 }
              
              
  Rt_adj_12[1:N0,m] = Rt_12[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution12=0;
  for(j in 1:(i-1)) {
  convolution12 += prediction12[j, m] * SI[i-j];
  }
  cumm_sum12[i,m] = cumm_sum12[i-1,m] + prediction12[i-1,m];
  Rt_adj_12[i,m] = ((pop[m]-cumm_sum12[i,m]) / pop[m]) * Rt_12[i,m];
  prediction12[i, m] = Rt_adj_12[i,m] * convolution12;
  }
  
  
  
  E_cases12[1, m]= 1e-15 * prediction12[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases12[i,m] += prediction12[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi12[1, m]= 1e-15 * prediction12[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi12[i,m] += prediction12[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso12[1, m]= 1e-15 * E_hospitalsi12[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso12[i,m] += E_hospitalsi12[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals12[1, m] =  E_hospitalsi12[1,m]-E_hospitalso12[1,m];
  for (i in 2:N2){
  E_hospitals12[i,m]+=  E_hospitals12[i-1,m]+  E_hospitalsi12[i,m]-E_hospitalso12[i,m];
  }
  
  E_icusi12[1, m]= 1e-15 * E_hospitalsi12[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi12[i,m] += E_hospitalsi12[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso12[1, m]= 1e-15 * E_icusi12[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso12[i,m] += E_icusi12[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus12[1, m] = E_icusi12[1,m]-E_icuso12[1,m];
  for (i in 2:N2){
  E_icus12[i,m]+=  E_icus12[i-1,m]+  E_icusi12[i,m]-E_icuso12[i,m];
  }
  
  E_deathsc12[1, m]=  1e-15 * prediction12[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc12[i,m] +=  prediction12[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh12[1, m]=  1e-15*E_hospitalsi12[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh12[i,m] +=  E_hospitalsi12[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths12[i,m] += E_deathsh12[i,m]+E_deathsc12[i,m];
  }
  
  }
  




   for (m in 1:M){
 
  for (i in 2:N0){
  cumm_sum13[i,m] = cumm_sum13[i-1,m] + y[m];
  }
  prediction13[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
  
               for (i in 1:N2){
                 lp_13[i,m]=grocery_5[i,m]*(-(alpha_1))+work_5[i,m]*(-(alpha_2))+home_5[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_13[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
           Rt_13[i,m] = mu[m]*2*inv_logit(lp_13[i,m])*exp(lp2_13[i,m]);
                 }
              
              
  Rt_adj_13[1:N0,m] = Rt_13[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution13=0;
  for(j in 1:(i-1)) {
  convolution13 += prediction13[j, m] * SI[i-j];
  }
  cumm_sum13[i,m] = cumm_sum13[i-1,m] + prediction13[i-1,m];
  Rt_adj_13[i,m] = ((pop[m]-cumm_sum13[i,m]) / pop[m]) * Rt_13[i,m];
  prediction13[i, m] = Rt_adj_13[i,m] * convolution13;
  }
  
  
  
  E_cases13[1, m]= 1e-15 * prediction13[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases13[i,m] += prediction13[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }
  
  
  E_hospitalsi13[1, m]= 1e-15 * prediction13[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi13[i,m] += prediction13[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }
  
  E_hospitalso13[1, m]= 1e-15 * E_hospitalsi13[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso13[i,m] += E_hospitalsi13[j,m] * f4[i-j,m] ;
  }
  }
  
  
  E_hospitals13[1, m] =  E_hospitalsi13[1,m]-E_hospitalso13[1,m];
  for (i in 2:N2){
  E_hospitals13[i,m]+=  E_hospitals13[i-1,m]+  E_hospitalsi13[i,m]-E_hospitalso13[i,m];
  }
  
  E_icusi13[1, m]= 1e-15 * E_hospitalsi13[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi13[i,m] += E_hospitalsi13[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }
  
  E_icuso13[1, m]= 1e-15 * E_icusi13[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso13[i,m] += E_icusi13[j,m] * f6[i-j,m] ;
  }
  }
  
  
  E_icus13[1, m] = E_icusi13[1,m]-E_icuso13[1,m];
  for (i in 2:N2){
  E_icus13[i,m]+=  E_icus13[i-1,m]+  E_icusi13[i,m]-E_icuso13[i,m];
  }
  
  E_deathsc13[1, m]=  1e-15 * prediction13[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc13[i,m] +=  prediction13[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }
  
  E_deathsh13[1, m]=  1e-15*E_hospitalsi13[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh13[i,m] +=  E_hospitalsi13[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }
  
  
  for (i in 1:N2){
  E_deaths13[i,m] += E_deathsh13[i,m]+E_deathsc13[i,m];
  }
  
  }




   for (m in 1:M){

  for (i in 2:N0){
  cumm_sum14[i,m] = cumm_sum14[i-1,m] + y[m];
  }
  prediction14[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days

               for (i in 1:N2){
                 lp_14[i,m]=grocery_5[i,m]*(-(alpha_1))+work_5[i,m]*(-(alpha_2))+home_5[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_14[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_14[i,m] = mu[m]*2*inv_logit(lp_14[i,m])*exp(lp2_14[i,m]);
                 }


  Rt_adj_14[1:N0,m] = Rt_14[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution14=0;
  for(j in 1:(i-1)) {
  convolution14 += prediction14[j, m] * SI[i-j];
  }
  cumm_sum14[i,m] = cumm_sum14[i-1,m] + prediction14[i-1,m];
  Rt_adj_14[i,m] = ((pop[m]-cumm_sum14[i,m]) / pop[m]) * Rt_14[i,m];
  prediction14[i, m] = Rt_adj_14[i,m] * convolution14;
  }



  E_cases14[1, m]= 1e-15 * prediction14[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases14[i,m] += prediction14[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }


  E_hospitalsi14[1, m]= 1e-15 * prediction14[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi14[i,m] += prediction14[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }

  E_hospitalso14[1, m]= 1e-15 * E_hospitalsi14[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso14[i,m] += E_hospitalsi14[j,m] * f4[i-j,m] ;
  }
  }


  E_hospitals14[1, m] =  E_hospitalsi14[1,m]-E_hospitalso14[1,m];
  for (i in 2:N2){
  E_hospitals14[i,m]+=  E_hospitals14[i-1,m]+  E_hospitalsi14[i,m]-E_hospitalso14[i,m];
  }

  E_icusi14[1, m]= 1e-15 * E_hospitalsi14[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi14[i,m] += E_hospitalsi14[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }

  E_icuso14[1, m]= 1e-15 * E_icusi14[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso14[i,m] += E_icusi14[j,m] * f6[i-j,m] ;
  }
  }


  E_icus14[1, m] = E_icusi14[1,m]-E_icuso14[1,m];
  for (i in 2:N2){
  E_icus14[i,m]+=  E_icus14[i-1,m]+  E_icusi14[i,m]-E_icuso14[i,m];
  }

  E_deathsc14[1, m]=  1e-15 * prediction14[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc14[i,m] +=  prediction14[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }

  E_deathsh14[1, m]=  1e-15*E_hospitalsi14[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh14[i,m] +=  E_hospitalsi14[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }


  for (i in 1:N2){
  E_deaths14[i,m] += E_deathsh14[i,m]+E_deathsc14[i,m];
  }

  }




   for (m in 1:M){

  for (i in 2:N0){
  cumm_sum15[i,m] = cumm_sum15[i-1,m] + y[m];
  }
  prediction15[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days

               for (i in 1:N2){
                 lp_15[i,m]=grocery_6[i,m]*(-(alpha_1))+work_6[i,m]*(-(alpha_2))+home_6[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_15[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]+contacts[i,m]*beta_4[m];
           Rt_15[i,m] = mu[m]*2*inv_logit(lp_15[i,m])*exp(lp2_15[i,m]);
                 }


  Rt_adj_15[1:N0,m] = Rt_15[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution15=0;
  for(j in 1:(i-1)) {
  convolution15 += prediction15[j, m] * SI[i-j];
  }
  cumm_sum15[i,m] = cumm_sum15[i-1,m] + prediction15[i-1,m];
  Rt_adj_15[i,m] = ((pop[m]-cumm_sum15[i,m]) / pop[m]) * Rt_15[i,m];
  prediction15[i, m] = Rt_adj_15[i,m] * convolution15;
  }



  E_cases15[1, m]= 1e-15 * prediction15[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases15[i,m] += prediction15[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }


  E_hospitalsi15[1, m]= 1e-15 * prediction15[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi15[i,m] += prediction15[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }

  E_hospitalso15[1, m]= 1e-15 * E_hospitalsi15[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso15[i,m] += E_hospitalsi15[j,m] * f4[i-j,m] ;
  }
  }


  E_hospitals15[1, m] =  E_hospitalsi15[1,m]-E_hospitalso15[1,m];
  for (i in 2:N2){
  E_hospitals15[i,m]+=  E_hospitals15[i-1,m]+  E_hospitalsi15[i,m]-E_hospitalso15[i,m];
  }

  E_icusi15[1, m]= 1e-15 * E_hospitalsi15[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi15[i,m] += E_hospitalsi15[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }

  E_icuso15[1, m]= 1e-15 * E_icusi15[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso15[i,m] += E_icusi15[j,m] * f6[i-j,m] ;
  }
  }


  E_icus15[1, m] = E_icusi15[1,m]-E_icuso15[1,m];
  for (i in 2:N2){
  E_icus15[i,m]+=  E_icus15[i-1,m]+  E_icusi15[i,m]-E_icuso15[i,m];
  }

  E_deathsc15[1, m]=  1e-15 * prediction15[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc15[i,m] +=  prediction15[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }

  E_deathsh15[1, m]=  1e-15*E_hospitalsi15[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh15[i,m] +=  E_hospitalsi15[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }


  for (i in 1:N2){
  E_deaths15[i,m] += E_deathsh15[i,m]+E_deathsc15[i,m];
  }

  }



   for (m in 1:M){

  for (i in 2:N0){
  cumm_sum16[i,m] = cumm_sum16[i-1,m] + y[m];
  }
  prediction16[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days

               for (i in 1:N2){
                 lp_16[i,m]=grocery_6[i,m]*(-(alpha_1))+work_6[i,m]*(-(alpha_2))+home_6[i,m]*(-(alpha_3));//+schools[i,m]*(-(alpha_4));
                 lp2_16[i,m]=spline1[i,m]*beta_1[m]+spline2[i,m]*beta_2[m]+spline3[i,m]*beta_3[m]-0.5*beta_4[m];
           Rt_16[i,m] = mu[m]*2*inv_logit(lp_16[i,m])*exp(lp2_16[i,m]);
                 }


  Rt_adj_16[1:N0,m] = Rt_16[1:N0,m];
  for (i in (N0+1):N2) {
  real convolution16=0;
  for(j in 1:(i-1)) {
  convolution16 += prediction16[j, m] * SI[i-j];
  }
  cumm_sum16[i,m] = cumm_sum16[i-1,m] + prediction16[i-1,m];
  Rt_adj_16[i,m] = ((pop[m]-cumm_sum16[i,m]) / pop[m]) * Rt_16[i,m];
  prediction16[i, m] = Rt_adj_16[i,m] * convolution16;
  }



  E_cases16[1, m]= 1e-15 * prediction16[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_cases16[i,m] += prediction16[j,m] * f1[i-j,m] * ifr_noise_1[m];
  }
  }


  E_hospitalsi16[1, m]= 1e-15 * prediction16[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalsi16[i,m] += prediction16[j,m] * f3[i-j,m] * ifr_noise_3[m];
  }
  }

  E_hospitalso16[1, m]= 1e-15 * E_hospitalsi16[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_hospitalso16[i,m] += E_hospitalsi16[j,m] * f4[i-j,m] ;
  }
  }


  E_hospitals16[1, m] =  E_hospitalsi16[1,m]-E_hospitalso16[1,m];
  for (i in 2:N2){
  E_hospitals16[i,m]+=  E_hospitals16[i-1,m]+  E_hospitalsi16[i,m]-E_hospitalso16[i,m];
  }

  E_icusi16[1, m]= 1e-15 * E_hospitalsi16[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icusi16[i,m] += E_hospitalsi16[j,m] * f5[i-j,m] * ifr_noise_5[m];
  }
  }

  E_icuso16[1, m]= 1e-15 * E_icusi16[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_icuso16[i,m] += E_icusi16[j,m] * f6[i-j,m] ;
  }
  }


  E_icus16[1, m] = E_icusi16[1,m]-E_icuso16[1,m];
  for (i in 2:N2){
  E_icus16[i,m]+=  E_icus16[i-1,m]+  E_icusi16[i,m]-E_icuso16[i,m];
  }

  E_deathsc16[1, m]=  1e-15 * prediction16[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsc16[i,m] +=  prediction16[j,m] * f2[i-j,m] * ifr_noise_2[m];
  }
  }

  E_deathsh16[1, m]=  1e-15*E_hospitalsi16[1,m];
  for (i in 2:N2){
  for(j in 1:(i-1)){
  E_deathsh16[i,m] +=  E_hospitalsi16[j,m] * f7[i-j,m] * ifr_noise_7[m];
  }
  }


  for (i in 1:N2){
  E_deaths16[i,m] += E_deathsh16[i,m]+E_deathsc16[i,m];
  }

  }


}
}
