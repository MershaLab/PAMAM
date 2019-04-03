
### install 'rootSolve' package 

install.packages("rootSolve")
library(rootSolve)


####################### Three-way admixture analysis #########################

##Discrete Phenotype: 
##
############## Power and Sampe size: #########################
## Generic functions for case-only study: PowerCase_3way, SampleCase_3way
## Generic functions for case-control study: PowerCaseControl_3way, SampleCaseControl_3way



### Key parameters #####

###  grr  = genotype risk ratio
###  null_prop = c(\theta1, \theta2, \theta3) admixture proportion under the null
###  risk_freq = c(p1, p2, p3), population risk allele frequency of the three ancestry populations 
###  samp_size = sample size for case-only study
###  sample_case, sample_control = cases and control sample sizes, respectively. 
###  alpha = Type I error
###  beta = Type II error
#############################



#### Estimating admixture proportion under alternate hypothesis ####

alter_hypo <- function(null_vector, risk_freq, grr){
  ave_freq <- sum(null_vector*risk_freq)
  fac_nomi <- risk_freq*(grr-1) +1
  fac_deno <- ave_freq*(grr-1) +1
  fact_ratio <- fac_nomi/fac_deno
  alt_value <- null_vector*fact_ratio
  return(alt_value)
}



##################### Case-only study ##########

##################### Power calculation #######

PowerCase_3way <- function(null_vector, risk_freq, grr, alpha, sample_size){
  alt_vector <- alter_hypo(null_vector, risk_freq, grr)
  effect <- sum(((null_vector - alt_vector)^2)/null_vector)
  noncp <- 2*sample_size*effect
  dfreedom <- length(null_vector)-1
  null_critical <- qchisq(p = alpha,df = dfreedom, ncp = 0, lower.tail = F)
  test_power <- pchisq(null_critical, df = dfreedom, ncp = noncp, lower.tail = F)
  return(test_power)
}


################# Sample size estimtion #############

SampleCase_3way <- function(null_vector, risk_freq, grr, alpha, beta){
  dfreedom <- length(null_vector) - 1
  alt_vector <- alter_hypo(null_vector, risk_freq, grr)
  effect <- sum(((null_vector - alt_vector)^2)/null_vector)
  null_critical <- qchisq(p = alpha, df = dfreedom, ncp = 0, lower.tail = F)
  noncp_fun <- function(x) pchisq(null_critical, dfreedom, x, lower.tail = T) - beta
  if(noncp_fun(0) <= 0) samplesize <- 0 else {
    #### null_critical is very small and beta is large so no sample size analyses
    #### Need to adjust the alpha and beta
    uplimit <- 10
    while(noncp_fun(uplimit) >= 0) uplimit <- uplimit*10
    #### Find the upper limit to make noncp_fun < 0 so a unique real root can be found
    solve_noncp <- uniroot(noncp_fun, c(0,uplimit))
    noncp <- solve_noncp$root
    samplesize <- 2*ceiling(0.25*noncp/chi_effect)
  }
  return(samplesize)
}


################# Case control study #####################

####### We use Pearson's Chi-square goodness of fit test to ##########
### estimate power and sample size for three-way admixture mapping.###

################# Power Calculation ######################

case_control_power_3way <- function(null_prop, risk_freq, grr,
                                    sample_case, sample_control, alpha){
alt_prop <- alter_hypo(null_prop, risk_freq, grr)
case_allele_count <- alt_prop* 2*sample_case
control_allele_count <- null_prop*2*sample_control
chi_matrix <- matrix(c(case_allele_count,control_allele_count), nrow = 2, byrow = T)
chi_test <- chisq.test(chi_matrix)
noncp <- chi_test$statistic
dfreedom <- chi_test$parameter
null_critical <- qchisq(p = alpha,df = dfreedom, ncp = 0, lower.tail = F)
chi_power <- pchisq(null_critical, df = dfreedom, ncp = noncp, lower.tail = F)
return(chi_power)
}

################# Sample Size Calculation ######################

case_control_sample_3way <- function(null_prop, risk_freq, grr,
                                     alpha, beta){
  alt_prop <- alter_hypo(null_prop, risk_freq, grr)
  chi_matrix <- matrix(c(alt_prop,null_prop)*10000, nrow = 2, byrow = T)
  #### chi_matrix is based on hypothetical matrix of 10000 cases and 10000 controls
  chi_test <- chisq.test(chi_matrix)
  effect_stat <- chi_test$statistic
  chi_effect <- as.vector((effect_stat/20000))
  #### chi_effect is adjusted back to single sample
  dfreedom <- chi_test$parameter
  null_critical <- qchisq(p = alpha,df = dfreedom, ncp = 0, lower.tail = F)
  noncp_fun <- function(x) pchisq(null_critical,dfreedom, x, lower.tail = T) - beta
  if(noncp_fun(0) <= 0) samplesize <- 0 else {
    #### null_critical is very small and beta is large so no sample size analyses
    uplimit <- 10
    while(noncp_fun(uplimit) >= 0) uplimit <- uplimit*10
    #### Find the upper limit to make noncp_fun < 0 so a unique real root can be found
    solve_noncp <- uniroot(noncp_fun, c(0,uplimit))
    noncp <- solve_noncp$root
    samplesize <- 2*ceiling(0.25*noncp/chi_effect)
  }
  return(samplesize)
}

###################################
####### Quantitative trait ########

#### Key parameters #########
### nullR2 = R-square under null model
### alternateR2 = R-square under alternate model
### u = Numerator degree of freedom = # of ancestral populations - 1
### n_covariates = # of covariates used
### alpha = Type I error rate
### power = 1 - Type II error rate
### N = Sample size


####### Power Calculation #########

quantitative_three_way_power <- function(nullR2, alternateR2, u, alpha, N, n_covariates){
  effectsize <- (alternateR2 - nullR2)/(1-alternateR2)
  v <- N - u - n_covariates - 1
  f_critical <- qf(p = alpha, ncp = 0, df1 = u, df2 = v, lower.tail = F)
  non_central <- effectsize*(u + v + 1)
  powerF <- pf(f_critical, ncp = non_central, df1 = u, df2 = v, lower.tail = F)
  return(round(powerP,4))
}


################# Sample Size Calculation ###############

quantitative_three_way_samplesize <- function(nullR2, alternateR2, u, alpha, power, n_covariates){
  effectsize <- (alternateR2 - nullR2)/(1-alternateR2)
  type2 <- 1-power
  f_critical <- function(x) qf(p = alpha, ncp = 0, df1 = u, df2 = x, lower.tail = F)
  deno_df <- function(x) pf(f_critical(x), ncp = effectsize*(u + x +1), df1 = u, df2 = x, lower.tail = T) - type2
  if(deno_df(1e7) > 0) v <- 1e7 else
    if(deno_df(1e7) < 0) v <- uniroot(deno_df, c(1,1e7))$root
  samplesize <- ceiling(v + u + n_covariates + 1)
 return(samplesize)
}

#################################### Have a nice day #####################3