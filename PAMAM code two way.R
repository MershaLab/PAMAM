
## PAMAM power analysis R code for two-way admixture mapping ##
## The codes are based on the following manuscript - 
### Gautam et al. (2017), AdmixPower: Statistical power and sample size ###
######  estimation for mapping genetic loci in admixed populations.########
################## Genetics,vol. 207 no. 3 873-882. ####################### 

####################### Two-way admixture analysis #########################
##
## Discrete Phenotype: 
##
############## Power and Sampe size: #########################
## Generic functions for case-only study: PowerCase, SampleCase
## Generic functions for case-control study: PowerCaseControl, SampleCaseControl
#
## Inputs: null.prop, alt.prop, case.n, control.n, type1.error, type2.error, side
## null.prop = Admixture proportion under null hypothesis
## alt.prop = Admixture proportion under alternate hypothesis
## case.n = # of cases; control.n = # of controls
## type1.error = type I error rate adjusted for multiple testing
## type2.error = type II error rate = 1-power
## side = 1 if one-sided  or  = 2 if two-sided hypothesis

####################################
##          Function: PowerCase
####################################
PowerCase <- function(null.prop, alt.prop, case.n, type1.error, side){
  V0 <- null.prop*(1 - null.prop)/(2*case.n)
  V1 <- alt.prop*(1 - alt.prop)/(null.prop - null.prop^2)
  if (side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE) 
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  nume <- sqrt(V0)*z.type1.error - abs(alt.prop - null.prop)
  deno <- sqrt(V0*V1)
  power.q <- nume/deno
  power <- pnorm(power.q, 0, 1, lower.tail = FALSE, log.p = FALSE)
  return(power)
}

####################################
##          Function: SampleCase
####################################
SampleCase <- function(null.prop, alt.prop, type1.error, type2.error, side){
  z.type2.error <- qnorm(type2.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  if (side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE) 
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  c1 <- sqrt(alt.prop - alt.prop^2); c2 <- sqrt(null.prop - null.prop^2)
  nume <- z.type2.error*c1 + z.type1.error*c2
  deno <- alt.prop - null.prop
  size <- ceiling(0.5*(nume/deno)^2)
  return(size)
}

#########################################
##          Function: PowerCaseControl
#########################################
PowerCaseControl <- function(null.prop, alt.prop, case.n, control.n, type1.error, side){
  V01 <- null.prop*(1 - null.prop)/(2*case.n); V02 <- null.prop*(1 - null.prop)/(2*control.n)
  V11 <- alt.prop*(1 - alt.prop)/(2*case.n); V12 <- null.prop*(1 - null.prop)/(2*control.n)
  V0 <- V01 + V02; V1 <- (V11 + V12)/V0
  if (side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  nume <- sqrt(V0)*z.type1.error - abs(alt.prop - null.prop)
  deno <- sqrt(V0*V1)
  power.q <- nume/deno
  power <- pnorm(power.q, 0, 1, lower.tail = FALSE, log.p = FALSE)
  return(power)
}

#########################################
##          Function: SampleCaseControl
#########################################

SampleCaseControl <- function(null.prop, alt.prop, type1.error, type2.error, side){
  z.type2.error <- qnorm(type2.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  if (side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  c1 <- sqrt((alt.prop - alt.prop^2) + (null.prop - null.prop^2))
  c2 <- sqrt(2*(null.prop - null.prop^2))
  nume <- z.type2.error*c1 + z.type1.error*c2
  deno <- alt.prop - null.prop
  size <- 2*ceiling(0.5*(nume/deno)^2)
  return(size)
}

################ Estimating null.prop and alt.prop ##################
####    The estimate of null.prop and alt.prop depend of    ####
####       different population specific parameters,        ####
####          disease model, and admixture model.           ####
##################################################################### 

############# Methpd by Montana and Pritchard(2004): ################
## Multiplicative mode, Hybrid-isolation
## Population specific parameters: admix.prop, grr, px, py
## grr  = genotypr risk ratio
## px = risk allele frequency in population X
## py = risk allele frequency in population Y
## Functions:
## PowerDiscreteGRR
## SampleDiscreteGRR
####################

AdmixPropGRR <- function(admix.prop, grr, px, py){
  p.bar <- px*admix.prop + py*(1 - admix.prop)
  alt.prop <- admix.prop*(1 + px*grr - px)/(1 + p.bar*grr - p.bar)
  return(c(admix.prop,alt.prop))
}

PowerDiscreteGRR <- function(admix.prop, grr, px, py, case.n, control.n, type1.error, side, study.design){
  prop <- AdmixPropGRR(admix.prop, grr, px, py)
  if (study.design == 1) {
    power <- PowerCase(prop[1], prop[2], case.n, type1.error, side)
  } else {
    power <- PowerCaseControl(prop[1], prop[2], case.n, control.n, type1.error, side)
  }
  return(power)
}

SampleDiscreteGRR <- function(admix.prop, grr, px, py, type1.error, type2.error, side, study.design){
  prop <- AdmixPropGRR(admix.prop, grr, px, py)
  if (study.design == 1 ) {
    sample.size <- SampleCase(prop[1], prop[2], type1.error, type2.error, side)
  } else {
    sample.size <- SampleCaseControl(prop[1], prop[2], type1.error, type2.error, side)
  }
  return(sample.size)
}


################ Method by Zhu et al. (2004): ################## 

## Admixture Process - Hybrid-Isolation and Continuous Gene Flow
## Multiplicative, Additive, Reccessive, and Dominance model of inheritance 
## Number of generation since admixture = g, recombination rate = recom.rate
## generation = # of generation; recom.rate = recombination rate; 
## admix.prop = admixture proportion; prr = parental risk ratio
## Extented for Case-Control study.design
## Functions:
## PowerDiscretePRR
## SampleDiscretePRR
##################################################################### 
## Given the penetrances and ancestral allele frequency, compute the parental risk ratio
## Penetrance function = (f0, f1, f2) ## a vector object of size 3

PrrPenetrance <- function(penetrance, px, py){
  f0 <- penetrance[1]; f1 <- penetrance[2]; f2 <- penetrance[3]
  fac1 <- c(f2*px^2, f1*2*px(1 - px), f0*(1 - px)^2)
  fac2 <- c(f2*py^2, f1*2*py(1 - py), f0*(1 - py)^2)
  r <- sum(fac1)/sum(fac2)
  return(r)
}

## Given the disease prevalences in the parental population, compute the parental risk ratio

PrrPrevalence <- function(prevalence.x, prevalence.y){
  rr <- prevalence.x/prevalence.y
  return(rr)
}

###############################################################
##
## Hybrid-Isolation (HI):
## Contribution from the population X at the beginning of admixture process = admix.prop  
##
##############################################################
##
## HI1. HI Multiplicative mode:
##
##########################################

AdmixPropMultiHI <- function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2 - 2*admix.prop
  a1 <- (1 - recom.rate)^(g - 2)
  fac1 <- 0.5*(2 - gamma - 2*recom.rate)*gamma*a1
  fac2 <- (sqrt(prr) - 1)/((2 - gamma)*sqrt(prr) + gamma)
  est.prop <- admix.prop + fac1*fac2
  return(est.prop)
}

###########################################
##
## HI2. HI additive mode:
##
###########################################

AdmixPropAddHI <- function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2 - 2*admix.prop
  fac1 <- 0.25*(2 - gamma - 2*recom.rate)*gamma*(1 - recom.rate)^(g - 2)
  fac2 <- (2 - gamma)*prr + gamma
  est.prop <- admix.prop + fac1*(prr - 1)/fac2
  return(est.prop)
}

###############################################
##
## HI3. HI recessive mode:Penetrance f0 = f1 and py = 0
##
###############################################

AdmixPropRecHI<-function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2 - 2*admix.prop
  fac1 <- 0.5*(2 - gamma - 2*recom.rate)*gamma*(1 - recom.rate)^(g - 2)
  fac2 <- ((2 - gamma)^2)*(prr - 1) + 4
  est.prop <- admix.prop + (2 - gamma)*(prr - 1)*fac1/fac2
  return(est.prop)
}

#################################################
##
## HI4. HI Dominant mode:Penetrance f1 = f2, px = 1
##
#################################################

AdmixPropDomHI <- function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2 - 2*admix.prop
  fac1 <- 0.5*(2 - gamma - 2*recom.rate)*gamma*(1 - recom.rate)^(g - 2)
  fac2 <- (1 - prr)*gamma^2 + 4*prr
  est.prop <- admix.prop + gamma*(prr - 1)*fac1/fac2
  return(est.prop)
}

##############################################################################
##
## CGF admixture process:Contribution from Y is continuous in each generation
## Contribution from the Y in the admixted population per generation: gamma/2
## After g generation, the contribution from population X  = admix.prop = (1- gamma/2)^g
## Contribution from population Y  = 1 - (admix.prop)
## gamma = 2*(1 - admix.prop^(1/g))
##
###############################################################################
##
## CGF1. CGF Multiplicative Mode:
##
#######################################

AdmixPropMultiCGF <- function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2*(1 - admix.prop^(1/g))
  F1 <- (1 - gamma/2)^g
  rr <- sqrt(prr) - 1
  F2 <- (1 - gamma)*F1*rr
  F3.1 <- recom.rate*(1 - gamma)*F1; F3.2 <- 0.5*gamma*(1 - recom.rate - 0.5*gamma)*(1 - recom.rate)^g
  F3.3 <- (F3.1 - F3.2)/(recom.rate - gamma/2)
  F3 <- F3.3*rr/(1 - 0.5*gamma)
  F4 <- (1 - gamma)*F1*rr
  F5 <- F1*rr
  num <- F1*(F2 + 1 - 0.5*gamma)*(1 + F3)
  deno <- (F4 + 1)*(F5 + 1)
  allele.prop <- num/deno
  return(allele.prop)
}

#####################################
##
## CGF2. CGF Additive Mode:
##
#####################################

AdmixPropAddCGF <- function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2*(1 - admix.prop^(1/g))		
  F1 <- (1 - gamma/2)^g
  F2 <- (1 - gamma)*F1*(prr - 1)
  F3.1 <- recom.rate*(1 - gamma)*F1; F3.2 <- 0.5*gamma*(1 - recom.rate - 0.5*gamma)*(1 - recom.rate)^g
  F3.3 <- (F3.1 - F3.2)/(recom.rate - gamma/2)
  F3 <- F3.3*(prr - 1)
  F4 <- (2 - gamma)*F1*(prr - 1)
  num <- F1*(0.5*F2 + 0.5*F3 + 1 - 0.5*gamma)
  deno <- 0.5* F4 + 1
  allele.prop <- num/deno
  return(allele.prop)
}

###############################################
##
## CGF3. CGF recessive mode:Penetrance f0 = f1 and py = 0
##
###############################################

AdmixPropRecCGF <- function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2*(1 - admix.prop^(1/g))
  F1 <- (1 - gamma/2)^g
  F2 <- (1 - gamma)*F1
  F3.1 <- recom.rate*(1 - gamma)*F1; F3.2 <- 0.5*gamma*(1 - recom.rate - 0.5*gamma)*(1 - recom.rate)^g
  F3.3 <- (F3.1 - F3.2)/(recom.rate - gamma/2)
  F3 <- F3.3*(prr - 1)
  F4 <- (1 - gamma)*(F1^2)*(prr - 1)
  num <- F1*(F2*F3 + 1 - 0.5*gamma)
  deno <- (F4 + 1)
  allele.prop <- num/deno
  return(allele.prop)
}

#################################################
##
## CGF4. GCF Dominant mode:Penetrance f1 = f2, px = 1
##
#################################################

AdmixPropDomCGF <- function(generation, recom.rate, admix.prop, prr){
  g <- generation; gamma <- 2*(1 - admix.prop^(1/g))
  F1 <- (1 - gamma/2)^g
  F1.1 <- (1-gamma/2)
  F2 <- 1-(1-gamma)*(F1.1)^(g-1)
  F3.1 <- recom.rate*(1 - gamma)*F1; F3.2 <- 0.5*gamma*(1 - recom.rate - 0.5*gamma)*(1 - recom.rate)^g
  F3.3 <- (F3.1 - F3.2)/(recom.rate - gamma/2)
  F3 <- (F1.1 - F3.3)*(1 - prr)
  F4 <- (1 - F1)
  F5 <- 1-prr - (1-gamma)*F1* (1 - prr)
  num <- F1*(F2*F3 + F1.1*prr)
  deno <- F4*F5 + prr
  allele.prop <- num/deno
  return(allele.prop)
}

##### One generic function to  compute null.prop and alt.prop under Zhu et al. (2004) approach ########
### We need to specify the recom.rate under alternate hypothesis.
### Under Null,  recom.rate = 0.5
### Under alternate, recom.rate < 0.5 
###
#######################################################################

AdmixPropPRR <- function(generation, recom.rate, admix.prop, prr, mode, process){
  if (process == 'HI') {
    if (mode == 'multi') {
      null.prop <- AdmixPropMultiHI(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropMultiHI(generation, recom.rate, admix.prop, prr)
    } else if (mode == 'add') {
      null.prop <- AdmixPropAddHI(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropAddHI(generation, recom.rate, admix.prop, prr)
    } else if (mode == 'rec') {
      null.prop <- AdmixPropRecHI(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropRecHI(generation, recom.rate, admix.prop, prr)
    } else if (mode == 'dom') {
      null.prop <- AdmixPropDomHI(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropDomHI(generation, recom.rate, admix.prop, prr)
    }
  } else if (process == 'CGF') {
    if (mode == 'multi') {
      null.prop <- AdmixPropMultiCGF(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropMultiCGF(generation, recom.rate, admix.prop, prr)
    } else if (mode == 'add') {
      null.prop <- AdmixPropAddCGF(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropAddCGF(generation, recom.rate, admix.prop, prr)
    } else if (mode == 'rec') {
      null.prop <- AdmixPropRecCGF(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropRecCGF(generation, recom.rate, admix.prop, prr)
    } else if (mode == 'dom') {
      null.prop <- AdmixPropDomCGF(generation, 0.5, admix.prop, prr)
      alt.prop <- AdmixPropDomCGF(generation, recom.rate, admix.prop, prr)
    }
  }
  return(c(null.prop, alt.prop))
}

##################################################
## Power and sample size under Zhu et al. (2004)
##################################################

PowerDiscretePRR <- function(generation, recom.rate, admix.prop, prr, case.n, control.n, type1.error, side, study.design, mode, process){
  prop <- AdmixPropPRR(generation, recom.rate, admix.prop, prr, mode, process) 
  if(study.design == 1) {
    power <- PowerCase(prop[1], prop[2], case.n, type1.error, side) 
  } else {
    power <- PowerCaseControl(prop[1], prop[2], case.n, control.n, type1.error, side)
  }
  return(power)
}

SampleDiscretePRR <- function(generation, recom.rate, admix.prop, prr, type1.error, type2.error, side, study.design, mode,process){
  prop <- AdmixPropPRR(generation, recom.rate, admix.prop, prr, mode, process)
  if(study.design == 1) {
    sample.size <- SampleCase(prop[1], prop[2], type1.error,type2.error, side)
  } else {
    sample.size <- SampleCaseControl(prop[1], prop[2], type1.error,type2.error, side)
  }
  return(sample.size)
}		


################ Power and sample size based on ancestry odds ratio:##############
## AOR = ancestry odds ratio
## admix.prop = ancestry proportion
## Multiplicative, case-control study design only
## Fuctions:
## PowerDiscreteAOR
## SampleDiscreteAOR
################################################################################

AdmixPropAOR <- function(admix.prop, aor){
  alt.prop <- (aor*admix.prop)/(1 - admix.prop + aor*admix.prop)
  return(c(admix.prop, alt.prop))
}

PowerDiscreteAOR <- function(admix.prop, aor, case.n, control.n, type1.error, side){
  prop <- AdmixPropAOR(admix.prop, aor)
  power <- PowerCaseControl(prop[1], prop[2], case.n, control.n, type1.error, side)
  return(power)
}

SampleDiscreteAOR <- function(admix.prop, aor, type1.error, type2.error, side){
  prop <- AdmixPropAOR(admix.prop, aor)
  sample.size <- SampleCaseControl(prop[1], prop[2], type1.error,type2.error, side)
  return(sample.size)
}



############## Quantitative traits: Linear regression model: ################
##
## Power and Sample size when the effect size (alpha1 = slope of the model) is known  
## 
## u = locus specific ancesrty proportion  - global ancestry proportion
## ui = excess ancestry at a locus in individual i
## vi = phenotype measurement of ith individual
## Multiple Linear regression model: v = alpha0 + alpha1 *u + zeta*W + e
## W = Covariates; e = random error term, 
## v = (normalized) qunatitative phenotype, u = excess in admixture proportion
## Testing: alpha1 = 0 vs alpha1 != 0
##  
## Need the follwoing information from the population (or data):
## Under H1 --> alpha1 (defined as 'coeff.u' in the function)
## SD of u -->  sigma.u
## SD of error --> sigma.error
## TypeI error --> type1.error
## TypeII error --> type2.error
## r2.covariates = multiple R-square of the variable u with W (from linear regression of u as a dependent of W) 
##	r2.covariates = 0 if no covariates in the model
## Side: side == 1 for one-sided test, == 2 for two-sided test
## 
############ Functions ################
## PowerQTraitCoeff(coeff.u, sigma.u, r2.covariates, sigma.error, type1.error, sample.size, side)
## SampleQTraitCoeff(coeff.u, sigma.u, r2.covariates, sigma.error, type1.error, type2.error, side)
#######################################

PowerQTraitCoeff <- function(coeff.u, sigma.u, r2.covariates, sigma.error, type1.error, sample.size, side){
  inflat <- sample.size*(1 - r2.covariates)		
  se.coeff <- sigma.error/(sigma.u*sqrt(inflat))
  if (side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  power.q <- z.type1.error - coeff.u/se.coeff
  power.qtl <- pnorm(power.q, 0, 1, lower.tail = TRUE, log.p = FALSE)
  return(round(power.qtl, 4))
}

SampleQTraitCoeff <- function(coeff.u, sigma.u, r2.covariates, sigma.error, type1.error, type2.error, side){
  if(side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  z.type2.error <- qnorm(type2.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  num <- sigma.error^2*(z.type1.error + z.type2.error)^2
  deno <- coeff.u^2*sigma.u^2*(1 - r2.covariates)
  sample.size <- ceiling(num/deno)
  return(sample.size)
}



#########################################################
##
## Power and Samle size with No covariates: r-square between v (phenoype) and  u (ancestry) is given
##
############# Functions ###############
## PowerQTraitRSquare(r.square, sample.size, type1.error, side)
## SampleQTraitRSquare(r.square, type1.error, type2.error, side)
##
#######################################

PowerQTraitRSquare <- function(r.square, sample.size, type1.error, side){
  if(side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  delta <- sqrt(sample.size*r.square/(100 - r.square))
  power.q <- z.type1.error - delta
  power <- pnorm(power.q, 0, 1, lower.tail = TRUE, log.p = FALSE)
  return(round(power, 4))
}

SampleQTraitRSquare <- function(r.square, type1.error, type2.error, side){
  if(side == 1) {
    z.type1.error <- qnorm(type1.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  } else {
    z.type1.error <- qnorm(type1.error/2, 0, 1, lower.tail = FALSE, log.p = FALSE)
  }
  z.type2.error <- qnorm(type2.error, 0, 1, lower.tail = FALSE, log.p = FALSE)
  fac1 <- (100 - r.square)/r.square
  fac2 <- (z.type1.error + z.type2.error)^2
  sample.size <- ceiling(fac1*fac2)
  return(sample.size)
}

################################## The End #############################################