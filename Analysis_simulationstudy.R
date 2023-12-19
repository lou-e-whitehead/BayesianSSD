# Code to reproduce simulation & analysis results in Performance Evaluation (Section 5)
# of 'Bayesian sample size determination using robust commensurate priors with
# interpretable discrepancy weights'

# Lou E. Whitehead 19th December 2023, l.whitehead2@newcastle.ac.uk

################################################################################
# Config 1 settings (calculated in SS_simulationstudy.R) 
# ** NB - replace: n , CPPmu_wdash, CPPvar_wdash for other configs **
################################################################################
n = 38 # sample size for Config 1

# parameters for the commensurate predictive prior for Config 1:
# mu_delta ~ N(CPPmu_wdash*, CPPvar_wdash*)
CPPmu_wdash=0.3078138
CPPvar_wdash=0.06619903

################################################################################
# Generate fake data, theta = 0.4 ** change to 0 to reproduce results for theta = 0 **
################################################################################

true.trt = c(rep(1, n/2), rep(0, n/2)) # treatment assignment
Nobs = n # number of observations
theta.alt = 0.4 # ** change to 0 to reproduce results for theta = 0 **
true.theta = theta.alt # true treatment effect
true.g0 = 0

f=function(true.g0,true.theta,true.trt, Nobs){true.g0 + c(true.theta*true.trt) + rnorm(Nobs, mean = 0, sd = sqrt(0.4))}

set.seed(452)
nsims=10000
y=matrix(nrow = n, ncol = nsims)
for (i in 1:nsims){
  y[,i] <- f(true.g0,true.theta,true.trt, Nobs) # outcome data
}

################################################################################
## Analyse fake data according to Eqns (8) and (9) in manuscript (using CPP(wdash) for config 1)
################################################################################

R=0.5 # proportion of participants randomized to treatment

mean_tx=vector()
mean_ctl=vector()
tx_eff=vector()
var_outcome=vector()
d_theta=vector()
var_theta=vector()
prob_theta_geq_0=vector()
prob_theta_leq_d=vector()

for (i in 1:nsims){
  n=length(y[,i])
  mean_tx[i]=mean(y[,i][which(true.trt==1)]) # mean of observed data treatment group
  mean_ctl[i]=mean(y[,i][which(true.trt==0)]) # mean of observed data control group
  tx_eff[i]=mean_tx[i]-mean_ctl[i] # observed difference in means
  var_outcome[i]=0.4 # variance of outcome (assume known)
  d_theta[i]=(CPPvar_wdash^-1*CPPmu_wdash+tx_eff[i]*n*R*(1-R)/var_outcome[i])/(CPPvar_wdash^-1+n*R*(1-R)/var_outcome[i]) # posterior mean
  var_theta[i]=(CPPvar_wdash^-1+(n*R*(1-R)/var_outcome[i]))^-1 # posterior variance
  prob_theta_geq_0[i]=pnorm(d_theta[i]/sqrt(var_theta[i]))# p(theta>0)
  prob_theta_leq_d[i]=pnorm((0.4-d_theta[i])/sqrt(var_theta[i])) #p(theta<0.4)
  
}

################################################################################
# Verify pre-specified statistical properties are upheld 
# i.e., decision making is possible in 100% of trials using following criteria:
# p(theta>0) >= 0.95 -> Go, p(theta<=delta) >= 0.80 -> No-Go
################################################################################

(a=sum(prob_theta_geq_0>=0.95)/nsims) # p(Go)
(b=sum(prob_theta_leq_d[which(prob_theta_geq_0<0.95)]>=0.80)/nsims) # p(No-Go)
a+b # = 100% as required
