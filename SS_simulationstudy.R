# Code to reproduce sample size calculations used in Performance Evaluation (Section 5)
# in 'Bayesian sample size determination using robust commensurate priors with
# interpretable discrepancy weights'

# Lou E. Whitehead 19th December 2023, l.whitehead2@newcastle.ac.uk

################################################################################
# *** Configuration 1 ***
# (relevant parameter values can be altered accordingly to reproduce results for other configs)
################################################################################

################################################################################
# Proposed sample size function borrowing from Q sources
################################################################################

ss.fun <- function(w,tau2){
  
  R = 0.5; # proportion of trial participants allocated to the new treatment
  dw =c(1.1, 1.1); br =c(1e3, 1); # parameters for the Gamma mixture prior
  targEff = 0.4; # target MCID
  eta = 0.95;zeta = 0.80; # posterior decision thresholds
  newtrial_sig02=0.4; # assumed variance in outcomes in the new trial
  
  a=(qnorm(eta) + qnorm(zeta))^2/targEff^2
  b=sum(1/(tau2+w*dw[2]/(dw[1]-1) + (1-w)*br[2]/(br[1]-1)))
  c=newtrial_sig02/(R*(1-R))
  
  n=(a-b)*c # minimum sample size required in the new trial
}

################################################################################
# Config 1 sample size (replace theta_q , tau2_q, w_q accordingly for other configs)
################################################################################

# Historical data, summarized as lambda_q ~ N(theta_q, tau2_q)
theta=c(0.260, 0.240, 0.370, 0.340, 0.320) # Historical means, theta_q
tau2=c(0.250,0.230,0.220,0.360,0.260) # Historical variances , tau^2_q

w=c(0.252,0.319,0.140,0.306,0.149) # w_q values for borrowing

################################################################################
# Process to transform w -> wdash
################################################################################
# predictive precision, g=xi_q^-2 [Equation (14)]
m1=0.1/1.1;m2=1e3/1
g=function(w,tau2,m1,m2){(tau2+w/m1+(1-w)/m2)^-1} 

# h = linear interpolation on g=xi_q^-2 [Equation (15)]
h=function(w, tau2, m1, m2) {g(0,tau2,m1,m2)+w*(g(1,tau2,m1,m2)-g(0,tau2,m1,m2))}

# g_inv = inverse of g=xi_q^2 [Equation (16)]
g_inv=function(p, tau2, m1, m2){((1/p - tau2)*m1*m2 -m1)/(m2-m1)}

# Implementation (Transform w -> wdash)
p=h(w, tau2, m1, m2)
wdash = g_inv(p, tau2, m1, m2)

# Sample size for configuration 1
(n=ss.fun(wdash, tau2))

################################################################################
# Formulation of CPP*, mu_delta ~ N(mu_CPP*, var_CPP*)
# (used for analysis, see Analysis_simulationstudy.R)
################################################################################
# synthesis weights
pstar=g(wdash,tau2,m1,m2)/sum(g(wdash,tau2,m1,m2)) 

# CPP* used for analysis
(mu_CPPwdash=sum(pstar*theta))
(var_CPPwdash=1/sum((tau2+wdash/m1+(1-wdash)/m2)^-1))


