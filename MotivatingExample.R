# Code to reproduce a Motivating Example (Section 4) from manuscript:
# 'Bayesian sample size determination using robust commensurate priors with
# interpretable discrepancy weights'

# Lou E. Whitehead 19th December 2023, l.whitehead2@newcastle.ac.uk

################################################################################
# Motivating Example borrowing from 5 historical sources
################################################################################

################################################################################
# Proposed sample size function borrowing from Q sources
################################################################################

ss.fun <- function(w,tau2){
  
  R = 0.5; # proportion of trial participants allocated to the new treatment
  dw =c(1.1, 1.1); br =c(1e6, 1); # parameters for the Gamma mixture prior
  targEff = 0.1; # target MCID
  eta = 0.95;zeta = 0.80; # posterior decision thresholds
  newtrial_sig02=0.24^2; # assumed variance in outcomes in the new trial
  
  a=(qnorm(eta) + qnorm(zeta))^2/targEff^2
  b=sum(1/(tau2+w*dw[2]/(dw[1]-1) + (1-w)*br[2]/(br[1]-1)))
  c=newtrial_sig02/(R*(1-R))
  
  n=(a-b)*c # minimum sample size required in the new trial
}

################################################################################
# Specify w values for borrowing
w=c(0.15, 0.13, 0.15, 0.17, 0.28)

# Historical variances
tau2=c(0.15^2, 0.10^2, 0.10^2, 0.12^2, 0.13^2)

################################################################################
# Sample size before transformation
################################################################################
(ss_new_trial_before = ss.fun(w=w, tau2=tau2))

################################################################################
# Process to transform w -> wdash
################################################################################
# predictive precision, g=xi_q^-2 [Equation (14)]
m1=0.1/1.1;m2=1e6/1
g=function(w,tau2,m1,m2){(tau2+w/m1+(1-w)/m2)^-1} 

# h = linear interpolation on g=xi_q^-2  [Equation (15)]
h=function(w, tau2, m1, m2) {g(0,tau2,m1,m2)+w*(g(1,tau2,m1,m2)-g(0,tau2,m1,m2))}

# g_inv = inverse of xi_q^2 [Equation (16)]
g_inv=function(p, tau2, m1, m2){((1/p - tau2)*m1*m2 -m1)/(m2-m1)}

# Implementation (Transform w -> wdash)
p=h(w, tau2, m1, m2)
wdash = g_inv(p, tau2, m1, m2)

################################################################################
# Sample size after transformation
################################################################################

(ss_new_trial_after = ss.fun(w=wdash, tau2=tau2))
