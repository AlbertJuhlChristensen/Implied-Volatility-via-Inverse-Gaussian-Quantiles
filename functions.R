library(pracma)
library(ggplot2)
library(statmod)
library(scales)
library(dplyr)
library(knitr)
library(tidyr)
BS_call_sim = function(S_t, t,K,BigT,r,delta,sigma){
  d1 = (log(S_t/K)+(r-delta+sigma^2/2) *(BigT-t))/(sigma*sqrt(BigT-t))
  d2 = d1 - sigma*sqrt(BigT-t)
  call_price = S_t*exp(-delta*(BigT-t))*pnorm(d1)-exp(-r*(BigT-t))*K*pnorm(d2)
  return(call_price)
}

BS_put_sim = function(S_t, t,K,BigT,r,delta,sigma){
  d1 = (log(S_t/K)+(r-delta+sigma^2/2) *(BigT-t))/(sigma*sqrt(BigT-t))
  d2 = d1 - sigma*sqrt(BigT-t)
  put_price = K*exp(-r*(BigT-t))*pnorm(-d2)-exp(-delta*(BigT-t))*S_t*pnorm(-d1)
  return(put_price)
}

IV_func = function(sigma,obs,S_t,t,K,BigT,r,delta,opttype){
  if (opttype=="call"){
    price = BS_call_sim(S_t,t,K,BigT,r,delta,sigma)
  }
  else {
    price = BS_put_sim(S_t,t,K,BigT,r,delta,sigma)
  }
  return(price-obs)
}

IV_uniroot = function(obs,S_t,t,K,BigT,r,delta,opttype){
  sigma_hat = uniroot(IV_func, interval = c(0.001,10), obs=obs, S_t=S_t,
                      t=t,K=K,BigT=BigT,r=r,delta=delta,opttype=opttype)$root
  return(sigma_hat)
}

Vega = function(sigma,S_t,t,K,BigT,r,delta,...){
  d1 = (log(S_t/K)+(r-delta+sigma^2/2)*(BigT-t)) /
    (sigma*sqrt(BigT-t))
  
  vega = S_t*exp(-delta*(BigT-t))*dnorm(d1)*sqrt(BigT-t)
  return(vega)
}

IV_NR = function(obs,S_t,t,K,BigT,r,delta,opttype){
  sigma0 = runif(1,0.05,0.55)
  sigma_hat = newtonRaphson(IV_func,x0=sigma0,dfun=Vega,maxiter = 500,
                            obs=obs, S_t=S_t, t=t,K=K,BigT=BigT,
                            r=r,delta=delta,opttype=opttype, tol=1e-14)$root
  return(sigma_hat)
}

IV_IG = function(obs,S_t,t,K,BigT,r,delta,opttype){
  Forward=exp((r-delta)*(BigT-t))*S_t
  D = exp(-r*(BigT-t))
  k=log(K/Forward)
  m=pmin(1,exp(k))
  opt=obs/(D*Forward)
  if (opttype=="call"){
  sigma_hat = 2/(sqrt(BigT-t)) * 1 /sqrt(qinvgauss(p=(1-opt)/m,mean=2/abs(k),shape=1))
  }
  else {
    sigma_hat = 2/(sqrt(BigT-t)) * 1 /sqrt(qinvgauss(p=(K/Forward-opt)/m,
                                              mean=2/abs(k),shape=1))
  }
  return(sigma_hat)
}


