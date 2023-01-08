library(ggamma)
library(SuppDists)
library(caret)
library(hydroGOF)
library(quantreg)

ga<-function(theta,data){
  sd<-theta[1]
  beta<-theta[2]
  n<-length(data)
  return(-n*log(beta^sd*gamma(sd))+(sd-1)*sum(log(data))-beta^(-1)*sum(data))
}
gammaqujian<-function(p,pre_mu,pre_sig){
  #shape<-(pre_mu)^2/(pre_sig)^2
  #rate<-pre_mu/(pre_sig)^2
  shape<-1/(pre_sig)^2;
  rate<-1/(pre_mu*pre_sig^2);
  feiweishu_up<-qgamma(p,shape,rate,lower.tail = FALSE);
  feiweishu_down<-qgamma(p,shape,rate,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
gouzaogamma<-function(n,pre_mu,pre_sig){
  shape<-1/(pre_sig)^2;
  rate<-1/(pre_mu*pre_sig^2);
  p<-runif(n)
  ngamma<-qgamma(p,shape,rate,lower.tail = TRUE)
  #ngamma<-rgamma(n,shape,rate)
  return(ngamma)
}
p<-runif(100)
levenejianyan<-function(ngamma1,ngamma2){
  p<-leveneTest(ngamma1~factor(ngamma2),center=median)
  return(p)
}
ggammaqujian<-function(p,pre_mu,pre_sig,pre_nu){
  a<-1/(pre_sig^2);
  b<-pre_nu;
  d<-1/(pre_mu*(pre_sig^2));
  k<-d/b;
  #shape<-1/(pre_sig)^2;
  #rate<-1/(pre_mu*pre_sig^2);
  feiweishu_up<-qggamma(p,a,b,k,lower.tail = FALSE);
  feiweishu_down<-qggamma(p,a,b,k,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
igqujian<-function(p,pre_mu,pre_sig){
  lambda<-(pre_mu^3)/(pre_sig^2);
  feiweishu_up<-qinvGauss(p,pre_mu,lambda,lower.tail = FALSE);
  feiweishu_down<-qinvGauss(p,pre_mu,lambda,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
expqujian<-function(p,pre_mu){
  lambda<-1/pre_mu;
  feiweishu_up<-qexp(p,lambda,lower.tail = FALSE);
  feiweishu_down<-qexp(p,lambda,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
lnoqujian<-function(p,pre_mu,pre_sig){
  omiga<-exp(pre_sig^2);
  ppred_mu<-exp(pre_mu)*omiga^0.5;
  ppred_sig<-exp(pre_mu*2)*omiga*(omiga-1);
  feiweishu_up<-ppred_mu+1.96*ppred_sig^0.5;
  feiweishu_down<-ppred_mu-1.96*ppred_sig^0.5;
  return(list(Down=feiweishu_down,Mu=ppred_mu,Up=feiweishu_up))
}
noqujian<-function(p,pre_mu,pre_sig){
  feiweishu_up<-qnorm(p,pre_mu,pre_sig,lower.tail = FALSE);
  feiweishu_down<-qnorm(p,pre_mu,pre_sig,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
gblqujian<-function(p,pre_mu,pre_sig){
  feiweishu_up<-qgumbel(p,pre_mu,pre_sig,lower.tail = FALSE);
  feiweishu_down<-qgumbel(p,pre_mu,pre_sig,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
loqujian<-function(p,pre_mu,pre_sig){
  pre_sig<-sqrt(3*pre_sig)/pi
  feiweishu_up<-qlogis(p,pre_mu,pre_sig,lower.tail = FALSE);
  feiweishu_down<-qlogis(p,pre_mu,pre_sig,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
tfqujian<-function(p,pre_mu,pre_sig){
  feiweishu_up<-qt(p,pre_mu,pre_sig,lower.tail = FALSE);
  feiweishu_down<-qt(p,pre_mu,pre_sig,lower.tail = TRUE);
  return(list(Down=feiweishu_down,Up=feiweishu_up))
}
