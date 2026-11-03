# Power-Lindley Generalized Pareto (PL-GPD) composite distribution
# Density, distribution and quantile functions, with example plots.
# Source: Appendix 2 of N. G. Okoacha, M.Sc. Statistics thesis,
# University of Lagos (supervisor: Dr E. B. Nkemnole). Code reproduced unchanged.
# Parameters: u = threshold, a = alpha and b = beta (Power Lindley body),
#             r = gamma (GPD shape), c = GPD scale.
#
#=====================================================================
#=====================================================================
library(fitdistrplus)#for MLE estimation
library(ggplot2)#for generating plots
library(KScorrect)#for calculating KS and its P-Value
library(mc2d)# support ggplot
library(rgenoud)#support ggplot
library(extraDistr)# contains some univariate distributions
library(Matrix)# for evaluation of matrixe
library(numDeriv)#for numerical derivatives
library(fBasics) #for calculating summary statistics
library(sn)# contains the skew families
library(EnvStats)# contains some univariate distribution
library(stats)
library(LambertW)
library(lamW)
library (tseries)
#===============================================================================
####################### F1-F2-Composite Model ########################### 
#=========================================================================
#The Generalized Pareto Distribution
#γ=r #α
dGPD<-function(x,u,r,c) (1/c)*(1+(r*(x-u)/c))^(-1-(1/r))
pGPD<-function(q,u,r,c) 1-(1+(r*(q-u)/c))^(-1/r)
qGPD<-function(p,u,r,c) (c/r)*((1-p)^(-r)-1)+u

#The Power Lindley Distribution
#α=a,β=b
dPLD<-function(x,a,b) ((a*b^2)/(b+1))*(1+x^a)*x^(a-1)*exp(-b*x^a)
pPLD<-function(q,a,b) 1-(1+((b*q^a)/(b+1)))*exp(-b*q^a)
qPLD<-function(p,a,b) (-1-(1/b)-(1/b)*lambertWm1(-(1+b)*
                                    exp(-(1+b))*(1-p)))^(1/a)

#==============================================================================
#derivative of the Power Lindley density
#r=r,tau=mu,phi=sigma, alpha=alpha
dPLDprime<-function(x,a,b){
  dPLD(x,a,b)*( ((a-1)/x)+((a*x^a)/(1+x^a))-(a*b*x^(a-1)) )

}
#======================================================================
#====================================================================== 
#============Power Lindley- Generalized Pareto Distribution ==========
#=======================================================================
#=======================================================================
#PLGPD pdf 
dPLGPD<-function(x,u,a,b,r,c){ 
  d1=rep(0,length(x))
  
  w<- (pPLD(u,a,b)*dGPD(u,u,r,c))/(dPLD(u,a,b)+dGPD(u,u,r,c)*pPLD(u,a,b))
  
  #c<- (-(1-w)*(1+r)*pPLD(u,a,b)/w*dPLDprime(u,a,b))^(1/2)
  d1[which(x<=u)]<-w*dPLD(x[which(x<=u)],a,b)/pPLD(u,a,b) 
  d1[which(x>=u)]<-(1-w)*dGPD(x[which(x>=u)],u,r,c) 
  return(d1) 
} 
#=========================================================================
#PLGPD cdf 
pPLGPD<-function(q,u,a,b,r,c){ 
  p1=rep(0,length(q)) 
  
  w<- (pPLD(u,a,b)*dGPD(u,u,r,c))/(dPLD(u,a,b)+dGPD(u,u,r,c)*pPLD(u,a,b))
  
  #c<- (-(1-w)*(1+r)*pPLD(u,a,b)/w*dPLDprime(u,a,b))^(1/2)
  
  
  p1[which(q<=u)]<-w*pPLD(q[which(q<=u)],a,b)/pPLD(u,a,b)
  p1[which(q>=u)]<-w+(1-w)*pGPD(q[which(q>=u)],u,r,c)
  return(p1) 
} 
#============================================================
#PLGPD quantile function 
qPLGPD<-function(p,u,a,b,r,c){ 
  q1=rep(0,length(p))
  
  w<- (pPLD(u,a,b)*dGPD(u,u,r,c))/(dPLD(u,a,b)+dGPD(u,u,r,c)*pPLD(u,a,b))
  
  #c<- (-(1-w)*(1+r)*pPLD(u,a,b)/w*dPLDprime(u,a,b))^(1/2)
  
  
  q1[which(p<=w)]<-qPLD(p[which(p<=w)]*pPLD(u,a,b)/w,a,b) 
  q1[which(p>w)]<-qGPD((p[which(p>w)]-w)/(1-w),u,r,c)
  return(q1) 
} 
#==============================================================
#PLGPD PDF
x<-seq(0,5, 0.001)

y1<-dPLGPD(x,4.5,0.5,1,0.5,1.5)
y2<-dPLGPD(x,4,1.5,1,0.8,1.5)
y3<-dPLGPD(x,3.5,2.5,1,2,1.5)
y4<-dPLGPD(x,3,4,1,2.5,1.5)

plot(c(0,5),c(0,1.5),type="n",ylab="f(x)",xlab="x
     β=1, c=1.5")
lines(x,y1,col=1,lty=1)
lines(x,y2,col=2,lty=3)
lines(x,y3,col=3,lty=5)
lines(x,y4,col=4,lty=7)

legend("topright",c(
  "α=0.5 u=4.5 γ=0.5",
  "α=1.5 u=4    γ=0.8",
  "α=2.5 u=3.5 γ=2",
  "α=4    u=3    γ=2.5"),
  col=1:4, bty="n",lty=c(1,3,5,7))
#==============================================================

#==============================================================
#PLGPD CDF
x<-seq(0,5, 0.001)

y1<-pPLGPD(x,4.5,0.5,1,0.5,1.5)
y2<-pPLGPD(x,4,1.5,1,0.8,1.5)
y3<-pPLGPD(x,3.5,2.5,1,2,1.5)
y4<-pPLGPD(x,3,4,1,2.5,1.5)

plot(c(0,5),c(0,1),type="n",ylab="F(x)",xlab="x
     β=1, c=1.5")
lines(x,y1,col=1,lty=1)
lines(x,y2,col=2,lty=3)
lines(x,y3,col=3,lty=5)
lines(x,y4,col=4,lty=7)

legend("bottomright",c(
  "α=0.5 u=4.5 γ=0.5",
  "α=1.5 u=4    γ=0.8",
  "α=2.5 u=3.5 γ=2",
  "α=4    u=3    γ=2.5"),
  col=1:4, bty="n",lty=c(1,3,5,7))
#==============================================================
#===========================================================
