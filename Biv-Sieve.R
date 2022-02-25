##-------------------------------------------------------------------------------------
## R code to conduct simulation experiment in the paper:
## "Sieve estimation of the additive hazards model with bivariate current status data"
##  by 
## 
## Ce Zhang, Riyadh Rustam Al-Mosawi, Dipankar Bandyopadhyay, Haiwu Huang and Xuewen Lu
## Author: Riyadh R. Al-Mosawi, PhD
## Edited by: Dipankar Bandyopadhyay, PhD
## Date: 02/25/2022
##-------------------------------------------------------------------------------------


rm(list=ls())

library(numDeriv)     # To compute the numerical derivatives
library(alabama)      # To use constrOptim.nl function
library(copula)       # To compute and generate different copula models



options(width=100,length=200,max.print = 10000,digits = 3,scipen=999)


### -------- User-Defined Functions

log.like=function(para){ 
  # the negative of log-likelihood function
  beta.para=para[beta.pos] 
  b1.para  =para[b1.pos] 
  b2.para  =para[b2.pos] 
  phi.para =para[phi.pos]
  
  # marginal cumulative hazard functions of T1 and T2
  min.s=1e-10; max.s=1-min.s
  
  # marginal survival function of T1 
  B1=apply(as.matrix(0:q),1,function(y) bern(y,q,l1,u1,c1)) 
  s1=exp(-B1%*%b1.para-as.matrix(X)%*%beta.para*c1) 
  s1=s1*(s1>=min.s)*(s1<=max.s)+min.s*(s1<min.s)+max.s*(s1>max.s)
  
  # marginal survival function of T2 
  B2=apply(as.matrix(0:q),1,function(y) bern(y,q,l2,u2,c2))   
  s2=exp(-B2%*%b2.para-as.matrix(X)%*%beta.para*c2)
  s2=s2*(s2>=min.s)*(s2<=max.s)+min.s*(s2<min.s)+max.s*(s2>max.s)
  
  k1=apply(as.matrix(0:q),1,function(x) bern(x,q,0,1,s1))
  k2=apply(as.matrix(0:q),1,function(x) bern(x,q,0,1,s2))
  
  s11=phi.val(k1,k2,phi.para,q)
  
  s11=pmax(s11,min.s)
  s01=pmax(s2-s11,min.s) 
  s10=pmax(s1-s11,min.s) 
  s00=pmax(1-s1-s2+s11,min.s)
  l.like=sum(del1*del2*log(s11)+(1-del1)*del2*log(s01)+del1*(1-del2)*log(s10)+(1-del1)*(1-del2)*log(s00))
  return(-l.like)
}


bern=function(j,p,l,u,t){       # Bernstein polynomial
  return(choose(p,j)*((t-l)/(u-l))^j*(1-(t-l)/(u-l))^(p-j))
}

hin.b=function(y){     # set of constraints on b1 (or b2) of function H1( or H2)
  return(c(y[1],diff(y)))
}

hin.b1b2=function(y){ # set of constraints on b1 and b2 of function H1 and H2
  return(c(hin.b(y[1:(q+1)]),hin.b(y[(q+2):(2*(q+1))])))
}

hin.phi=function(y){    # set of constraints on phi of copula function
  if(q==2) return(c(y[1],1/2-y[1]))  
  if(q==3) return(c(y[1],1/3-y[2],1/3-y[3],2/3-y[4],y[2]-y[1],y[3]-y[1],y[4]-y[2],y[4]-y[3]))
  if(q==4) return(c(y[1],1/4-y[3],2/4-y[6],1/4-y[7],2/4-y[8],3/4-y[9],y[2]-y[1],
          y[3]-y[2],y[5]-y[4],y[6]-y[5],y[8]-y[7],y[9]-y[8],y[4]-y[1],y[7]-y[4],y[5]-y[2],y[8]-y[5],y[6]-y[3],y[9]-y[6])) 
}

hin.b1b2phi=function(y){     # set of constraints on b1,b2 and phi
  return(c(hin.b(y[1:(q+1)]),hin.b(y[(q+2):(2*(q+1))]),hin.phi(y[(2*(q+1)+1):length(y)])))
}


phi.val=function(k1=k1,k2=k2,y,q=q){  # Compute the value of the function phi 
  if(q==2){
    fit.val=y[1]*k1[,2]*k2[,2]+(1/2)*k1[,2]*k2[,3]+(1/2)*k1[,3]*k2[,2]+k1[,3]*k2[,3]
  }
  if(q==3){
    fit.val=y[1]*k1[,2]*k2[,2]+y[2]*k1[,2]*k2[,3]+y[3]*k1[,3]*k2[,2]+y[4]*k1[,3]*k2[,3]
    fit.val=fit.val+(1/3)*(k1[,2]*k2[,4]+k1[,4]*k2[,2])+(2/3)*(k1[,3]*k2[,4]+k1[,4]*k2[,3])+k1[,4]*k2[,4]
  }
  if(q==4){
    fit.val=y[1]*k1[,2]*k2[,2]+y[2]*k1[,2]*k2[,3]+y[3]*k1[,2]*k2[,4]
    fit.val=fit.val+y[4]*k1[,3]*k2[,2]+y[5]*k1[,3]*k2[,3]+y[6]*k1[,3]*k2[,4]
    fit.val=fit.val+y[7]*k1[,4]*k2[,2]+y[8]*k1[,4]*k2[,3]+y[9]*k1[,4]*k2[,4]    
    fit.val=fit.val+(1/4)*k1[,2]*k2[,5]+(1/4)*k1[,5]*k2[,2]+(2/4)*k1[,3]*k2[,5]
    fit.val=fit.val+(2/4)*k1[,5]*k2[,3]+(3/4)*k1[,4]*k2[,5]
    fit.val=fit.val+(3/4)*k1[,5]*k2[,4]+k1[,5]*k2[,5]
  }
  return(fit.val)  
}


Gen.Data=function(n,pr,mn,std,a,l1,u1,l2,u2){ # Data generation function
  # generate Covariates
  X1=rbinom(n,1,pr)             # Bernoulli covariate with probability=pr 
  X2=rnorm(n,mean=mn,sd=std)    # normal covariate with mean=mn and sigma=std

  ## generate bivariate failure times using copula package
  w=cop.gen(n) 
  u=w[,1]
  v=w[,2] 

  x.beta=cbind(X1,X2)%*%beta.true
  T1=-log(u)/(1+x.beta)                    # derived from s(t)=exp(-0.3t-x.beta*t)
  T2=-(x.beta-(x.beta^2-4*log(v))^(1/2))/2 # derived from s(t)=exp(-t^2-x.beta*t)
  
  c1=runif(n,l1,u1)         # The censoring times if the T1
  c2=runif(n,l2,u2)         # The censoring times if the T2
  Lam1=Lambda1(c1)          # The true cumulative hazards function H1
  Lam2=Lambda2(c2)          # The true cumulative hazards function H2
  
  #-- Generating the survival times
  del1=as.numeric(I(T1>=c1))     # censoring indicator of the first failure time, T1
  del2=as.numeric(I(T2>=c2))     # censoring indicator of the second failure time, T2
  
  Data=list(X1=X1,X2=X2,c1=c1,c2=c2,del1=del1,del2=del2,Lam1=Lam1,Lam2=Lam2)
  return(Data)
}


Init.Val=function(b1.0,b2.0,phi.0){  # function to obtain initial values of the Bernstein coefficients
  B1=apply(as.matrix(0:q),1,function(y) bern(y,q,l1,u1,c1))
  B2=apply(as.matrix(0:q),1,function(y) bern(y,q,l2,u2,c2))   # Bernstein polynomials 
  
  s1=exp(-Lam1.true-x.beta*c1)  
  s2=exp(-Lam2.true-x.beta*c2) # true survival functions
  
  s=cbind(s1,s2)
  cop.true=as.matrix(apply(as.matrix(s),1,copula))
  
  #-- Computing the initial values based on least square method
  
  #-- Initial values of the Bernstein coefficients, b1
  b1.ini=try(auglag(par=b1.0,fn=function(y) sum((Lam1.true-B1%*%y)^2),
        hin=hin.b,control.outer=list(eps=1e-10,itmax=500000,trace=FALSE))$par,silent = TRUE)
  if(is.character(b1.ini)) b1.ini=b1.0
  
  #-- Initial values of the Bernstein coefficients, b2
  b2.ini<-try(auglag(par=b2.0,fn=function(y) sum((Lam2.true-B2%*%y)^2),
       hin=hin.b,control.outer=list(eps=1e-10,itmax=500000,trace=FALSE))$par,silent = TRUE)
  if(is.character(b2.ini)) b2.ini=b2.0
  
  #-- Initial values of the Bernstein coefficients, phi
  
  k1=apply(as.matrix(0:q),1,function(y) bern(y,q,0,1,s1))
  k2=apply(as.matrix(0:q),1,function(y) bern(y,q,0,1,s2))
  
  phi.ini=try(auglag(par=phi.0,fn=function(w) sum((cop.true-phi.val(k1,k1,w,q))^2),
        hin=hin.phi,control.outer=list(eps=1e-10,itmax=500000,trace=FALSE))$par,silent = TRUE)
  if(is.character(phi.ini)) phi.ini=phi.0
  
  return(list(b1.ini=b1.ini,b2.ini=b2.ini,phi.ini=phi.ini))
}

Copula.fun=function(name,tau=tau){
  #----------------------------------------------
  # Function to define the copula function
  # a : association parameter 
  # tau: tau of the Kendall's tau measure
  #----------------------------------------------
  copula.name=copula=cop.gen=a=NULL
  
  if(name=="clayton"){
    ###------------ Clayton copula
    copula.name="clayton copula"
    copula     =function(y) pCopula(y,claytonCopula(param=a, dim=2))
    cop.gen    =function(m) rCopula(m, claytonCopula(param=a, dim=2))
    if(tau==0.25) a=0.666
    if(tau==0.50) a=2
    #tau= tau(claytonCopula(param=a, dim=2))
  }
  
  if(name=="fgm"){
    copula.name="FGM copula"
    copula     =function(y) pCopula(y,fgmCopula(param=a, dim=2))
    cop.gen    =function(m) rCopula(m, fgmCopula(param=a, dim=2))
    if(tau==0.22) a=1 
  }
  
  ###------------ Gaussian copula
  if(name=="gaussian"){
    copula.name="Gaussian copula"
    copula     =function(y) pCopula(y,normalCopula(param=a, dim=2))
    cop.gen    =function(m) rCopula(m, normalCopula(param=a, dim=2))
    if(tau==0.25) a=0.383
    if(tau==0.50) a=0.707
  }
  
  if(name=="gumbel"){
    ###------------ Gumbel copula
    copula.name="Gumbel copula"
    copula     =function(y) pCopula(y,gumbelCopula(param=a, dim=2))
    cop.gen    =function(m) rCopula(m, gumbelCopula(param=a, dim=2))
    if(tau==0.25) a=1.33
    if(tau==0.50) a=2
  }
  
  if(name=="frank"){
    ###------------ Frank copula
    copula.name="Frank copula"
    copula     =function(y) pCopula(y,frankCopula(param=a, dim=2))
    cop.gen    =function(m) rCopula(m, frankCopula(param=a, dim=2))
    if(tau==0.25) a=2.375 
    if(tau==0.50) a=5.74
  }
  return(list(copula.name=copula.name,copula=copula,copula.gen=cop.gen,a=a))
} 


#---------------------------------------------
#--------- The Main Code ---------------------
#---------------------------------------------

q=3           # degree of Bernstein polynomials 3
n=100         # sample size
            # lower and upper bounds 
l1=0.06     # the lower bound of the censoring times C1
u1=1        # the upper bound of the censoring times C1
l2=0.06     # the lower bound of the censoring times C2
u2=1        # the upper bound of the censoring times C2

Sim.No=50 #500    # number of replications 

beta.true=c(0.5,0.9)    # true value of beta


# true cumulative hazards functions
Lambda1=function(y) y        
Lambda2=function(y) y^2



#----- defining the position of the parameters

beta.dim=length(beta.true)
b1.dim  =q+1
b2.dim  =q+1
phi.dim =(q-1)^2
para.dim=beta.dim+b1.dim+b2.dim+phi.dim     
beta.pos=1:beta.dim
b1.pos  =(beta.dim+1):(beta.dim+b1.dim)
b2.pos  =(beta.dim+b1.dim+1):(beta.dim+b1.dim+b2.dim)
phi.pos =(beta.dim+b1.dim+b2.dim+1):para.dim


## Defining the initial values of beta, b1,b2 and phi
beta.0=rep(0,beta.dim)
b1.0  =seq(0.2,50,length=q+1) 
b2.0  =seq(20,30,length=q+1)
phi.0 =rep(0,phi.dim) 
ij    =0
for(i in 1:(q-1))
  for(j in 1:(q-1)){
    ij=ij+1
    if(ij==1) phi.0[ij]=0.04 else phi.0[ij]=phi.0[ij-1]+0.03
  }

## These codes are useful for plotting the curves
par(mfrow=c(1,2))
c.axis.1=seq(l1,u1,(u1-l1)/100)
c.axis.2=seq(l2,u2,(u2-l2)/100)
Lambda1.true=Lambda1(c.axis.1)      # true cumulative hazards function1
Lambda2.true=Lambda2(c.axis.2)      # true cumulative hazards function2
BB1=apply(as.matrix(0:q),1,function(x) bern(x,q,l1,u1,c.axis.1)) 
BB2=apply(as.matrix(0:q),1,function(x) bern(x,q,l2,u2,c.axis.2)) 


Copula.Name="clayton"   # You can use "frank" or "gumbel" or "gaussian" or"fgm",
Copula.tau=0.25        # You can use 0.5

copula.inf=Copula.fun(name=Copula.Name,tau=Copula.tau)
copula.name=copula.inf[[1]]
copula     =copula.inf[[2]]
cop.gen    =copula.inf[[3]]
a          =copula.inf[[4]]

#-----------------------------------------
#--- Main code implementation
#-----------------------------------------
for(i in 1:Sim.No) {
  cat("\n--Results of the ", i,"th Simulation-------Results upto the ", i,"th Simulation", "\n")
  # save results of each simulation
  if(i==1){
    Start.date=date()
    Beta.Est=b1.Est=b2.Est=phi.Est=CP.Prof=ESD.Prof=SSD=Results=NULL 
    set.seed(1234)
  }
  #-- data Generation Part
  Data=suppressWarnings(try(Gen.Data(n,pr=0.6,mn=1.2,std=0.6,a,l1,u1,l2,u2),silent=TRUE))
  if(is.character(Data)) next
  X1=Data$X1
  X2=Data$X2
  X=cbind(X1,X2) 
  c1=Data$c1
  c2=Data$c2     
  del1=Data$del1
  del2=Data$del2
  Lam1.true=Data$Lam1
  Lam2.true=Data$Lam2
  x.beta=X%*%beta.true
  
  B1=apply(as.matrix(0:q),1,function(y) bern(y,q,l1,u1,c1))
  B2=apply(as.matrix(0:q),1,function(y) bern(y,q,l2,u2,c2))  
  
  #-- Compting the initial values based on least square method
  init.val=Init.Val(b1.0,b2.0,phi.0)
  b1.ini  =init.val$b1.ini
  b2.ini  =init.val$b2.ini
  phi.ini =init.val$phi.ini
  beta.ini=beta.true
  #-- Computing the MLE 
  est.para.b1b2=try(auglag(par=c(b1.ini,b2.ini),fn=function(y) log.like(c(beta.ini,y,phi.ini)),
                           hin=hin.b1b2,control.outer=list(eps=1e-10,itmax=5000000,trace=FALSE)))
  if(is.character(est.para.b1b2)) next
  
  b1.u =est.para.b1b2$par[b1.pos-beta.dim]
  b2.u =est.para.b1b2$par[b2.pos-beta.dim]
  
  est.para.phi=try(auglag(par=phi.ini,fn=function(y) log.like(c(beta.ini,b1.u,b2.u,y)),
                          hin=hin.phi,control.outer=list(eps=1e-10,itmax=5000000,trace=FALSE)))
  if(is.character(est.para.phi)) next
  
  phi.u=est.para.phi$par
  
  est.para.beta=try(optim(par=beta.ini,fn=function(y) log.like(c(y,b1.u,b2.u,phi.u)),
                          method="L-BFGS-B",control=list(maxit=5000000)),silent = TRUE)
  if(is.character(est.para.beta)) next
  
  beta.u=est.para.beta$par 
  
  # Computing SE using profile information matrix and observed information matrix 
  Beta.Est=rbind(Beta.Est,beta.u)
  b1.Est  =rbind(b1.Est,b1.u) 
  b2.Est  =rbind(b2.Est,b2.u) 
  phi.Est =rbind(phi.Est,phi.u)
  
  beta.se.prof =beta.se.obse=numeric(beta.dim)
  beta.cov.prof=beta.cov.obse=rep(-1,beta.dim)
  
  prof.inf=suppressWarnings(try(sqrt(diag(solve(hessian(function(y) log.like(c(y,b1.u,b2.u,phi.u)), beta.u)))),silent=TRUE))
  if(!is.character(prof.inf))
    if(!any(is.nan(prof.inf)))
      if(!any(is.na(prof.inf))){
        beta.se.prof =prof.inf
        beta.cov.prof=as.numeric((beta.true>=beta.u-1.96*beta.se.prof)&(beta.true<=beta.u+1.96*beta.se.prof))
        ESD.Prof     =rbind(ESD.Prof,beta.se.prof) 
        CP.Prof      =rbind(CP.Prof,beta.cov.prof) 
      }
  
  
  cum.beta.u=cum.beta.ssd=cum.beta.se.prof=cum.beta.se.obse=numeric(beta.dim)
  cum.beta.cov.prof=cum.beta.cov.obse=rep(-1,beta.dim)
  #--- cp=-1 means se does not converge, 
  #--- cp=0 means the true beta lies outside the confidence interval
  #--- cp=1 means the true beta lies inside the confidence interval
  if(!is.null(Beta.Est)) if(length(Beta.Est[,1])==1) cum.beta.u=beta.u else cum.beta.u=round(apply(Beta.Est,2,mean),3)
  if(!is.null(Beta.Est)) if(length(Beta.Est[,1])==1) cum.beta.ssd=numeric(beta.dim) else cum.beta.ssd=round(apply(Beta.Est,2,sd),3)
  if(!is.null(ESD.Prof)) if(length(ESD.Prof[,1])==1) cum.beta.se.prof=beta.se.prof else cum.beta.se.prof=round(apply(ESD.Prof,2,mean),3)
  if(!is.null(CP.Prof)) if(length(CP.Prof[,1])==1) cum.beta.cov.prof=beta.cov.prof*100 else cum.beta.cov.prof=apply(CP.Prof,2,mean)*100
  Estim.Value=NULL
  Estim.Value=rbind(Estim.Value,c(beta.u,cum.beta.u))
  Estim.Value=rbind(Estim.Value,c(beta.u-beta.true,cum.beta.u-beta.true))
  Estim.Value=rbind(Estim.Value,c(rep(0,beta.dim),cum.beta.ssd))
  Estim.Value=rbind(Estim.Value,c(beta.se.prof,cum.beta.se.prof))
  Estim.Value=rbind(Estim.Value,c(beta.cov.prof,cum.beta.cov.prof))
  Est.Val=data.frame(Estim.Value)
  colnames(Est.Val)=c("Beta1","Beta2","tot. Beta1","tot. Beta2")
  row.names(Est.Val)= c("Estimates of Beta                                   =",
                        "Estimates of Bias                                   =",
                        "Sample standrad error (SSE)                         =",
                        "Estimates of se using profile  likelihood           =",
                        "Estimates of CP using profile likelihood            =")
  print(round(Est.Val,3))
  
  #-- Plots of the estimated curves of Lambda1 and Lambda2 till the i-th iteration
  if(i==1) next
  
  Lambda1.median=BB1%*%apply(b1.Est,2,median)# fitted cumulative hazards function using medians
  
  matplot(c.axis.1, cbind(Lambda1.true,Lambda1.median), type = "l",main=expression(paste("Fitted ",hat(H)[1] (c))),
          xlab = "c", ylab=expression(paste("baseline cumulative hazards function, ", H[1])),
          lwd = 2, lty =c(1,2), col = c("black","red"),xlim=c(0,1),ylim=c(0,1.2))
  legend("topleft", c(expression(paste("fitted ",hat(H)[1])),expression(paste("true ",hat(H)[1]))),lty =c(1,2), lwd = 1, col =c("black","red"), inset = 0.00001)
  
  Lambda2.median=BB2%*%apply(b2.Est,2,median) # fitted cumulative hazards function using medians
  matplot(c.axis.2, cbind(Lambda2.true,Lambda2.median), type = "l",main=expression(paste("Fitted ",hat(H)[2] (c))),
          xlab = "c", ylab=expression(paste("baseline cumulative hazards function, ", H[2])),
          lwd = 2, lty =c(1,2), col = c("black","red"),xlim=c(0,1),ylim=c(0,1.2))
  legend("topleft", c(expression(paste("fitted ",hat(H)[2])),expression(paste("true ",hat(H)[2]))),lty =c(1,2), lwd = 1, col =c("black","red"), inset = 0.00001)
  
  
  if(i==Sim.No){
    End.date=date()
    ##---- Computing the estimators
    MLE.Beta.Est =round(apply(Beta.Est,2,mean),3)
    MLE.Beta.Bias=round(MLE.Beta.Est-beta.true,3)
    MLE.Beta.SSD =round(apply(Beta.Est,2,sd),3)
    MLE.Beta.ESD.Prof=round(apply(ESD.Prof,2,mean),3)
    Cov.Pro.Prof=round(apply(CP.Prof,2,mean),3)*100
    
    Results=NULL
    Results=rbind(Results,beta.true)
    Results=rbind(Results,MLE.Beta.Est)
    Results=rbind(Results,MLE.Beta.Bias)
    Results=rbind(Results,MLE.Beta.SSD)
    Results=rbind(Results,MLE.Beta.ESD.Prof)
    Results=rbind(Results,Cov.Pro.Prof)
    RES=data.frame(Results)
    colnames(RES)=c("Beta1   ","Beta2   ")
    row.names(RES)=c("True Beta","Beta Estimates","Avg. Bias","sampled standard error (SSE)","profile ESE","profile 95% CP")
    print(round(RES,4))
    
    MLE.beta=apply(Beta.Est[,1:beta.dim],2,mean)
    MLE.b1=apply(b1.Est,2,mean)
    MLE.b2=apply(b2.Est,2,mean)
    MLE.phi=apply(phi.Est,2,mean)
    MLE.Est=c(MLE.beta,MLE.b1,MLE.b2,MLE.phi)
    BIC=log.like(MLE.Est)+2*log(n)*para.dim #-- compute the BIC value
    
    cat("copula type                    : ",Copula.Name, "\n")
    cat("degree of Bernstein polynomial : ",q, "\n")
    cat("b1 estimates                   : ",MLE.b1, "\n")
    cat("b2 estimates                   : ",MLE.b2, "\n")
    cat("phi estimates                  : ",MLE.phi, "\n")
    cat("BIC                            : ",BIC,"\n")
    cat("interval of generating c1      : ","(",l1,",",u1,")\n")  
    cat("interval of generating c2      : ","(",l2,",",u2,")\n")  
    cat("sample size (n)                : ",n,"\n")
    cat("no. of simulations (Sim.No)    : ",Sim.No,"\n")
    cat("true association parameter (a) : ",a,"\n")
    cat("true Kendall's parameter (tau) : ",Copula.tau,"\n")
    
    #--- Saving the esults in a text file
    sink("Results.txt",append=TRUE)
    cat("============================================","\n")
    cat("Results of                      : ",Copula.Name, "\n")
    cat("degree of Bernstein polynomial  : ",q, "\n")
    cat("b1 estimates                    : ",MLE.b1, "\n")
    cat("b2 estimates                    : ",MLE.b2, "\n")
    cat("phi estimates                   : ",MLE.phi, "\n")
    cat("BIC                             : ",BIC,"\n")
    cat("interval of generating c1       : ","(",l1,",",u1,")\n")  
    cat("interval of generating c2       : ","(",l2,",",u2,")\n")  
    cat("sample size (n)                 : ",n,"\n")
    cat("no. of simulations              : ",Sim.No,"\n")
    cat("true association parameter (a)  : ",a,"\n")
    cat("true Kendall's parameter (tau)  : ",Copula.tau,"\n")
    cat("Start time                      : ",Start.date,"\n")
    cat("End time                        : ",End.date,"\n")
    cat("--------------------------------------------","\n")
    print(RES)
    sink()
    
    ####---- Save as Text files
    write.table(round(Beta.Est,3),file=paste("beta_n",n,"_q",q,"_tau",round(Copula.tau*100),"_",Copula.Name,".txt",sep=""))
    write.table(round(ESD.Prof,3),file=paste("esd_prof_n",n,"_q",q,"_tau",round(Copula.tau*100),"_",Copula.Name,".txt",sep=""))
    write.table(round(CP.Prof,3),file=paste("cp_prof_n",n,"_q",q,"_tau",round(Copula.tau*100),"_",Copula.Name,".txt",sep=""))
    write.table(round(b1.Est,3),file  =paste("b1_n",n,"_q",q,"_tau",round(Copula.tau*100),"_",Copula.Name,".txt",sep=""))
    write.table(round(b2.Est,3),file  =paste("b2_n",n,"_q",q,"_tau",round(Copula.tau*100),"_",Copula.Name,".txt",sep=""))
    write.table(round(phi.Est,3),file =paste("phi_n",n,"_q",q,"_tau",round(Copula.tau*100),"_",Copula.Name,".txt",sep=""))
    #---- Ploting the cumulative hazard functions of T1 and T2
    
    dev.new()
    postscript(paste("Graph_n",n,"_q",q,"_tau",round(Copula.tau*100),"_",Copula.Name,".eps",sep=""), paper="letter", onefile=FALSE, horizontal=FALSE)
    
    par(mfrow=c(1,2))
    
    Lambda1.fit=BB1%*%apply(b1.Est,2,median)    # fitted cumulative hazards function
    matplot(c.axis.1, cbind(Lambda1.fit,Lambda1.true), type = "l",main=expression(paste("Fitted ",hat(H)[1] (c))),
            xlab = "c", ylab=expression(paste("baseline cumulative hazards function, ", H[1])),
            lwd = 2, lty =c(1,2), col = c("black","red"),ylim=c(0,1.2))
    legend("topleft", c(expression(paste("fitted ",hat(H)[1])),expression(paste("true ",hat(H)[1]))),lty =c(1,2), lwd = 1, col =c("black","red"), inset = 0.00001)
    
    Lambda2.fit=BB2%*%apply(b2.Est,2,median)    # fitted cumulative hazards function
    
    matplot(c.axis.2, cbind(Lambda2.fit,Lambda2.true), type = "l",main=expression(paste("Fitted ",hat(H)[2] (c))),
            xlab = "c", ylab=expression(paste("baseline cumulative hazards function, ", H[2])),
            lwd =2, lty = c(1,2), col = c("black","red"),xlim=c(l2,u2),ylim=c(0,1.2))
    legend("topleft",c(expression(paste("fitted ",hat(H)[2])),expression(paste("true ",hat(H)[2]))),lty =c(1,2), lwd = 1, col =c("black","red"), inset = 0.00001)
    
    dev.off()
    
  }
}  # end of the main loop


