##### generate data ####
set.seed(20)
###parameter setting##
n<-200
simu<-100
###parameter for cure rate##
intercept<-1 
a1<-1
a2<--1
#### parameter for survival model###
b1<--0.5
b2<--1
#### examination times ###
p<-10   # decide the right point
v<-10  # decide the left point
order<-3  #decide the order of i spline
MaxIter=200
generate<-function(n,intercept,a1,a2,b1,b2,p,v)
{
  #n sample size
  #intercept cofficients of cure rate
  #a1,a2     cofficients of cure rate
  #b1,b2     cofficients of PH model
  #p         examination times
  #v         first generate uniform range
  #status    censording indicator, 0-left censored,1-interval censored, 2-right censored
  #L         left point of all objects
  #R         right point of all objects
  #Y         possible examination time
  #cure rate indicator
  z1<-runif(n,0,2)  #or normal 
  z2<-rbinom(n,1,0.5)  
  linpred<-cbind(1,z1,z2) %*% c(intercept,a1,a2)
  prob<-exp(linpred)/(1+exp(linpred))
  y<-rbinom(n=n, size=1, prob=prob)
  #survival probabilities exponentional(1) and time calculation
  u<-runif(n,0,1)
  x<-cbind(z1,z2)
  time<-qweibull(1-exp(log(u)*exp(-x%*%c(b1,b2))),shape=6,scale=8)
  #time<-rexp(n,exp(x%*%b1))
  ###generate examination times####
  C<-matrix(rep(0,n),ncol = 1)   # all right point of subject
  status<-matrix(rep(0,n),ncol=1)
  L<-matrix(rep(0,n),ncol=1)
  R<-matrix(rep(0,n),ncol=1)
  for(j in 1:n)
  { 
    P<-rpois(1,p)+1   # examination times 
    Y<-matrix(rep(0,P),ncol =1)
    Y[1]<-runif(1,0,v)
    len<-runif(1,1,4)  # decide the length of interval 
    for (i in 2:P) 
    {
      Y[i]<-Y[1]+len*(i-1)  
    }               
    C[j]<-Y[P]
    T<-time*y+C*(1-y) # generate exact failure data
    if(y[j]==0) #cure fraction
    {
      status[j]<-2
      L[j]<-Y[P]
      R[j]<-NA    #(L_i,R_i)=(Y_p, infinite)
    } 
    if(y[j]==1) #uncure fraction
    {
      if(T[j]<Y[1]) 
      {
        status[j]<-0
        L[j]<-0
        R[j]<-Y[1]
      } #left censored, (L_j,R_j)=(0, Y_1)
      
      if(T[j]>Y[P])
      {
        status[j]<-2
        L[j]<-Y[P]
        R[j]<-NA
      }  #right censored,(L_j,R_j)=(Y_p,NA)
      if(T[j]<Y[P]&T[j]>Y[1])
      { status[j]<-1 
      i<-1
      repeat
      {
        i<-i+1
        if(T[j]<Y[i])
        {break}
      }
      L[j]<-Y[i-1]
      R[j]<-Y[i]  # find the interval contain T,(L_j,R_j)=(Y_i-1,Y_i)
      }
    }
  }
  data<-data.frame(T,L,R,status,z1,z2)
  return(data)
}

#library(HI)
Ispline<-function(x,order,knots)
{
  # M Spline function with order k=order+1. or I spline with order
  # x is a row vector
  # k is the order of I spline
  # knots are a sequence of increasing points
  # the number of free parameters in M spline is the length of knots plus 1.
  
  k=order+1
  m=length(knots)
  n=m-2+k # number of parameters
  t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots
  
  yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
  for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
  }
  
  yytem1=yy1
  for (ii in 1:order){
    yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
    for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
    }
    yytem1=yytem2
  }
  
  index=rep(0,length(x))
  for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
  }
  
  yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))
  
  if (order==1){
    for (i in 2:n){
      yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
    }
  }else{
    for (j in 1:length(x)){
      for (i in 2:n){
        if (i<(index[j]-order+1)){
          yy[i-1,j]=1
        }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
          yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
        }else{
          yy[i-1,j]=0
        }
      }
    }
  }
  return(yy)
}
positivepoissonrnd<-function(lambda){ 
  samp = rpois(1, lambda)
  while (samp==0) {
    samp = rpois(1, lambda)
  }
  return(samp)
}
lambda <- function(x) 
{
  nzi <- (1:length(x))[x!=0]   ## avoid the 0 in the x, if it is 0,return 0.125
  ans <- rep(0.125,length(x))
  ans[nzi] <- -tanh(x[nzi]/2)/(4*x[nzi]) ### negative?
  return(ans)
}
phi <- function(x) 
{
  nzi <- (1:length(x))[x!=0]
  ans <- rep(0.125,length(x))
  ans[nzi] <- (x[nzi]/2) - log(1+exp(x[nzi])) + (0.25*x[nzi]*tanh(x[nzi]/2))
  return(ans)
} 

####objective function h for optimization beta step
h.beta<-function(beta,X,u,gamma,bisu,bisv,V,W,status,Hyper)
{
  D=(status==1)+(status==2)
  first=sum(u*(bisu%*%gamma)*exp(X%*%beta))+sum(u*D*((bisv-bisu)%*%gamma)*exp(X%*%beta))
  second=sum(u*V*(X%*%beta))+sum(u*D*W*(X%*%beta))
  third=crossprod(beta)*(1/Hyper$sigma.beta)/2
  h.beta=first-second+third
  return(h.beta)
}


## Jacobian of h.beta
grad.h.beta<-function(beta,X,status,u,bisu,bisv,gamma,V,W,Hyper)
{
  D=(status==1)+(status==2)
  First=crossprod(X,u*(bisu%*%gamma)*exp(X%*%beta))+crossprod(X,u*D*((bisv-bisu)%*%gamma)*exp(X%*%beta))
  Second=crossprod(X,u*V)+crossprod(X,u*D*W)
  Third=beta*(1/Hyper$sigma.beta)
  grad.h.beta=First-Second+Third
  return(grad.h.beta)
}

h.extr<-function(extr,X,beta,mu.q.u,bisu,bisv,status,V,W,Hyper)
{
  D=(status==1)+(status==2)
  mu.q.gamma=exp(extr)
  C=exp(X%*%beta+rowSums(X^2*Hyper$sigma.beta/2))
  modi=bisu%*%mu.q.gamma+10^-5
  modi1=D*((bisv-bisu)%*%mu.q.gamma)+10^-5
  First=sum(mu.q.u*C*(bisu%*%mu.q.gamma))+sum(mu.q.u*C*D*((bisv-bisu)%*%mu.q.gamma))
  Second=sum(mu.q.u*V*log(modi))+sum(mu.q.u*W*D*log(modi1))
  Third=(1/Hyper$gamma.shape-1)*sum(log(mu.q.gamma))-1/Hyper$gamma.shape*sum(mu.q.gamma)
  h.gamma=First-Second-Third
  return(h.gamma)
}

grad.h.extr<-function(extr,X,beta,mu.q.u,bisu,bisv,status,V,W,Hyper)
{
  D=(status==1)+(status==2)
  mu.q.gamma=exp(extr)
  C=exp(X%*%beta+rowSums(X^2*Hyper$sigma.beta/2))
  modi=bisu%*%mu.q.gamma+10^-5
  modi1=D*((bisv-bisu)%*%mu.q.gamma)+10^-5
  First=crossprod(bisu,mu.q.u*C)+crossprod(bisv-bisu,mu.q.u*C*D)
  Second=crossprod(bisu,modi*mu.q.u*V)+crossprod(diag(as.numeric(D))%*%(bisv-bisu),modi1*mu.q.u*W*D)
  prior= (1/Hyper$gamma.shape-1)*(mu.q.gamma^-1)-1/Hyper$gamma.shape
  grad.h.extr=(First-Second-prior)*exp(extr)
  return(grad.h.extr)
}

h.u<-function(X,mu.q.alpha,mu.q.beta,bisu,bisv,status,Hyper,mu.q.gamma,V,W,Z)
{
  D=(status==1)+(status==2)
  C=exp(X%*%mu.q.beta+rowSums(X^2*Hyper$sigma.beta/2))
  c_1=-(bisu%*%mu.q.gamma)*C+V*(bisu%*%mu.q.gamma-1+X%*%mu.q.beta)-log(factorial(V))
  c_2=-D*(bisv-bisu)%*%mu.q.gamma*C+W*((bisv-bisu)%*%mu.q.gamma-1+X%*%mu.q.beta)-log(factorial(W))
  mu=2*Z%*%mu.q.alpha+c_1+c_2
  r=(1+exp(-mu))^-1
  return(r)
}

data<-generate(n,intercept,a1,a2,b1,b2,p,v)
tol<-0.01

### main routine ###
L=matrix(data$L,ncol=1)
R=matrix(data$R,ncol=1)
R2=ifelse(is.na(R),0,R)
xcov=matrix(c(data$z1,data$z2),ncol=2)
p=ncol(xcov)
status=matrix(data$status,ncol=1)
n=nrow(L)
u<-L*as.numeric(status==1)+R2*as.numeric(status==0)     #t_i1        
v<-L*as.numeric(status==2)+R2*as.numeric(status==1)     #t_i2
obs=cbind(u,v)


## generate basis functions
knots<-seq(min(obs),max(obs),length=10)
k=length(knots)-2+order
#grids<-seq(min(obs),max(obs),length=100)
#kgrids=length(grids)
#G<-length(x_user)/p     # number of survival curves

#bibl=Ispline(t(bl.Li),order,knots)
bisu=Ispline(t(obs[,1]),order,knots)
bisv=Ispline(t(obs[,2]),order,knots)
#bgs=Ispline(grids,order,knots)
sum((bisv-bisu)<0)

X<-cbind(data$z1,data$z2)
Z<-cbind(1,X)  
bisu<-t(bisu)  # transport to a N*k matrix
bisv<-t(bisv)

#MFVB<-function(X,Z,bisu,bisv,MaxIter=200)
#{
## Argument 
##         T   Failture time
##         R   right point of censoring time
##         L   left point of censoring time
##    status   censoring indicator
##         X   matrix of covariates 
##   MaxIter   Maximum number of iterations for the variational approximation steps  
##      bisu   Ispline function of t.ij1, matrix with dimension N*K
##      bisv   Ispline fucntion of t.ij2, matrix wiht dimension N*K
##         N   total number of observision, m*n.i
##         V   the count of occurence of subject j in cluster i up to time t_ij1
##         W   the count of occurence of subject j in cluster i up to time t_ij2

#### Define Hyper Structure
Hyper<-list(alpha0=rep(0,ncol(Z)),    ## alpha0 should be a vector?
            sigma.alpha=1,         ## the choose of Hyper Structure
            beta0=rep(0,ncol(X)),
            sigma.beta=1,
            gamma.shape=10)         ## the prior of gamma?

####Define dimension structure
Dim<-list(N=nrow(Z),
          P=ncol(X),
          K=ncol(bisu)   ## bisu: a matrix with dim K*N 
          
)

####Build parameter tracks
MaxIter=200
mu.q.alpha.track<-matrix(NA,nrow = MaxIter,ncol = Dim$P+1)
mu.q.beta.track<-matrix(NA,nrow = MaxIter,ncol = Dim$P)
extr.track<-matrix(NA,nrow=MaxIter,ncol=Dim$K)
mu.q.gamma.track<-matrix(NA,nrow=MaxIter,ncol=Dim$K)
mu.q.u.track<-matrix(NA,nrow = MaxIter,ncol = Dim$N)   ##????
Sigma.q.alpha.track<-vector("list",length = MaxIter)
Sigma.q.beta.track<-vector("list",length= MaxIter)
Sigma.q.gamma.track<-vector("list",length= MaxIter)


###convergence trackers
log.p.track = rep(NA,MaxIter)
StepConv = rep(NA,MaxIter)
StepConvBound = 0.01 #relconv 1% relative convergence for now


### intializing parameters  ?????? the choice of intial value
cat("Intializing parameters \n")    
Xi=rep(1,Dim$N)                   # N dimension 
mu.q.u=rep(1,Dim$N)
mu.q.beta=c(-0.5,-1)
extr=rep(0,Dim$K)
mu.q.gamma=exp(extr)  ### how to choose a suit gamma intial value?

##
Lambdatu=bisu%*%mu.q.gamma# n x 1 t_i1
Lambdatv=bisv%*%mu.q.gamma



###Initialize first values of trackers
mu.q.alpha.track[1,]=rep(1,Dim$P+1)
mu.q.beta.track[1,]=mu.q.beta
extr.track[1,]=extr
mu.q.gamma.track[1,]=mu.q.gamma
mu.q.u.track[1,]=mu.q.u
Sigma.q.alpha.track[[1]]=diag(rep(1,Dim$P+1))
Sigma.q.beta.track[[1]]=diag(rep(1,Dim$P))
Sigma.q.gamma.track[[1]]=diag(rep(1,Dim$K))

###Define convergence objects 
KeepGoing=T
conv = rep(NA,MaxIter) #tracks convergence of optim




###Get into the variational loop
StartTime = proc.time()
i=1
while(KeepGoing){
  i=i+1 
  
  ###sample V and W from distribution
  N<-Dim$N
  temp_beta<-mu.q.beta.track[i-1,]
  V=array(rep(0,N),dim=c(N,1)); W=V
  for (j in 1:N){
    if (status[j]==0){
      templam1=Lambdatu[j]*exp(X[j,]%*%temp_beta)
      V[j]=positivepoissonrnd(templam1)
      #zz[i,]=rmultinom(1,z[i],gamcoef*t(bisu[,i]))
    }else if (status[j]==1){
      templam1=(Lambdatv[j]-Lambdatu[j])*exp(X[j,]%*%temp_beta)
      W[j]=positivepoissonrnd(templam1)
      #ww[i,]=rmultinom(1,w[i],gamcoef*t(bisv[,i]-bisu[,i]))
    }
  }
  
  for(j in 1:N){
    if(V[j]>100){
      V[j]<-100
    }
  }
  for(j in 1:N){
    if(W[j]>100){
      W[j]<-100
    }
  }
  ##avoid the overwhelming of factorial(z) and factorial(w)
  
  ###update mu.q.alpha and Sigma.q.alpha
  
  wtVec <- lambda(Xi)
  AVec <- phi(Xi) 
  ridgeVec <- rep(1/Hyper$sigma.alpha,Dim$P+1)
  M.q.Sigma<-4*crossprod(Z,wtVec*Z) 
  
  ### Estimate the Sigma.q.alpha
  Htmp<-diag(ridgeVec)-M.q.Sigma  
  
  ### check singularity
  solve_check = try(solve(Htmp),T)
  if(class(solve_check)=="matrix"){
    #solve worked fine
    Sigma.q.alpha = solve(Htmp)
    singular_warning = F
  }else{
    #Use Moore-Penrose pseudoinverse
    Sigma.q.alpha = ginv(Htmp)
    singular_warning = T
  }
  
  
  ### check positive defination ###
  eigen_val = eigen(Sigma.q.alpha)$values
  min_eigen_val = min(eigen_val)
  if(all(eigen_val>0)){
    #positive definite
    Sigma.q.alpha.track[[i]]<- Sigma.q.alpha
    pd_warning = F
  }else{
    #use a ridge constant to force positive definite
    ridge_constant = abs(min_eigen_val)
    Sigma.q.alpha.track[[i]]=Sigma.q.alpha +diag(ridge_constant,nrow = nrow(Sigma.q.alpha))
    pd_warning = T
  }
  if(singular_warning){
    cat("** Warning: Using ginv for Sigma.q.alpha**\n")
  }
  
  if(pd_warning){
    cat("** Warning: Ridge Shift for Sigma.q.alpha **\n")
  }
  ## Estimate the mu.q.alpha 
  
  mu.q.alpha.track[i,]<-Sigma.q.alpha.track[[i]]%*%crossprod(Z,2*mu.q.u.track[i-1,]-1)  ## the expection of u ?
  
  ## updata Xi
  
  EsqMat <- Sigma.q.alpha.track[[i]] + tcrossprod(mu.q.alpha.track[i,])
  Xi <- sqrt(diag(Z%*%EsqMat%*%t(Z)))
  
  
  
  ### mu.q.beta  and Sigma.q.beta via optimization of start value
  #optim.starttime<-proc.time()  #??? the choice of start value 
  foo<-optim(par=mu.q.beta.track[i-1,],fn=h.beta,gr=grad.h.beta,X=X,u=mu.q.u.track[i-1,],gamma=mu.q.gamma.track[i-1,],
             bisu=bisu,bisv=bisv,V=V,W=W,status=status,Hyper=Hyper,
             lower = -Inf, upper = Inf,
             control = list(MaxIter), hessian = T)
  #optim.time=proc.time()-optim.starttime
  
  # Estimate mu.q.beta 
  mu.q.beta.track[i,]=foo$par
  
  # Estimate Sigma.q.beta
  Htmp.beta=foo$hessian
  
  #check Singularity 
  solve_check1=try(solve(Htmp.beta),"matrix")
  if(class(solve_check1)=="matrix")
  {
    #solve work well
    Sigma.q.beta=solve(Htmp.beta)
    singular_warning.beta=FALSE
  }else{
    Sigma.q.beta=ginv(Htmp.beta)
    singular_warinig.beta=TRUE
  }
  
  # Check positive
  eign_val1=eigen(Sigma.q.beta)$values
  min_eigen_val1=min(eign_val1)
  if(all(min_eigen_val1>0))
  {
    Sigma.q.beta.track[[i]]=Sigma.q.beta
    pd_warning1=FALSE
  }else{
    ridge_constant1=abs(min_eigen_val1)
    Sigma.q.beta.track[[i]]=Sigma.q.beta +diag(ridge_constant1,nrow = nrow(Sigma.q.beta))
    pd_warning1=TRUE
  }
  
  if(singular_warning.beta){
    cat("** Warning: Using ginv for Sigma.q.beta**\n")
  }
  
  if( pd_warning1){
    cat("** Warning: Ridge Shift for Sigma.q.beta **\n")
  }
  
  
  ## Estimate extr  and mu.q.gamma
   foo1<-optim(par=extr.track[i-1,],fn=h.extr,gr=grad.h.extr,X=X,beta=mu.q.beta.track[i,],mu.q.u=mu.q.u.track[i-1,],
   bisu=bisu,bisv=bisv,status=status,V=V,W=W,Hyper=Hyper,hessian=TRUE)   
  
    extr.track[i,]=foo1$par
    Htmp.extr=foo1$hessian
  
  
   mu.q.gamma.track[i,]<-exp(extr.track[i,])
   inv=diag(1/exp(extr.track[i,]))
   fir_dif=grad.h.extr(extr.track[i,],X=X,beta=mu.q.beta.track[i,],mu.q.u=mu.q.u.track[i-1,],
   bisu=bisu,bisv=bisv,status=status,V=V,W=W,Hyper=Hyper)
   Htmp.gamma=inv%*%(Htmp.extr-diag(as.vector(fir_dif*exp(extr.track[i,]))))%*%inv  ###???? first difference?

   
   #check singularity
   solve_check2=try(solve(Htmp.gamma),TRUE)
   if(class(solve_check2)=="matrix")
   {
     Sigma.q.gamma=solve(Htmp.gamma)
     singular_warning.gamma=FALSE
   }else{
     #Use Moore-Penrose pseudoinverse
     Sigma.q.gamma=ginv(Htmp.gamma)
     singular_warinig.gamma=TRUE
   }
   
   
   #Check positive 
   eign_val2=eigen(Sigma.q.gamma)$value
   min_eigen_val2=min(eign_val2)
   if(all(min_eigen_val2>0))
   {
     Sigma.q.gamma.track[[i]]=Sigma.q.gamma
     pd_warning2=FALSE
   }else{
     ridge_constant2=abs(min_eigen_val2)
     Sigma.q.gamma.track[[i]]=Sigma.q.gamma+diag(ridge_constant2,nrow = nrow(Sigma.q.gamma))
     pd_warning2=TRUE
   }#?? the dimension of Sigma.q.gamma
   
   if(singular_warning.gamma){
     cat("** Warning: Using ginv for Sigma.q.beta**\n")
   }
   
   if(pd_warning2){
     cat("** Warning: Ridge Shift for Sigma.q.beta **\n")
   }
   
   #Updata mu.q.u 
   
   mu.q.u.track[i,]<-h.u(X=X,mu.q.alpha=mu.q.alpha.track[i,],mu.q.beta=mu.q.beta.track[i,],bisu=bisu,bisv=bisv,status=status,Hyper=Hyper,mu.q.gamma=mu.q.gamma.track[i,],
                         V=V,W=W,Z=Z)


