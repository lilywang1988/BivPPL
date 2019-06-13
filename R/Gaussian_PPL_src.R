# library(mvtnorm) 	#>=version 1.0-6
# library(survival) #>=	version 2.42-3
# library(RcppArmadillo) #>=version 0.7.960.1.2	
# library(Rcpp) # >=version 0.12.13
# library(Matrix) #version >=1.2-11

BivPPL<-function(data,max.iter=1000,eps=1e-5,NP.max.iter=500,NP.eps=1e-6,huge=FALSE,independence=FALSE,D.initial=0.5,raw=TRUE,check.log.lik=FALSE){  
  #initializing regression parameters beta using coxph() from R package survival
  #initializing the D matrix
  
  if(raw){  # raw data generated from the gen.data ()
    N<-max(length(data$x),length(data$y))
    z1<-t(do.call(cbind,data$z1))
    z2<-t(do.call(cbind,data$z2))
    mi<-sapply(data$x,length)
    beta1_hat<-coxph(Surv(unlist(data$x),unlist(data$delta1))~z1)$coefficients
    beta2_hat<-coxph(Surv(unlist(data$y),unlist(data$delta2))~z2)$coefficients
    q1<-length(beta1_hat)
    q2<-length(beta2_hat)
    r_hat<-rep(0,2*N)
    D_hat<-diag(2)*D.initial
    ppl_est<-c(beta1_hat,beta2_hat,r_hat)
    q<-q1+q2
    data_new<-list(x=unlist(data$x),y=unlist(data$y),delta1=unlist(data$delta1),delta2=unlist(data$delta2),z1=z1,z2=z2,N=N,mi=mi)
  }else{
    if('z1' %in% names(data)){ # list of variables
      N=data$N
      z1<-as.matrix(data$z1)
      z2<-as.matrix(data$z2)
      beta1_hat<-coxph(Surv(data$x,data$delta1)~z1)$coefficients
      beta2_hat<-coxph(Surv(data$y,data$delta2)~z2)$coefficients
      q1<-length(beta1_hat)
      q2<-length(beta2_hat)
      r_hat<-rep(0,2*N)
      D_hat<-diag(2)*D.initial
      ppl_est<-c(beta1_hat,beta2_hat,r_hat)
      q<-q1+q2
      data_new<-data
    }else if('Z1' %in% names(data)){ # transformed by tocoxme() or tocoxme
      nobs=nrow(data$Z1)/2
      N=data$N
      z1<-as.matrix(data$Z1)
      z2<-as.matrix(data$Z2)
      beta1_hat<-coxph(Surv(data$time,data$delta)~z1)$coefficients
      beta2_hat<-coxph(Surv(data$time,data$delta)~z2)$coefficients
      q1<-length(beta1_hat)
      q2<-length(beta2_hat)
      r_hat<-rep(0,2*N)
      D_hat<-diag(2)*D.initial
      ppl_est<-c(beta1_hat,beta2_hat,r_hat)
      q<-q1+q2
      data_new<-list(z1=z1[1:nobs,],z2=z2[(nobs+1):(2*nobs),],x=data$time[1:nobs],y=data$time[(nobs+1):(2*nobs)],delta1=data$delta[1:nobs],delta2=data$delta[(nobs+1):(2*nobs)],N=N,mi=data$mi)
    } else{
      stop("Invalid definition on \'Raw\' or data input ")
    }
    
  }
  d=100 #an arbitrary but reasonably large value
  #call the PPL.est() for inner loop, iterate for outer loop 
  for(t in 1:max.iter){
    ppl_est.p<- ppl_est
    D_hat.p<- D_hat
    #####################
   {
     if(t==1) ppl_res<-PPL.est(data_new,ppl_est.p,D_hat.p,q1,q2,100,NP.eps*100,huge,check.log.lik)
     else ppl_res<-PPL.est(data_new,ppl_est.p,D_hat,q1,q2,NP.max.iter,NP.eps,huge,check.log.lik)
   }  
    ppl_est<-ppl_res$coefficients
    r_hat<-ppl_est[-(1:q)] #the predicted random intercepts
    r1_hat<-r_hat[(1:N)*2-1]
    r2_hat<-r_hat[(1:N)*2]
    r_cbind<-cbind(r1_hat,r2_hat)
    Kppl.inv<-ppl_res$Kppl.inv
    ASE<-ppl_res$ASE
    #huge=T: sparse the Kppl matrix before its inverse
    if(huge==F)D_hat<-(t(r_cbind)%*%r_cbind+Reduce('+',lapply(1:N,function(i) Kppl.inv[(2*i-1):(2*i),(2*i-1):(2*i)])))/N else 
      D_hat<-(t(r_cbind)%*%r_cbind/N+apply(Kppl.inv,1:2,mean))
    if(independence==T) D_hat<-diag(diag(D_hat)) #for the ones we want to assume indepencent two process, i.e. equivalent to fit two sepeately using coxme
    ####################
    {
      if(!check.log.lik)d=max(abs(ppl_est[1:q]-ppl_est.p[1:q]),max(abs(D_hat-D_hat.p))) #@only consider the convergence of interesting parameters
      else {
        if(huge==F) Kppl.inv.det.log<-Reduce('+',lapply(1:N,function(i) log(det(Kppl.inv[(2*(i-1)+1):(2*i),(2*(i-1)+1):(2*i)])))) else
        Kppl.inv.det.log<-sum(sapply(1:N,function(i) log(det(Kppl.inv[,,i]))))
        LPm1<--N/2*log(det(D_hat))+Kppl.inv.det.log/2+ppl_res$PPL1
        if(t>1) d=abs(LPm1.p-LPm1)
        LPm1.p=LPm1
        }
    }
    if(d<eps |t>=max.iter) break
  }
  #Preparation for the likelihood ratio test 
  if(!check.log.lik){
    if(huge==F) Kppl.inv.det.log<-Reduce('+',lapply(1:N,function(i) log(det(Kppl.inv[(2*(i-1)+1):(2*i),(2*(i-1)+1):(2*i)])))) else
      Kppl.inv.det.log<-sum(sapply(1:N,function(i) log(det(Kppl.inv[,,i]))))
     LPm1<--N/2*log(det(D_hat))+Kppl.inv.det.log/2+ppl_res$PPL1
  }
  return(list(beta_hat=ppl_est[1:q],beta_ASE=ASE,D_hat=D_hat,r_hat=r_hat,LogMargProb=LPm1,num.iter=t,dist=d,   
              Kppl.inv.det.log=Kppl.inv.det.log,PPL1=ppl_res$PPL1,PPL2=ppl_res$PPL2,Kppl.inv=ppl_res$Kppl.inv,Kppl=ppl_res$Kppl)) 
}
#The function of the Newton-Raphson for the inner loop
#NP.max.iter restrict the number of Newton-Raphson interations.
PPL.est<-function(data_new,ppl_est,D_hat,q1,q2,NP.max.iter=1000,NP.eps=1e-6,huge=F,check.log.lik=F){
  d=100
  for(t in 1:NP.max.iter){
    ppl_est.p<-ppl_est
    ppl.prep<-PPL.prep(data_new,ppl_est.p,D_hat,q1,q2,huge)
    ppl_est<-ppl_est.p+ppl.prep$gradient
    if(check.log.lik && t>1){
      d=abs(ppl.prep$PPL1-PPL1.p)
    }else if (!check.log.lik) d=max(abs(ppl.prep$gradient))
    PPL1.p<-ppl.prep$PPL1
    if(d<NP.eps) break
  }
  return(list(NP.iter=t,coefficients=ppl_est,ASE=ppl.prep$ASE,Kppl.inv=ppl.prep$Kppl.inv,Kppl=ppl.prep$Kppl,PPL1=ppl.prep$PPL1,PPL2=ppl.prep$PPL2))
}
# Each estimating step of the inner loop
PPL.prep<-function(data_new,ppl_est,D_hat,q1,q2,huge=F){
  if(huge==F) { 
    q=q1+q2
    N=data_new$N
    mi=data_new$mi
    beta1_hat<-ppl_est[1:q1]
    beta2_hat<-ppl_est[(q1+1):(q)]
    r_hat<-ppl_est[-(1:(q))]
    r1_hat<-r_hat[(1:N)*2-1]
    r2_hat<-r_hat[(1:N)*2]
    R<-diag(2*N)
    R_select<-rep(1:N,times=mi)
    r1_hat_<-r1_hat[R_select]
    r2_hat_<-r2_hat[R_select]
    R1<-R[2*R_select-1,]
    R2<-R[2*R_select,]
    beta2_scr<-rep(0,q2)
    beta1_scr<-rep(0,q1)
    Sigma_inv<-as.matrix(kronecker(diag(N),solve(D_hat)))
    r_inf=matrix(Sigma_inv,2*N,2*N) #@
    r_scr<-c(-Sigma_inv%*%r_hat)
    beta2_r_inf<-matrix(0,q2,2*N)
    beta1_r_inf<-matrix(0,q1,2*N)
    beta2_inf<-matrix(0,q2,q2)
    beta1_inf<-matrix(0,q1,q1)
    select1<-order(data_new$x)
    select2<-order(data_new$y)
    z1<-as.matrix(as.matrix(data_new$z1)[select1,])
    x<-data_new$x[select1]
    delta1<-data_new$delta1[select1]
    r1_hat_pool<-r1_hat_[select1]
    z2<-as.matrix(as.matrix(data_new$z2)[select2,])
    y<-data_new$y[select2]
    delta2<-data_new$delta2[select2]
    r2_hat_pool<-r2_hat_[select2]
    R1_pool<-R1[select1,]
    R2_pool<-R2[select2,]
    eta1<-z1%*%as.matrix(unlist(beta1_hat),q1,1)+r1_hat_pool
    eta2<-z2%*%as.matrix(unlist(beta2_hat),q2,1)+r2_hat_pool
    exp_eta1<-exp(eta1)
    exp_eta2<-exp(eta2)
    m=sum(mi) 
    beta_PPL_loop(R1_pool,R2_pool,z1,z2,delta1,delta2,exp_eta1,exp_eta2,beta1_scr,beta2_scr,r_scr,
                  beta1_inf,beta2_inf,r_inf,beta1_r_inf,beta2_r_inf,m,eta1,eta2)
    PPL2=-t(r_hat)%*%Sigma_inv%*%r_hat/2
    PPL1=PPL2+sum(delta1*(eta1-log(rev(cumsum(rev(exp_eta1))))))+sum(delta2*(eta2-log(rev(cumsum(rev(exp_eta2)))))) #@
    #print(-t(r_hat)%*%r_inf%*%r_hat/2)
    #get the information matrix
    inf.inv<-inf<-matrix(0,2*N+q,2*N+q)
    inf[1:q1,1:q1]<-beta1_inf
    inf[(q1+1):(q),(q1+1):(q)]<-beta2_inf
    inf[(q+1):(q+2*N),(q+1):(q+2*N)]<-r_inf
    inf[1:q1,(q+1):(q+2*N)]<-beta1_r_inf
    inf[(q1+1):(q),(q+1):(q+2*N)]<-beta2_r_inf
    inf[(q+1):(q+2*N),1:q1]<-t(beta1_r_inf)
    inf[(q+1):(q+2*N),(q1+1):(q)]<-t(beta2_r_inf)  
    #Kppl.inv<-qr.solve(r_inf) #qr.inf is believed to be faster and more stable than solve
    Kppl.inv<-chol2inv(chol(r_inf)) #Choleski decomposed inverse is believed to be even faster than solve 
    #get the inverse information matrix
    inf.inv[1:q,1:q]<-solve(inf[1:q,1:q]-inf[1:q,(q+1):(q+2*N)]%*%Kppl.inv%*%inf[(q+1):(q+2*N),1:q])
    inf.inv[(q+1):(q+2*N),(q+1):(q+2*N)]<-solve(r_inf-inf[(q+1):(q+2*N),1:q]%*%solve(inf[1:q,1:q])%*%inf[1:q,(q+1):(q+2*N)])
    inf.inv[1:q,(q+1):(q+2*N)]<--inf.inv[1:q,1:q]%*%inf[1:q,(q+1):(q+2*N)]%*%Kppl.inv
    inf.inv[(q+1):(q+2*N),1:q]<--inf.inv[(q+1):(q+2*N),(q+1):(q+2*N)]%*%inf[(q+1):(q+2*N),1:q]%*%solve(inf[1:q,1:q])
    return(list(gradient=inf.inv%*%c(beta1_scr,beta2_scr,r_scr),ASE=sqrt(diag(inf.inv)[1:q]),Kppl.inv=Kppl.inv,Kppl=r_inf,PPL1=PPL1,PPL2=PPL2))
  } else {
    q=q1+q2
    N=data_new$N
    mi=data_new$mi
    beta1_hat<-ppl_est[1:q1]
    beta2_hat<-ppl_est[(q1+1):(q)]
    r_hat<-ppl_est[-(1:(q))]
    r1_hat<-r_hat[(1:N)*2-1]
    r2_hat<-r_hat[(1:N)*2]
    R<-rep(1:N,each=2)
    R_select<-rep(1:N,times=mi)
    r1_hat_<-r1_hat[R_select]
    r2_hat_<-r2_hat[R_select]
    R1<-R[2*R_select-1]
    R2<-R[2*R_select]
    beta2_scr<-rep(0,q2)
    beta1_scr<-rep(0,q1)
    Sigma_inv<-replicate(N,solve(D_hat)) #dim= 2 2 N
    r_inf<-array(Sigma_inv,dim=c(2,2,N))#@
    r_scr<-c(-sapply(1:N,function(i) Sigma_inv[,,i]%*%r_hat[(2*i-1):(2*i)])) #already checked
    beta2_r_inf<-matrix(0,q2,2*N)
    beta1_r_inf<-matrix(0,q1,2*N)
    beta2_inf<-matrix(0,q2,q2)
    beta1_inf<-matrix(0,q1,q1)
    select1<-order(data_new$x)
    select2<-order(data_new$y)
    z1<-as.matrix(as.matrix(data_new$z1)[select1,])
    x<-data_new$x[select1]
    delta1<-data_new$delta1[select1]
    r1_hat_pool<-r1_hat_[select1]
    z2<-as.matrix(as.matrix(data_new$z2)[select2,])
    y<-data_new$y[select2]
    delta2<-data_new$delta2[select2]
    r2_hat_pool<-r2_hat_[select2]
    R1_pool<-R1[select1]
    R2_pool<-R2[select2]
    eta1<-z1%*%as.matrix(unlist(beta1_hat),q1,1)+r1_hat_pool
    eta2<-z2%*%as.matrix(unlist(beta2_hat),q2,1)+r2_hat_pool
    exp_eta1<-exp(eta1)
    exp_eta2<-exp(eta2)
    m=sum(mi) 
    beta_PPL_loop_huge(R1_pool,R2_pool,z1,z2,delta1,delta2,exp_eta1,exp_eta2,beta1_scr,beta2_scr,r_scr,
                       beta1_inf,beta2_inf,r_inf,beta1_r_inf,beta2_r_inf,m,N,eta1,eta2) 
    PPL2=-sum(sapply(1:N,function(i) t(r_hat[(2*i-1):(2*i)])%*%solve(D_hat)%*%r_hat[(2*i-1):(2*i)]))/2
    PPL1=PPL2+sum(delta1*(eta1-log(rev(cumsum(rev(exp_eta1))))))+sum(delta2*(eta2-log(rev(cumsum(rev(exp_eta2)))))) #@
    r_inf.inv<-array(apply(r_inf,3,solve),dim=c(2,2,N)) #I will ignore the covariance between beta and gamma to avoid calculating and saving a large matrix
    beta_inf<-bdiag(beta1_inf,beta2_inf)
    beta_r_inf<-rbind(beta1_r_inf,beta2_r_inf)
    # we use sepearte graidents for beta and gamma, instead of one from the whole Hessian matrix
    beta_gradient<-chol2inv(chol(beta_inf))%*%c(beta1_scr,beta2_scr) 
    r_gradient<-c(sapply(1:N,function(i) chol2inv(chol(r_inf[,,i]))%*%r_scr[(2*i-1):(2*i)]))
    ASE<-sqrt(diag(solve(beta_inf- Reduce('+',lapply(1:N,function(i) beta_r_inf[,(2*i-1):(2*i)]%*%r_inf.inv[,,i]%*%t(beta_r_inf[,(2*i-1):(2*i)]))) )))
    return(list(gradient=c(as.vector(beta_gradient),r_gradient),ASE=ASE,Kppl.inv=r_inf.inv,Kppl=r_inf,PPL1=PPL1,PPL2=PPL2))
  }
}

#Generate data function gen.data(N,beta1,beta2,theta,c, Ctype,const_cov, same_cov) 
# N: sample size
# beta1: effects for event 1
# beta2: effects for event 2
# theta: a vector of entries of the D matrix (theta=[D[1,1],D[2,2],D[1,2]])
# c: the end of the follow-up time
# Ctype: the censoring type. If Ctype==1, fixed censoring time at c; if Ctype!=1, uniformly distributed censoring time following U[0,c]
# const_cov: whether fixed the covariates unchanged by event times.
# same_cov: whether the two events are dependent on the same batch of covariates or identical covariate values

gen.data<-function(N,beta1,beta2,theta,lambda01,lambda02,c=10,Ctype=TRUE,const_cov=FALSE,same_cov=FALSE){ 
  beta1=as.matrix(beta1)
  beta2=as.matrix(beta2)
  x<-vector("list",N)
  y<-vector("list",N)
  z1<-vector("list",N)
  z2<-vector("list",N)
  q1<-length(beta1)
  q2<-length(beta2)
  r<-matrix(NA,N,2)
  delta2<-delta1<-vector("list",N)
  D<-matrix(theta[c(1,3,3,2)],2,2)
  
  #Define censoring function: Ctype=1 administrative censoring at fixed times; Ctype!=1 uniform(0,c) censoring.
  C<-function(c,Ctype){
    if(Ctype==TRUE) return(c)
    else return(runif(1,0,c))
  }
  
  for(i in 1: N){
    j=0
    t1ij<-0
    t2ij<-0
    ri<-rmvnorm(n=1,mean=rep(0,2),sigma=D)
    r[i,]<-ri
    Ci<-C(c,Ctype)
    if(const_cov==T){
      if(same_cov==T&&q1==q2) z1ij<-z2ij<-t(rmvnorm(n=1,mean=rep(0,q2),sigma=diag(rep(1,q2)))) else{
        z1ij<-t(rmvnorm(n=1,mean=rep(0,q1),sigma=diag(rep(1,q1))))
        z2ij<-t(rmvnorm(n=1,mean=rep(0,q2),sigma=diag(rep(1,q2))))
      }
    }
    repeat{
      j=j+1
      #if(j==1) z1ij<-rnorm(1) else z1ij<-rnorm(1,z1ij,0.5)
      #if(j==1) z2ij<-rnorm(1,z1ij,1) else z2ij<-rnorm(1,z2ij,0.5)
      #print(const_cov)
      if(const_cov==F){
        if(same_cov==T&&q1==q2) z1ij<-z2ij<-t(rmvnorm(n=1,mean=rep(0,q2),sigma=diag(rep(1,q2)))) else{
          z1ij<-t(rmvnorm(n=1,mean=rep(0,q1),sigma=diag(rep(1,q1))))
          z2ij<-t(rmvnorm(n=1,mean=rep(0,q2),sigma=diag(rep(1,q2))))
        }
      }
      
      z1[[i]]<-cbind(z1[[i]],z1ij)
      z2[[i]]<-cbind(z2[[i]],z2ij)
      haz1ij<-lambda01*exp(t(beta1)%*%z1ij+ri[1])
      haz2ij<-lambda02*exp(t(beta2)%*%z2ij+ri[2])
      x_temp<--log(runif(1))/haz1ij
      y_temp<--log(runif(1))/haz2ij
      t1ij<-t2ij+x_temp
      if(t1ij>Ci) {delta1[[i]][j]=0;x[[i]][j]=Ci-t2ij;delta2[[i]][j]=0;y[[i]][j]=0;break} else {delta1[[i]][j]=1;x[[i]][j]=x_temp;}
      t2ij<-t1ij+y_temp
      if(t2ij>Ci) {delta2[[i]][j]=0;y[[i]][j]=Ci-t1ij;break} else {delta2[[i]][j]=1;y[[i]][j]=y_temp;}
    }
  }
  return(list(x=x,y=y,z1=z1,z2=z2,r=r,delta1=delta1,delta2=delta2))
}

#change the data to fit in coxme function in the coxme package, only for bivariate events 
tocoxme<-function(data){
  xtime<-ytime<-NULL
  ydelta<-xdelta<-NULL
  b0<-b1<-b2<-NULL
  Z1<-Z2<-NULL
  N<-length(data$x)
  mi_vec=NULL
  for(i in 1:N){
    mi<-length(data$y[[i]])
    xtime<-c(xtime,data$x[[i]])
    ytime<-c(ytime,data$y[[i]])
    xdelta<-c(xdelta,data$delta1[[i]])
    ydelta<-c(ydelta,data$delta2[[i]])
    if(i==1){
      Z1<-t(data$z1[[1]])
      Z2<-t(data$z2[[1]])
    }else{
      Z1<-rbind(Z1,t(data$z1[[i]]))
      Z2<-rbind(Z2,t(data$z2[[i]]))
    }
    
    b0<-c(b0,rep(i,times=mi))
    b1<-c(b1,rep(i,times=mi))
    b2<-c(b2,rep(i,times=mi))
    mi_vec[i]=mi
  }
  Z1_new<-rbind(Z1,matrix(0,nrow=nrow(Z2),ncol=ncol(Z1)))
  Z2_new<-rbind(matrix(0,nrow=nrow(Z1),ncol=ncol(Z2)),Z2)
  time<-c(xtime,ytime)
  delta<-c(xdelta,ydelta)
  joint<-c(rep(1,length(xdelta)),rep(2,length(ydelta)))
  b0_new<-rep(b0,times=2)
  b1_new<-c(b1,b2*0)
  b2_new<-c(b1*0,b2)
  return(list(time=time,delta=delta,joint=joint,b0=b0_new,b1=b1_new,b2=b2_new,Z1=Z1_new,Z2=Z2_new,mi=mi_vec,N=N))
}
#Examples: 
#result<-coxme(Surv(time,delta)~Z1_new+Z2_new+(1|b0_new)+(1|b1_new)+(1|b2_new)+strata(joint))
#result$coefficients
#as.vector(unlist(result$vcoef))
#sqrt(diag(vcov(result)))

#for negative correlations
tocoxme_n<-function(data){
  xtime<-ytime<-NULL
  ydelta<-xdelta<-NULL
  b0<-b1<-b2<-NULL
  Z01<-Z02<-Z1<-Z2<-NULL
  N<-length(data$x)
  mi_vec=NULL
  for(i in 1:N){
    mi<-length(data$y[[i]])
    xtime<-c(xtime,data$x[[i]])
    ytime<-c(ytime,data$y[[i]])
    xdelta<-c(xdelta,data$delta1[[i]])
    ydelta<-c(ydelta,data$delta2[[i]])
    if(i==1){
      Z1<-t(data$z1[[1]])
      Z2<-t(data$z2[[1]])
    }else{
      Z1<-rbind(Z1,t(data$z1[[i]]))
      Z2<-rbind(Z2,t(data$z2[[i]]))
    }
    
    b0<-c(b0,rep(i,times=mi))
    Z01<-c(Z01,rep(1,mi))
    Z02<-c(Z02,rep(-1,mi))
    b1<-c(b1,rep(i,times=mi))
    b2<-c(b2,rep(i,times=mi))
    mi_vec[i]=mi
  }
  Z1_new<-rbind(Z1,matrix(0,nrow=nrow(Z2),ncol=ncol(Z1)))
  Z2_new<-rbind(matrix(0,nrow=nrow(Z1),ncol=ncol(Z2)),Z2)
  Z0_new<-c(Z01,Z02)
  time<-c(xtime,ytime)
  delta<-c(xdelta,ydelta)
  joint<-c(rep(1,length(xdelta)),rep(2,length(ydelta)))
  b0_new<-rep(b0,times=2)
  b1_new<-c(b1,b2*0)
  b2_new<-c(b1*0,b2)
  return(list(time=time,delta=delta,joint=joint,b0=b0_new,b1=b1_new,b2=b2_new,Z1=Z1_new,Z2=Z2_new,Z0=Z0_new,mi=mi_vec,N=N))
}


# Other trivial functions
assemble=function(v){
  matrix(c(v[1]+v[2],v[1],v[1],v[1]+v[3]),nrow=2,ncol=2)
}

assemble_n=function(v){
  matrix(c(v[1]+v[2],-v[1],-v[1],v[1]+v[3]),nrow=2,ncol=2)
}

spread=function(m){
  c(m[1,1],m[2,2],m[1,2])
}

Bootstrap=function(data.complete,B=50,size.limit=200,raw=T,...){
  if(raw){
    q=nrow(as.matrix(data.complete$z1[[1]]))+nrow(as.matrix(data.complete$z2[[1]]))
    N=length(data.complete$z1)
    ### Bootstrap part 
    beta_b<-matrix(NA,nrow=B,ncol=q)
    D_b<-matrix(NA,nrow=B,ncol=3)
    b=1
    select.n<-min(N,size.limit)
    w=sqrt(select.n/N) # size adjustment, control the bootstrap sample to be <=200
    while(b<=B){
      resample<-sample.int(N,size=select.n,replace=T) #limit the resampled size to be 200
      data_b=list(x=data.complete$x[resample],y=data.complete$y[resample],z1=data.complete$z1[resample],z2=data.complete$z2[resample],
                  delta1=data.complete$delta1[resample],delta2=data.complete$delta2[resample])
      ppl_result_b<-tryCatch(BivPPL(data_b,raw=T,...), error=function(e) NULL)
      if(!is.null(ppl_result_b)){
        beta_b[b,]=ppl_result_b$beta_hat
        D_b[b,]=spread(ppl_result_b$D_hat)
        b=b+1
      }
    }
    return(list(beta_B_mean=colMeans(beta_b,na.rm = T),beta_B_SE=apply(beta_b,2,sd,na.rm=T)*w,
                D_B_mean=colMeans(D_b,na.rm = T),D_B_SE=apply(D_b,2,sd,na.rm=T)*w))
    ### 
  }else if('id' %in% names(data.complete)){
    id.ls<-unique(data.complete$id)
    N=length(unique(data.complete$id))
    if('z1' %in% names(data.complete)){
      q=ncol(data.complete$z1)+ncol(data.complete$z2)
      beta_b<-matrix(NA,nrow=B,ncol=q)
      D_b<-matrix(NA,nrow=B,ncol=3)
      b=1
      select.n<-min(N,size.limit)
      w=sqrt(select.n/N)
      while(b<=B){
        resample<-sample(id.ls,size=select.n,replace=T) #limit the resampled size to be 200
        select<-Reduce(c,lapply(1:select.n,function(i) which(resample[i]==id) ))
        mi_b<-sapply(1:select.n,function(i) data.complete$mi[which(resample[i]==id.ls)] )
        data_b=list(x=data.complete$x[select],y=data.complete$y[select],z1=data.complete$z1[select,],z2=data.complete$z2[select,],
                    delta1=data.complete$delta1[select],delta2=data.complete$delta2[select],N=select.n,mi=mi_b)
        ppl_result_b<-BivPPL(data_b,raw=F,max.iter=1000,eps=1e-6,huge=F,independence=F,D.initial=0.5)
        ppl_result_b<-tryCatch(PPL(data_b,raw=F,...), error=function(e) NULL)
        if(!is.null(ppl_result_b)){
          beta_b[b,]=ppl_result_b$beta_hat
          D_b[b,]=spread(ppl_result_b$D_hat)
          b=b+1
          print(b)
        }
      }
      return(list(beta_B_mean=colMeans(beta_b,na.rm = T),beta_B_SE=apply(beta_b,2,sd,na.rm=T)*w,
                  D_B_mean=colMeans(D_b,na.rm = T),D_B_SE=apply(D_b,2,sd,na.rm=T)*w))      
      
    }else if('Z1' %in% names(data.complete)){
      q=ncol(data.complete$Z1)+ncol(data.complete$Z2)
      beta_b<-matrix(NA,nrow=B,ncol=q)
      D_b<-matrix(NA,nrow=B,ncol=3)
      b=1
      select.n<-min(N,size.limit)
      w=sqrt(select.n/N)
      while(b<=B){
        resample<-sample(id.ls,size=select.n,replace=T) #limit the resampled size to be 200
        select<-match(resample,id)
        mi_b<-sapply(1:select.n,function(i) data.complete$mi[which(resample[i]==id.ls)] )
        data_b=list(time=data.complete$time[select],Z1=data.complete$Z1[select,],Z2=data.complete$Z2[select,],
                    delta=data.complete$delta[select],N=select.n,mi=mi_b)
        ppl_result_b<-tryCatch(BivPPL(data_b,raw=F,...), error=function(e) NULL)
        if(!is.null(ppl_result_b)){
          beta_b[b,]=ppl_result_b$beta_hat
          D_b[b,]=spread(ppl_result_b$D_hat)
          b=b+1
        }
      }
      return(list(beta_B_mean=colMeans(beta_b,na.rm = T),beta_B_SE=apply(beta_b,2,sd,na.rm=T)*w,
                  D_B_mean=colMeans(D_b,na.rm = T),D_B_SE=apply(D_b,2,sd,na.rm=T)*w))      
      
    }
  }else{
    stop("Invalid setting of \'raw\' or the input data")
  }
  
}

