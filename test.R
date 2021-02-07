rhoLimit<-function(p){
  #The limit for pairwise correlations with marginal frequencies p
  if(!is.atomic(p) || typeof(p)!='double') return(NaN)
  if(sum((p<0) | (p>1))!=0) return(NaN)
  minPair<-outer(p, p, FUN='pmin')
  maxPair<-outer(p, p, FUN='pmax')
  maxRho<-sqrt((minPair/(1-minPair))/(maxPair/(1-maxPair)))
  minDev<-pmin(p, 1-p)	
  minDevOdds<-minDev/(1-minDev)
  minRho<-(-sqrt(minDevOdds%o%minDevOdds))
  return(list(maxRho=maxRho, minRho=minRho))	
}

rhoMaxEx<-function(p){
  #The upper bound for the pairwise correlation in exchange correlated binary distribution with marginal frequencies p
  if(!is.atomic(p) || typeof(p)!='double') return(NaN)
  if(sum((p<0) | (p>1))!=0) return(NaN)
  minP<-min(p)
  maxP<-max(p)
  rhoMax<-sqrt((minP/(1-minP))/(maxP/(1-maxP)))
  return(rhoMax)
}

rcbern0<-function(n, p, rho){
  #Generate Correlated Bernoulli Distribution
  #See "A Note on Generating Correlated Binary Variables"
  #A. D. Lunn and S. J. Davies
  #Biometrika
  #Vol. 85, No. 2 (Jun., 1998), pp. 487-490 
  if(!is.atomic(p) || typeof(p)!='double') return(NaN)
  if(length(p)!=1) return(NaN)
  if((p<0) | (p>1)) return(NaN)
  
  if((rho<0) || (rho>1)) return(NaN)
  Z<-rbinom(1, 1, p)
  Y<-rbinom(n, 1, p)
  U<-rbinom(n, 1, sqrt(rho))
  X<-(1-U)*Y+U*Z
  
  return(X)
}

rcbern<-function(n, p, rho){
  #Generate Correlated Bernoulli Distribution
  #See "A Note on Generating Correlated Binary Variables"
  #A. D. Lunn and S. J. Davies
  #Biometrika
  #Vol. 85, No. 2 (Jun., 1998), pp. 487-490 
  if(!is.atomic(p) || typeof(p)!='double') return(NaN)
  if(sum((p<0) | (p>1))!=0) return(NaN)
  
  m<-length(p)
  minP<-min(p)
  maxP<-max(p)
  rhoLimit<-sqrt((minP/(1-minP))/(maxP/(1-maxP)))
  if((rho<0) || (rho>rhoLimit)){
    cat(paste('The range of rho is [',0,',',round(rhoLimit,3),']\n', sep=''))
    cat('rho is out-of-range\n')
    return(NaN)
  }
  Pc<-sqrt(minP*maxP)/(sqrt(minP*maxP)+sqrt((1-minP)*(1-maxP)))
  Pa<-sqrt(rho*p*(1-p)/(Pc*(1-Pc)))
  if(max(Pa)==1){
    X<-replicate(n, rep(rbinom(1, 1, p[1]), m))
  }else{
    Pb<-(p-Pa*Pc)/(1-Pa)
    X<-replicate(n, {
      U<-rbinom(m, 1, Pa)
      Y<-rbinom(m, 1, Pb)
      Z<-rbinom(1, 1, Pc)
      (1-U)*Y+U*Z
    })
  }
  X<-t(X)
  
  return(X)
}


simGenoData<-function(n, f, rho){
  m<-length(f)
  X<-matrix(0, n, m)
  for(i in 1:m){
    if(rho!=0){
      #X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))		
      X[,i]<-rcbern0(n, f[i], rho)+rcbern0(n, f[i], rho)
    }else{
      X[,i]<-rbinom(n,2,f[i])
    }
  }
  return(X)
}


##########################################################
# Generating a random positive-definite matrix with user-specified positive eigenvalues
# If eigenvalues are not specified, they are generated from a uniform distribution
genCorrMatrix <- function (n){
  Z <- matrix(ncol=n, runif(n^2))
  Z <- t(Z)%*%Z
  Z <- Z+100*diag(n)
  return(cov2corr(Z))
}

cov2corr<-function(A){
  sqrtDiagA<-sqrt(diag(A))
  return(A/(sqrtDiagA%o%sqrtDiagA))
}

simGenoDataR<-function(n, f){
  library(bindata)
  m<-length(f)
  X<-matrix(0, n, m)
  rhoM<-genCorrMatrix(n)
  for(i in 1:m){
    X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=rhoM))		
  }
  rho<-mean(rhoM[lower.tri(rhoM)])
  return(list(X=X, rhoM=rhoM, rho=rho))
}

#startTime<-Sys.time()
#X<-simGenoData(20, runif(5,0.05, 0.5), 0.2)
#elapsedTime<-Sys.time()-startTime
#print(elapsedTime)

#Both cases and controls are correlated
simGwasData<-function(n0, n1, f0, f1, rho){
  n<-n0+n1
  m<-length(f0)
  X<-matrix(0, n, m)
  for(i in 1:m){
    if(rho!=0){
      #X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))		
      f<-c(rep(f0[i],n0), rep(f1[i], n1))
      X[,i]<-rcbern(1, f, rho)+rcbern(1, f, rho)  
    }else{
      X[,i]<-c(rbinom(n0, 2, f0[i]), rbinom(n1, 2, f1[i]))
    }
  }
  y<-c(rep(0, n0), rep(1, n1))
  return(list(X=X, y=y))
}

#No correlations between cases and controls
simGwasDataI<-function(n0, n1, f0, f1, rho){
  n<-n0+n1
  m<-length(f0)
  X<-matrix(0, n, m)
  for(i in 1:m){
    if(rho!=0){
      #X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))
      f0i<-rep(f0[i], n0)
      f1i<-rep(f1[i], n1)
      X[,i]<-c(rcbern(1, f0i, rho)+rcbern(1, f0i, rho), rcbern(1, f1i, rho)+rcbern(1, f1i, rho))
    }else{
      X[,i]<-c(rbinom(n0, 2, f0[i]), rbinom(n1, 2, f1[i]))
    }
  }
  y<-c(rep(0, n0), rep(1, n1))
  return(list(X=X, y=y))
}

#Different correlations used among cases and controls, but no correlations between a case and a control
simGwasDataII<-function(n0, n1, f0, f1, rho0, rho1){
  n<-n0+n1
  m<-length(f0)
  
  X<-matrix(0, n, m)
  for(i in 1:m){
    if(rho0!=0 || rho1!=0){
      #X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))
      f0i<-rep(f0[i], n0)
      f1i<-rep(f1[i], n1)
      X[,i]<-c(rcbern(1, f0i, rho0)+rcbern(1, f0i, rho0), rcbern(1, f1i, rho1)+rcbern(1, f1i, rho1))
    }
    else{
      X[,i]<-c(rbinom(n0, 2, f0[i]), rbinom(n1, 2, f1[i]))
    }
  }
  y<-c(rep(0, n0), rep(1, n1))
  #	return(list(X=X, y=y))
  return(X)
}


