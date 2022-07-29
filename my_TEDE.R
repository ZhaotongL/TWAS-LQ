my_TEDE_Sc=function(effect_GX=NULL,LDcov,YY,GY,Gmu,G,N,rcx1,rcx2=NULL,GG1=NULL,GG2=NULL,sigY2=NULL){
  GG = LDcov
  if(!is.null(effect_GX)){
  w2 = effect_GX
  XY=t(w2) %*% GY
  XX=t(w2) %*% GG %*% w2

  beta_TWAS2=ginv(XX) %*% XY  #checked
  sigY2=c(sqrt((YY+t(beta_TWAS2) %*% XX %*% beta_TWAS2 - 2* t(beta_TWAS2) %*% XY)/(N-1)))
  control=1
  } else{
      sigY2=sigY2
      control=0
  }
  UU2=GY-Gmu
  UUnew2=UU2/sigY2^2   #correct scores
  cov_UUnew2 = LDcov
  cov_UUnew1=cov_UUnew2/sigY2^2
  

  wal=t(UUnew2)%*%ginv(cov_UUnew1)%*%UUnew2 #summary
  #1-pchisq(wal,df=J-1)

  mow=c(1-pchisq(wal,df=qr(cov_UUnew2)$rank))

  #control=1
  if(control==1){
    cov_UU=GG/sigY2^2   #old
    if(is.null(nrow(effect_GX))){
    qusi=GG%*%rcx1%*%GG
    cov_UU2=cov_UU+as.numeric(beta_TWAS2^2)*qusi/sigY2^4   #consider var(w)
    }else{
        A1 = GG1 %*% rcx1 %*% t(GG1)
        A2 = GG2 %*% rcx2 %*% t(GG2)
        cov_UU2=(GG*sigY2^2 + as.numeric(beta_TWAS2[1]^2)*A1 + as.numeric(beta_TWAS2[2]^2)*A2)/sigY2^4
    }
    UUnew2=UU2/sigY2^2
    cov_UUnew2=cov_UU2
    wal2=t(UUnew2)%*%ginv(cov_UUnew2)%*%UUnew2
    mow=c(1-pchisq(wal,df=qr(cov_UUnew1)$rank),1-pchisq(wal2,df=qr(cov_UUnew2)$rank))
  }

  return(mow)
}

my_TEDE_aSPU=function(effect_GX=NULL,LDcov,YY,GY,Gmu,G,N,rcx1,rcx2=NULL,GG1=NULL,GG2=NULL,distribution_based=FALSE,sigY2=NULL){

  GG = LDcov
  if(!is.null(effect_GX)){
  w2 = effect_GX
  XY=t(w2) %*% GY
  XX=t(w2) %*% GG %*% w2

  beta_TWAS2=ginv(XX) %*% XY  #checked
  sigY2=c(sqrt((YY+t(beta_TWAS2) %*% XX %*% beta_TWAS2 - 2* t(beta_TWAS2) %*% XY)/(N-1)))
  control=1
  }else{
      sigY2=sigY2
      control=0
  }

  UU2=GY-Gmu
  UUnew2=UU2/sigY2^2   #correct scores
  cov_UUnew2 = LDcov
  cov_UUnew1=cov_UUnew2/sigY2^2

  aa=UUnew2/sqrt(diag(cov_UUnew1))
  bb=cov2cor(cov_UUnew1)

  if(distribution_based==FALSE){
    resultA0b=aSPUs(aa,bb,pow=c(1,2,4,8,Inf),n.perm=10000,prune=FALSE)
    ka=resultA0b$pvs[length(resultA0b$pvs)]
  }else{
    resultA0b=aSPUsD(aa,bb)
    ka=resultA0b[length(resultA0b)]
  }

  mow = c(ka)
  if(control==1){
    cov_UU=GG/sigY2^2   #old
    if(is.null(nrow(effect_GX))){
    qusi=GG%*%rcx1%*%GG
    cov_UU2=cov_UU+as.numeric(beta_TWAS2^2)*qusi/sigY2^4   #consider var(w)
    }else{
        A1 = GG1 %*% rcx1 %*% t(GG1)
        A2 = GG2 %*% rcx2 %*% t(GG2)
        cov_UU2=(GG*sigY2^2 + as.numeric(beta_TWAS2[1]^2)*A1 + as.numeric(beta_TWAS2[2]^2)*A2)/sigY2^4
    }
    UUnew2=UU2/sigY2^2
    cov_UUnew2=cov_UU2
  aa=UUnew2/sqrt(diag(cov_UUnew2)) #standardize -> Z
  bb=cov2cor(cov_UUnew2)
  if(distribution_based==FALSE){
    resultA0c=aSPUs(aa,bb,pow=c(1,2,4,8,Inf),n.perm=10000,prune=FALSE)
    ka2=resultA0c$pvs[length(resultA0c$pvs)]
  }else{
    resultA0c=aSPUsD(aa,bb)
    ka2=resultA0c[length(resultA0c)]
  }
  mow=c(ka,ka2)
  }

  #names(mow)=c("TEDE-aSPU","TEDE-aSPU2")
  return(mow)
}

