library(MASS)
#assumes that no. of subjects and no. of items is divisible by 2.
gen_fake_lnorm<-function(nitem=16,
                         nsubj=40,
                         beta=c(6,-0.07),
                         ranefsd=c(0.25,0.12,0.18,0.0004),
                         corr=c(-.6,-.6),
                         sigma_e=0.51){
  ## prepare data frame for two condition latin square:
  g1<-data.frame(item=1:nitem,
                 cond=rep(letters[1:2],nitem/2))
  g2<-data.frame(item=1:nitem,
                 cond=rep(letters[2:1],nitem/2))
  
  ## assemble data frame:
  fakedat<-rbind(g1[rep(seq_len(nrow(g1)), 
                        nsubj/2),],
                 g2[rep(seq_len(nrow(g2)), 
                        nsubj/2),])
  
  ## add subjects:
  fakedat$subj<-rep(1:nsubj,each=nitem)
  
  ## add contrast coding:
  fakedat$x<-ifelse(fakedat$cond=="a",-1/2,1/2)
  
  ## Define variance covariance matrices:
  Sigma.u<-matrix(c(ranefsd[1]^2,
                    corr[1]*ranefsd[1]*ranefsd[2],
                    corr[1]*ranefsd[1]*ranefsd[2],
                    ranefsd[2]^2),nrow=2)
  
  Sigma.w<-matrix(c(ranefsd[3]^2,
                    corr[2]*ranefsd[3]*ranefsd[4],
                    corr[2]*ranefsd[3]*ranefsd[4],
                    ranefsd[4]^2),nrow=2)
  
  ## subj ranef
  u<-mvrnorm(n=length(unique(fakedat$subj)),
             mu=c(0,0),Sigma=Sigma.u)
  # item ranef
  w<-mvrnorm(n=length(unique(fakedat$item)),
             mu=c(0,0),Sigma=Sigma.w)

  ## generate data:  
  N<-dim(fakedat)[1]
  rt<-rep(NA,N)
  for(i in 1:N){
    rt[i] <- rlnorm(1,beta[1] + 
                      u[fakedat[i,]$subj,1] +
                      w[fakedat[i,]$item,1] + 
                      (beta[2]
                       +u[fakedat[i,]$subj,2]
                       +w[fakedat[i,]$item,2])*fakedat$x[i],sigma_e) 
  }   
  
  fakedat$rt<-rt
  fakedat
}
