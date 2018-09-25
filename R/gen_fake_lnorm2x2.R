library(MASS)
#assumes that no. of subjects and no. of items is divisible by 2.
gen_fake_lnorm2x2<-function(nitem=16,
                         nsubj=40,
                         beta=c(6,-0.07,0.05,0.05),
                         subj_ranefsd=c(0.25,0.12,0.05,0.01),
                         item_ranefsd=c(0.18,0.0004,0.0001,0.0001),
                         subj_corr=rep(-0.6,6),
                         item_corr=rep(-0.6,6),
                         sigma_e=0.51){
  ## prepare data frame for four condition latin square:
  g1<-data.frame(item=1:nitem,
                 cond=rep(letters[1:4],nitem/4))
  g2<-data.frame(item=1:nitem,
                 cond=rep(letters[c(2,3,4,1)],nitem/4))
  g3<-data.frame(item=1:nitem,
                 cond=rep(letters[c(3,4,1,2)],nitem/4))
  g4<-data.frame(item=1:nitem,
                 cond=rep(letters[c(4,1,2,3)],nitem/4))
  
  
  ## assemble data frame:
  gp1<-g1[rep(seq_len(nrow(g1)), 
              nsubj/4),]
  gp2<-g2[rep(seq_len(nrow(g2)), 
              nsubj/4),]
  gp3<-g3[rep(seq_len(nrow(g3)), 
              nsubj/4),]
  gp4<-g4[rep(seq_len(nrow(g4)), 
              nsubj/4),]
  fakedat<-rbind(gp1,gp2,gp3,gp4)
  
  ## add subjects:
  fakedat$subj<-rep(1:nsubj,each=nitem)
  
  ## add contrast coding:
  ## main effect 1:
  fakedat$c1<-ifelse(fakedat$cond%in%c("a","b"),-1/2,1/2)
  ## main effect 2: 
  fakedat$c2<-ifelse(fakedat$cond%in%c("a","c"),-1/2,1/2)
  ## interaction:
  fakedat$c3<-ifelse(fakedat$cond%in%c("a","d"),-1/2,1/2)
  
  ## Define variance covariance matrices:
  Sigma.u <- diag(subj_ranefsd^2)
  Sigma.u[1,2]<-subj_corr[1]*subj_ranefsd[1]*subj_ranefsd[2]
  Sigma.u[1,3]<-subj_corr[2]*subj_ranefsd[1]*subj_ranefsd[3]
  Sigma.u[1,4]<-subj_corr[3]*subj_ranefsd[1]*subj_ranefsd[4]
  Sigma.u[2,3]<-subj_corr[4]*subj_ranefsd[2]*subj_ranefsd[3]
  Sigma.u[2,4]<-subj_corr[5]*subj_ranefsd[2]*subj_ranefsd[4]
  Sigma.u[3,4]<-subj_corr[6]*subj_ranefsd[3]*subj_ranefsd[4]
  
  Sigma.u[2,1]<-subj_corr[1]*subj_ranefsd[1]*subj_ranefsd[2]
  Sigma.u[3,1]<-subj_corr[2]*subj_ranefsd[1]*subj_ranefsd[3]
  Sigma.u[4,1]<-subj_corr[3]*subj_ranefsd[1]*subj_ranefsd[4]
  Sigma.u[3,2]<-subj_corr[4]*subj_ranefsd[2]*subj_ranefsd[3]
  Sigma.u[4,2]<-subj_corr[5]*subj_ranefsd[2]*subj_ranefsd[4]
  Sigma.u[4,3]<-subj_corr[6]*subj_ranefsd[3]*subj_ranefsd[4]

  # symmetric:  
  #Sigma.u==t(Sigma.u)
  
  Sigma.w<-diag(item_ranefsd^2) 
  Sigma.w[1,2]<-item_corr[1]*item_ranefsd[1]*item_ranefsd[2]
  Sigma.w[1,3]<-item_corr[2]*item_ranefsd[1]*item_ranefsd[3]
  Sigma.w[1,4]<-item_corr[3]*item_ranefsd[1]*item_ranefsd[4]
  Sigma.w[2,3]<-item_corr[4]*item_ranefsd[2]*item_ranefsd[3]
  Sigma.w[2,4]<-item_corr[5]*item_ranefsd[2]*item_ranefsd[4]
  Sigma.w[3,4]<-item_corr[6]*item_ranefsd[3]*item_ranefsd[4]
    
  Sigma.w[2,1]<-item_corr[1]*item_ranefsd[1]*item_ranefsd[2]
  Sigma.w[3,1]<-item_corr[2]*item_ranefsd[1]*item_ranefsd[3]
  Sigma.w[4,1]<-item_corr[3]*item_ranefsd[1]*item_ranefsd[4]
  Sigma.w[3,2]<-item_corr[4]*item_ranefsd[2]*item_ranefsd[3]
  Sigma.w[4,2]<-item_corr[5]*item_ranefsd[2]*item_ranefsd[4]
  Sigma.w[4,3]<-item_corr[6]*item_ranefsd[3]*item_ranefsd[4]
  
  #Sigma.w==t(Sigma.w)
  
  ## subj ranef
  u<-mvrnorm(n=length(unique(fakedat$subj)),
             mu=c(0,0,0,0),Sigma=Sigma.u)
  # item ranef
  w<-mvrnorm(n=length(unique(fakedat$item)),
             mu=c(0,0,0,0),Sigma=Sigma.w)

  ## generate data:  
  N<-dim(fakedat)[1]
  rt<-rep(NA,N)
  for(i in 1:N){
    rt[i] <- rlnorm(1,beta[1] + 
                      u[fakedat[i,]$subj,1] +
                      w[fakedat[i,]$item,1] + 
                      (beta[2]+u[fakedat[i,]$subj,2]+w[fakedat[i,]$item,2])*fakedat$c1[i]+
                      (beta[3]+u[fakedat[i,]$subj,3]+w[fakedat[i,]$item,3])*fakedat$c2[i]+
                      (beta[4]+u[fakedat[i,]$subj,4]+w[fakedat[i,]$item,4])*fakedat$c3[i],sigma_e) 
  }   
  
  fakedat$rt<-rt
  fakedat
}
