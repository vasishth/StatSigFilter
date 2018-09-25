## ----setup,include=FALSE,cache=FALSE,echo=FALSE--------------------------
library(knitr)
library(coda)
library(plyr)
library(ggplot2)
library(xtable)
library(dplyr)
library(SIN)
library(papaja)

library(rstan)
library(parallel)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

library(loo)
library(lme4)

library(rstantools)

opts_chunk$set(fig.path='figures/fig-', 
               fig.align='center', 
               fig.show='hold',warning=FALSE,include=FALSE,message=FALSE)
options(replace.assign=TRUE,show.signif.stars=FALSE)
options(replace.assign=TRUE,width=75)
#opts_chunk$set(dev='postscript')

library(brms)

## ----loadfunctions, include=TRUE, echo=FALSE, warning=FALSE, error=TRUE, message=TRUE----
source("../R/createStanDat.R")
source("../R/createStanDatAcc.R")
source("../R/magnifytext.R")
source("../R/multiplot.R")
source("../R/plotresults.R")
source("../R/plotpredictions.R")
source("../R/stan_results.R")
source("../R/plotmeanSE.R")

## ----demotypeM,echo=FALSE,fig.width=7,fig.height=7,include=TRUE----------
set.seed(4321)
d<-15
sd<-100
lown<-power.t.test(d=d,sd=sd,power=.10,type="one.sample",alternative="two.sided",strict=TRUE)$n
highn<-power.t.test(d=d,sd=sd,power=.80,type="one.sample",alternative="two.sided",strict=TRUE)$n
nsim<-50
tlow<-thigh<-meanslow<-meanshigh<-CIuplow<-CIlwlow<-CIuphigh<-CIlwhigh<-NULL
critlow<-abs(qt(0.025,df=lown-1))
crithigh<-abs(qt(0.025,df=highn-1))

for(i in 1:nsim){
  x<-rnorm(lown,mean=d,sd=sd)
  meanslow[i]<-mean(x)
  tlow[i]<-t.test(x)$statistic
  CIuplow[i]<-mean(x)+critlow*sd(x)/sqrt(length(x))
  CIlwlow[i]<-mean(x)-critlow*sd(x)/sqrt(length(x))
  x<-rnorm(highn,mean=d,sd=sd)
  meanshigh[i]<-mean(x)
  thigh[i]<-t.test(x)$statistic
  CIuphigh[i]<-mean(x)+crithigh*sd(x)/sqrt(length(x))
  CIlwhigh[i]<-mean(x)-crithigh*sd(x)/sqrt(length(x))
}

 
siglow<-ifelse(abs(tlow)>abs(critlow),"p<0.05","p>0.05")
sighigh<-ifelse(abs(thigh)>abs(crithigh),"p<0.05","p>0.05")

summarylow<-data.frame(means=meanslow,significance=siglow, CIupper=CIuplow, CIlower=CIlwlow)
summaryhigh<-data.frame(index=1:nsim,means=meanshigh,significance=sighigh, CIupper=CIuphigh, CIlower=CIlwhigh)


# re-order data by mean effect size
summarylow<-summarylow[order(summarylow$means), ]
summarylow$index<-1:nrow(summarylow)
summaryhigh<-summaryhigh[order(summaryhigh$means), ]
summaryhigh$index<-1:nrow(summaryhigh)

p_low<-ggplot(summarylow, aes(y=means, x=index,
                              shape=significance,  ymax=CIupper, ymin=CIlower)) + 
  geom_pointrange()+
#  coord_flip()+
  geom_point(size=2.5)+
  scale_shape_manual(values=c(1, 17))+
  magnifytext(sze=22)+ 
  geom_hline(yintercept=15) +
  theme_bw() + 
    scale_x_continuous(name = "")+
  scale_y_continuous(name = "means",limits=c(-100,110))+
  labs(title="Effect 15 ms, SD 100, \n n=20, power=0.10")+
  theme(legend.position="none")+geom_hline(yintercept=0, linetype="dotted")

p_hi<-ggplot(summaryhigh, aes(y=means, x=index,
                              shape=significance, ymax=CIupper, ymin=CIlower)) + 
  geom_pointrange()+
#  coord_flip()+
  geom_point(size=2.5)+
  scale_shape_manual(values=c(1, 17))+
    scale_x_continuous(name = "Sample id")+
  magnifytext(sze=22)+ 
  geom_hline(yintercept=15) +
  theme_bw() + 
  scale_y_continuous(name = "means",limits=c(-100,110))+
  labs(title="Effect 15 ms, SD 100, \n n=350, power=0.80")+
  theme(legend.position=c(0.8,0.3))+geom_hline(yintercept=0, linetype="dotted")

multiplot(p_low,p_hi,cols=1)

## ----LKpredictions,cache=FALSE,echo=FALSE,include=TRUE,fig.width=8,fig.height=4----
LKsurprisal <- data.frame(
  pred = c(700,600,600,500),
  cond = c("a", "b", "c", "d")
)

LKmemory <- data.frame(
  pred = c(500,600,600,700),
  cond = c("a", "b", "c", "d")
)

plot1<-plotpredictions(dat=LKsurprisal,maintitle="Predictions of \n the expectation-based account")
plot2<-plotpredictions(dat=LKmemory,maintitle="Predictions of \n the memory account",ylabel="")
multiplot(plot1,plot2,cols=2)

## ----lkjvisual,cache=TRUE,echo=FALSE,include=FALSE-----------------------
fake_data <- list(x = rnorm(30000,0,1),N = 30000, R = 2) 

stancode <- "
data {
  int<lower=0> N; 
  real x[N]; 
  int R;
  }
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  x ~ normal(mu,sigma);  
}
generated quantities {
  corr_matrix[R] LKJ05;
  corr_matrix[R] LKJ1;
  corr_matrix[R] LKJ2;
  corr_matrix[R] LKJ4;
  LKJ05 = lkj_corr_rng(R,.1);
  LKJ1 = lkj_corr_rng(R,.5);
  LKJ2 = lkj_corr_rng(R,1);
  LKJ4 = lkj_corr_rng(R,2);
}
"

fitfake <- stan(model_code = stancode, pars = c("LKJ05","LKJ1","LKJ2","LKJ4"),
                data = fake_data, chains = 4, 
                iter = 2000)

corrs<-extract(fitfake,pars=c("LKJ05[1,2]","LKJ1[1,2]","LKJ2[1,2]","LKJ4[1,2]"))
corrs<-data.frame(corrs)
colnames(corrs)<-c("lkj05","lkj1","lkj2","lkj4")

## ----figpriors,echo=FALSE,include=TRUE,warning=FALSE,message=FALSE,fig.width=6,fig.height=4----

lkjplot1<-ggplot(corrs, aes(lkj05)) +
  geom_density(adjust=2)+xlab(expression(rho))+ggtitle("nu=0.1")+theme_bw()
lkjplot2<-ggplot(corrs, aes(lkj1)) +
  geom_density(adjust=3)+xlab(expression(rho))+ggtitle(expression(nu==0.5))+theme_bw()+ylim(0, 0.6)
lkjplot3<-ggplot(corrs, aes(lkj2)) +
  geom_density(adjust=3)+xlab(expression(rho))+ggtitle("nu=1")+theme_bw()
lkjplot4<-ggplot(corrs, aes(lkj2)) +
  geom_density(adjust=3)+xlab(expression(rho))+ggtitle(expression(nu==2))+theme_bw()+ylim(0, 0.6)
multiplot(lkjplot2,lkjplot4,cols=2)

## ----nonsigETdepmeasuresLKE2Replication,include=FALSE,eval=FALSE,echo=FALSE,message=FALSE,warning=FALSE----
## # Eyetracking Replication of LK Experiment 2
## # regression probability, skipping probability
## 
## # data ET Replication of LK Expt 2, reload to create a standalone chunk:
## dat<-read.table("../data/E4ETlevykellerExp2.txt",header=T)
## # subset data
## dat<-subset(dat,condition!="f" & condition!="p")
## 
## 
## # add column skipping probability (1 if word skipped, 0 otherwise)
## dat$SKP <- ifelse(dat$TFT>0,0,1)
## # add column regression probability (RBRC is no of regression from word n before continuing to the right.
## # 1 if regression occured, 0 otherwise.)
## dat$FPRP <- as.logical(dat$RBRC)
## 
## # contrast coding
## dat$dat<-ifelse(dat$condition%in%c("a","b"),1/2,-1/2)
## dat$adj<-ifelse(dat$condition%in%c("b","d"),-1/2,1/2)
## dat$int<-ifelse(dat$condition%in%c("b","c"),-1/2,1/2)
## 
## dat$roi<-factor(dat$roi)
## 
## region<-ifelse(dat$roi==23,"npacc",
##                ifelse(dat$roi==25,"verb",
##                       ifelse(dat$roi==27,"verb1","noncritical")))
## 
## dat$region<-region
## dat<-subset(dat,region!="noncritical")
## dat$region<-factor(dat$region)
## 
## # subset critical and postcritical region
## verb <-(subset(dat,region=="verb"))
## verb1 <-(subset(dat,region=="verb1"))
## 
## ## FPRP (first-pass regression prob) crit
## #summary(mFPRP <- glmer(FPRP~dat+adj+int+(1|subject)+(1|itemid), verb, family=binomial))
## 
## ## Skip prob crit
## #summary(mSKP <- glmer(SKP~dat+adj+int+(1|subject)+(1|itemid), verb, family=binomial))
## 
## 
## # FPRP crit
## subj <- as.integer(factor(verb$subject))
## N_subj <- length(unique(subj))
## item <- as.integer(factor(verb$itemid))
## N_items <- length(unique(item))
## 
## 
## X <- unname(model.matrix(~ 1 + dat + adj + int, verb))
## attr(X, which="assign") <- NULL
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb$FPRP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mFPRPStan <- stan(file = "StanModels/logitmaximal.stan",
##                data = stanDat,
##                iter = 2000,
##                chains = 4)
## 
## 
## # FPRP postcrit
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb1$FPRP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mFPRPpost <- stan(file = "StanModels/logitmaximal.stan",
##                   data = stanDat,
##                   iter = 2000,
##                   chains = 4)
## 
## 
## 
## # Skipping prob crit
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb$SKP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mSKPStan <- stan(file = "StanModels/logitmaximal.stan",
##                   data = stanDat,
##                   iter = 2000,
##                   chains = 4)
## 
## 
## 
## 
## # Skipping prob postcrit
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb1$SKP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mSKPpost <- stan(file = "StanModels/logitmaximal.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4,
##                  control = list(adapt_delta=0.99))
## 
## 

## ----originalLKdataExp1, include=FALSE, cache=FALSE,echo=FALSE,warning=FALSE----

##### ORIGINAL LEVY KELLER DATA EXPERIMENT 1 #####
# CRITICAL REGION (REGION 7)  #

reading_time <- read.table('../data/exp1_tt_r.res', header=TRUE)
#head(reading_time)

condition<-ifelse(reading_time$dat=="sub" & reading_time$adj=="sub","a",
                  ifelse(reading_time$dat=="sub" & reading_time$adj=="main","b",
                         ifelse(reading_time$dat=="main" & reading_time$adj=="sub","c", 
                                ifelse(reading_time$dat=="main" & reading_time$adj=="main","d","NA"))))

reading_time$condition<-factor(condition)


# contrast coding: 
reading_time$dat<-ifelse(reading_time$condition%in%c("a","b"),1/2,-1/2)
reading_time$adj<-ifelse(reading_time$condition%in%c("b","d"),-1/2,1/2)
reading_time$int<-ifelse(reading_time$condition%in%c("b","c"),-1/2,1/2)

                 ## ME DAT ## ME PP-ADJ ## INT
# a DAT-SC; PP-SC    0.5         0.5       0.5    
# b DAT-SC; PP-MC    0.5        -0.5      -0.5
# c DAT-MC; PP-SC   -0.5         0.5      -0.5
# d DAT-MC; PP-MC   -0.5        -0.5       0.5

# Thus, ME positive coefficient= longer RTs/slowdown when DAT bzw. PP ADJ in subordinate clause; 
# negative coefficient= faster RTs/speed-up when DAT bzw. PP ADJ in subordinate clause. 
# Interaction positve coefficient= longer RTs/slowdown when DAT/PP ADJ are in the same - subordinate OR main - clause. 

# remove zeros
reading_time_nozeros <- reading_time[reading_time$region7 != 0,]


# model (log tft) at the critical region (orgininally raw rts, see residuals)
#mLKE1crit  <- lmer(log(region7) ~ dat*adj + (dat*adj|subj) + (dat*adj|item), data=reading_time_nozeros)
#summary(mLKE1crit)

# I subset the data frame here (LK1 conds c and d only, such that I can bind it with conditions c and d from LK2).

dataLK1cd<-subset(reading_time, condition!="a" & condition!="b")

## ----AnalysisLK1critical,results='hide',include=FALSE,cache=TRUE,echo=FALSE,warning=FALSE----


# ANALYSIS of TRT original LK1 
# CRITICAL REGION (main verb)
stanDatLKE1<-createStanDat(d=reading_time_nozeros,
              rt=reading_time_nozeros$region7,
              form=as.formula("~1+dat+adj+int"))

LKE1 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE1,
                    iter = 2000, 
                    chains = 4)

LKE1_res<-stan_results(m=LKE1,params=c("Dat","Adj","DatxAdj"))

## original:
#interact  <- lmer(region7 ~ dat+adj+int + (dat+adj+int|subj) + (dat+adj+int|item), data=reading_time_nozeros)
#summary(interact)

## ----AnalysisLK1postcritical,include=FALSE,cache=TRUE,echo=FALSE,warning=FALSE----
# ANALYSIS of TRT original LK1 
# POSTCRITICAL REGION

reading_time_nozeros <- reading_time[reading_time$region8 != 0,]

#mLKE1post  <- lmer(log(region8) ~ dat*adj + (dat+adj|subj) + (dat+adj|item), data=reading_time_nozeros)
#model changed to above (failed to converge before (the original model))

stanDatLKE1post<-createStanDat(d=reading_time_nozeros,
              rt=reading_time_nozeros$region8,
              form=as.formula("~1+dat+adj+int"))

LKE1post <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE1post,
                    iter = 2000, 
                    chains = 4)

LKE1post_res<-stan_results(m=LKE1post,params=c("Dat","Adj","DatxAdj"))

## ----originalLKdataExp2,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE----

##### ORIGINAL LEVY KELLER DATA EXPERIMENT 2 #####

# CRITICAL REGION  (REGION 8) #

reading_time <- read.table('../data/exp3_tt_r.res', header=TRUE)


condition<-ifelse(reading_time$dat=="sub" & reading_time$adj=="sub","a",
                  ifelse(reading_time$dat=="sub" & reading_time$adj=="main","b",
                         ifelse(reading_time$dat=="main" & reading_time$adj=="sub","c", 
                                ifelse(reading_time$dat=="main" & reading_time$adj=="main","d","NA"))))

reading_time$condition<-factor(condition)


# contrast coding: 
reading_time$dat<-ifelse(reading_time$condition%in%c("a","b"),1/2,-1/2)
reading_time$adj<-ifelse(reading_time$condition%in%c("b","d"),-1/2,1/2)
reading_time$int<-ifelse(reading_time$condition%in%c("b","c"),-1/2,1/2)

                 ## ME DAT ## ME PP-ADJ ## INT
# a DAT-SC; PP-SC    0.5         0.5       0.5
# b DAT-SC; PP-MC    0.5        -0.5      -0.5
# c DAT-MC; PP-SC   -0.5         0.5      -0.5
# d DAT-MC; PP-MC   -0.5        -0.5       0.5


# exlude zeros from tfts
reading_time_nozeros <- reading_time[reading_time$region8 != 0,]

#mLKE2crit  <- lmer(log(region8) ~ dat*adj + (dat*adj|subj) + (dat*adj|item), data=reading_time_nozeros)
#summary(mLKE2crit)

# Subset the data frame (conds c and d only such that one can bind it with conditions c and d from LK1).
# LK2 c, d do not need to be relabelled

dataLK2cd<-subset(reading_time, condition!="a" & condition!="b")

## ----AnalysisLK2critical,results='hide',include=FALSE, cache=TRUE,echo=FALSE,warning=FALSE----

# ANALYSIS of TRT original LK2
# CRITICAL REGION (verb+aux)

stanDatLKE2<-createStanDat(d=reading_time_nozeros,
                       rt=reading_time_nozeros$region8,
                       form=as.formula("~1+dat+adj+int"))

#str(stanDat)
LKE2 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE2,
                    iter = 2000, 
                    chains = 4)

LKE2_res<-stan_results(m=LKE2,params=c("Dat","Adj","DatxAdj"))

## ----AnalysisLK2postcritical,include=FALSE,results='hide',cache=TRUE,echo=FALSE,warning=FALSE----

# EXP 2 REGION 9 POSTCRITICAL REGION #

reading_time_nozeros <- reading_time[reading_time$region9 != 0,]

#mLKE2post  <- lmer(log(region9) ~ dat*adj + (dat*adj|subj) + (dat+adj|item), data=reading_time_nozeros)
#summary(mLKE2post)


stanDatLKE2post<-createStanDat(d=reading_time_nozeros,
                       rt=reading_time_nozeros$region9,
                       form=as.formula("~1+dat+adj+int"))

LKE2post <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE2post,
                    iter = 2000, 
                    chains = 4)

LKE2post_res<-round(stan_results(m=LKE2post,
                           params=c("Dat","Adj","DatxAdj")))

## ----ResponseAccSPRLK1data,include=FALSE,echo=FALSE,warning=FALSE,message=FALSE----
dat<-read.table("../data/E1SPRlevykellerExp1.txt",header=T)

#head(dat)
#str(dat)

dat$item<-factor(dat$item)
dat$subj<-factor(dat$subj)

# subset acc data
datAcc<-subset(dat, dat$roi=="?")


# adjust column names:
datAcc <- plyr::rename(datAcc, c(word="response"))  # plyr

# "acc" column: 
datAcc <- plyr::rename(datAcc, c(region="acc"))
#str(datAcc)

# convert acc column to integer
datAcc$acc <- as.numeric(as.character(datAcc$acc))
#head(datAcc)

# subset acc data (filler items)
#datFillerAcc<-subset(datAcc, expt=="filler")
# mean acc for fillers
#mean(datFillerAcc$acc)

# subset acc data (exp items)
datAcc<-subset(datAcc, expt=="LKrep")
# mean accuracies exp items
#mean(datAcc$acc)

# compute means accuracies (condition)
#mean_acc<-round(tapply(datAcc$acc, datAcc$cond, mean),2)
#print(mean_acc)

# contrast coding
datAcc$dat<-ifelse(datAcc$cond%in%c("a","b"),1/2,-1/2)
datAcc$adj<-ifelse(datAcc$cond%in%c("b","d"),-1/2,1/2)
datAcc$int<-ifelse(datAcc$cond%in%c("b","c"),-1/2,1/2)
SPRLK1meanacc<-mean(100*with(datAcc,tapply(acc,cond,mean)),na.rm=TRUE)

## ----AnalysisAccSPRLK1,cache=TRUE,include=FALSE--------------------------

stanDatAccSPRLK1<-createStanDatAcc(d=datAcc,acc=datAcc$acc,form=as.formula("~1+dat+adj+int"))

SPRLK1Acc <- stan(file = "StanModels/responseaccmaximal.stan", 
                  data = stanDatAccSPRLK1,
                  iter = 2000, 
                  chains = 4)


AccSPRLK1_res<-stan_results(m=SPRLK1Acc,params=c("beta[2]","beta[3]","beta[4]"))

## ----DataSPRLK1,include=FALSE--------------------------------------------

# read in data SPR E1
dat<-read.table("../data/E1SPRlevykellerExp1.txt",header=T)

# subset experimental items (exlude filler, practice items)
E1spr<-subset(dat,roi!="?" & expt=="LKrep")

# in original LK1, postcrit region is "und so/und damit", hence, here sum rts of verb1+verb2

E1spr$region<-ifelse(E1spr$roi==13,"verb",
                     ifelse(E1spr$roi==14, "verb1", 
                            ifelse(E1spr$roi==15, "verb2", "noncritical")))


E1sprCRIT<-subset(E1spr, region=="verb")
E1sprPOST1<-subset(E1spr, region=="verb1")
E1sprPOST2<-subset(E1spr, region=="verb2")

E1sprCRIT<-E1sprCRIT[,c(1,3,4,7,8)]
E1sprPOST1<-E1sprPOST1[,c(1,3,4,7,8)]
E1sprPOST2<-E1sprPOST2[,c(1,3,4,7,8)]


E1sprPOST1$rt<-E1sprPOST1$rt + E1sprPOST2$rt

# inspect data frame/double-check the right rts were added
#E1sprPOST<-cbind(E1sprPOST1,E1sprPOST2)
#E1sprPOST$rtPOST<-E1sprPOST1$rt + E1sprPOST2$rt

E1spr<-rbind(E1sprCRIT,E1sprPOST1)

# contrast coding: 
E1spr$dat<-ifelse(E1spr$cond%in%c("a","b"),1/2,-1/2)
E1spr$adj<-ifelse(E1spr$cond%in%c("b","d"),-1/2,1/2)
E1spr$int<-ifelse(E1spr$cond%in%c("b","c"),-1/2,1/2)

                 ## ME DAT ## ME PP-ADJ ## INT
# a DAT-SC; PP-SC    0.5         0.5       0.5
# b DAT-SC; PP-MC    0.5        -0.5      -0.5
# c DAT-MC; PP-SC   -0.5         0.5      -0.5
# d DAT-MC; PP-MC   -0.5        -0.5       0.5


# subset critical, postcritical 
verb<-subset(E1spr,region=="verb")
verb1<-subset(E1spr,region=="verb1")

## ----AnalysisSPRLK1critical,cache=TRUE,include=FALSE---------------------

# E1 SPR rt at critical 
#head(verb)
stanDatSPRLK1<-createStanDat(d=verb,rt=verb$rt,
                             form=as.formula("~1+dat+adj+int"))

SPRLK1 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatSPRLK1,
                    iter = 2000, 
                    chains = 4)

# check
#m1 <- lmer(log(rt)~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb)
#summary(m1)

pars<-c("Dat","Adj","DatxAdj")
SPRE1_res<-round(stan_results(m=SPRLK1,params=pars))

## ----AnalysisSPRLK1criticalraw,echo=FALSE--------------------------------
#SPRLK1critm1 <- lmer(rt~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb)

## ----AnalysisSPRLK1postcritical,cache=TRUE,echo=FALSE,warning=FALSE------

# E1 SPR rt at postcritical
#head(verb1)
stanDatSPRE1post<-createStanDat(d=verb1,rt=verb1$rt,form=as.formula("~1+dat+adj+int"))

SPRE1post <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatSPRE1post,
                    iter = 2000, 
                    chains = 4)

# check
#m1post <- lmer(log(rt)~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb1)
#summary(m1post)

SPRE1post_res<-round(stan_results(m=SPRE1post,params=pars))

## ----AnalysisSPRLK1postcriticalraw,echo=FALSE----------------------------
#SPRLK1postcritm1post <- lmer(rt~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb1)

## ----AccE3data,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE--------
dat<-read.table("../data/E3SPRlevykellerExp2.txt",header=T)

dat$item<-factor(dat$item)
dat$subj<-factor(dat$subj)

# subset acc data
datAcc<-subset(dat, dat$roi=="?")

# adjust column names:
datAcc <- plyr::rename(datAcc, c(word="response"))  # plyr

# "acc" column: 
datAcc <- plyr::rename(datAcc, c(region="acc"))
#str(datAcc)

# convert acc column to integer
datAcc$acc <- as.numeric(as.character(datAcc$acc))
#head(datAcc)

# subset acc data (filler items)
#datFillerAcc<-subset(datAcc, expt=="filler")
# mean acc for fillers
#mean(datFillerAcc$acc)

# subset acc data (exp items)
datAcc<-subset(datAcc, expt=="LKrep")
# mean accuracies exp items
#mean(datAcc$acc)

# compute means accuracies (condition)
#mean_acc<-round(tapply(datAcc$acc, datAcc$cond, mean),2)
#print(mean_acc)

# contrast coding
datAcc$dat<-ifelse(datAcc$cond%in%c("a","b"),1/2,-1/2)
datAcc$adj<-ifelse(datAcc$cond%in%c("b","d"),-1/2,1/2)
datAcc$int<-ifelse(datAcc$cond%in%c("b","c"),-1/2,1/2)

SPRLK2meanacc<-100*mean(with(datAcc,tapply(acc,cond,mean)),na.rm=TRUE)

## ----AnalysisAccE3SPR,cache=TRUE,eval=FALSE,echo=FALSE,warning=FALSE-----
## stanDatAccSPRE2<-createStanDatAcc(d=datAcc,acc=datAcc$acc,form=as.formula("~1+dat+adj+int"))
## 
## E3Acc <- stan(file = "StanModels/responseaccmaximal.stan",
##                   data = stanDatAccSPRE2,
##                   iter = 2000,
##                   chains = 4)
## 
## 
## #print(mE3Acc, pars=c("beta","sigma_u","sigma_w"), probs=c(.025,.5,.975))

## ----DataRTE3crit,include=FALSE------------------------------------------

# # E3: ("direct replication of LK13 E2") Analysis
## Method: SPR (Linger)
## no participants or trials excluded
## N=28, items=24

dat<-read.table("../data/E3SPRlevykellerExp2.txt",header=T)

E3spr<-subset(dat,roi!="?" & expt=="LKrep")


# contrast coding (same as E1 above)
E3spr$dat<-ifelse(E3spr$cond%in%c("a","b"),1/2,-1/2)
E3spr$adj<-ifelse(E3spr$cond%in%c("b","d"),-1/2,1/2)
E3spr$int<-ifelse(E3spr$cond%in%c("b","c"),-1/2,1/2)

## subset critical region, postcritical region 
verb<-subset(E3spr,region=="verb")
verb1<-subset(E3spr,region=="verb1")
#head(verb)

# mean RTs at critical region:
#with(verb,round(tapply(rt,cond,mean)))
# mean RTs at postcritical region:
#with(verb1,round(tapply(rt,cond,mean)))


## ----AnalysisRTcritE3,cache=TRUE,echo=FALSE,include=FALSE----------------

stanDatSPRE2<-createStanDat(d=verb,rt=verb$rt,
              form=as.formula("~1+dat+adj+int"))

SPRE2 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatSPRE2,
                    iter = 2000, 
                    chains = 4)

# check
#m2 <- lmer(log(rt)~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb)
#summary(m2)

SPRE2_res<-round(stan_results(m=SPRE2,params=pars))

## ----AnalysisRTcritE3raw,echo=FALSE--------------------------------------
#mRTcritE3raw <- lmer(rt~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb)

## ----AnalysisRTE3postcrit,cache=TRUE,echo=FALSE,include=FALSE------------
stanDatSPRE2post<-createStanDat(d=verb1,
                                rt=verb1$rt,
              form=as.formula("~1+dat+adj+int"))

SPRE2post <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatSPRE2post,
                    iter = 2000, 
                    chains = 4)

SPRE2post_res<-round(stan_results(m=SPRE2post,params=pars))

## ----AnalysisE2Acc,include=FALSE-----------------------------------------
dat<-read.table("../data/E2ETlevykellerExp1.txt",header=T)

# subset Accuracy data
datAcc<-subset(dat,roi==1)
#summary(datAcc)
#head(datAcc)
#str(datAcc)
#levels(datAcc$condition)
#datAcc$RESPONSE_ACCURACY  
### quesitons only follow 50% of items: remove -2 (RESPONSE ACC) or 0 (answer) (i.e., no questions) and NAs


# subset data filler items (exclude LKrep and practice)
datFillerAcc<-subset(datAcc,condition=="f" & answer!=0)
datFillerAcc<-na.omit(datFillerAcc)

# subset data experimental items (exclude fillers and practice)
datAcc<-subset(datAcc,condition!="f" & condition!="p" & answer!=0)
datAcc<-na.omit(datAcc)

# sanity checks
#xtabs(~subject+condition,datAcc)
#length(unique(datAcc$itemid))


# "acc" column: 
datAcc <- plyr::rename(datAcc, c(RESPONSE_ACCURACY="acc"))
datFillerAcc <- plyr::rename(datFillerAcc, c(RESPONSE_ACCURACY="acc"))
# mean acc exp items
#mean(datAcc$acc) # [1] 0.6448598
# mean acc for fillers
#mean(datFillerAcc$acc) # [1] 0.898773

# compute means accuracies (condition)
#mean_acc<-round(tapply(datAcc$acc, datAcc$condition, mean),2)
#print(mean_acc)

# contrast coding
datAcc$dat<-ifelse(datAcc$cond%in%c("a","b"),1/2,-1/2)
datAcc$adj<-ifelse(datAcc$cond%in%c("b","d"),-1/2,1/2)
datAcc$int<-ifelse(datAcc$cond%in%c("b","c"),-1/2,1/2)

## ----E2AnalysisAcc,cache=TRUE,echo=FALSE,include=FALSE-------------------

subj <- as.numeric(as.factor(datAcc$subject))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(datAcc$itemid))
N_items <- length(unique(item))


X <- unname(model.matrix(~ 1 + dat + adj + int, datAcc))  
Z_u <- unname(model.matrix(~ 1, datAcc)) 
Z_w <- unname(model.matrix(~ 1, datAcc))

attr(X, which="assign") <- NULL              

# 2. Make Stan data (list)
stanDatAcc <- list(N = nrow(X), 
                P = ncol(X), 
                n_u = ncol(X),
                n_w = ncol(X),
                X = X,       
                Z_u = X,    
                Z_w = X,      
                J = N_subj, 
                K = N_items,
                acc = datAcc$acc,                  
                subj = as.integer(subj),  
                item = as.integer(item))  


# 3. Fit the model.

mE2Acc <- stan(file = "StanModels/responseaccmaximal.stan", 
                  data = stanDatAcc,
                  iter = 2000, 
                  chains = 4)


#print(mE2Acc, pars=c("beta","sigma_u","sigma_w"), probs=c(.025,.5,.975))

#Traceplot
#traceplot(mE2Acc, pars=c("beta","sigma_u","sigma_w"), inc_warmup=FALSE)

ETLK1meanacc<-100*mean(with(datAcc,tapply(acc,condition,mean)),na.rm=TRUE)

## ----DataE2,include=FALSE,echo=FALSE-------------------------------------
dat<-read.table("../data/E2ETlevykellerExp1.txt",header=T)

# exclude filler and practice items
dat<-subset(dat,condition!="f" & condition!="p")

# roi as factor column
dat$roi<-factor(dat$roi)
#unique(dat$roi)
#str(dat)

# contrast coding (same as SPR E1 and SPR E2)
dat$dat<-ifelse(dat$condition%in%c("a","b"),1/2,-1/2)
dat$adj<-ifelse(dat$condition%in%c("b","d"),-1/2,1/2)
dat$int<-ifelse(dat$condition%in%c("b","c"),-1/2,1/2)


#rois:
#precritical region (acc NP):
#a = 22 (was rois 22+23)
#b = 22 (was rois 22+23)
#c = 22 (was rois 22+23)
#d = 22 (was rois 22+23)

#critical region (verb):
#a = 24
#b = 24
#c = 24
#d = 24

#postcritical region ( "und"):
#a = 25 (was 25+26)
#b = 25 (was 25+26)
#c = 25 (was 25+26)
#d = 25 (was 25+26)

region<-ifelse(dat$roi==22,"npacc",
               ifelse(dat$roi==24,"verb",
                      ifelse(dat$roi==25,"verb1","noncritical")))

#length(region)  
#dim(dat)
#summary(dat)

dat$region<-region
dat<-subset(dat,region!="noncritical")
dat$region<-factor(dat$region)
#summary(dat)

dat$subj<-dat$subject
dat$item<-dat$itemid

# subset critical region, postcritical regionm
verb<-subset(dat,region=="verb")
verb1<-subset(dat,region=="verb1")

## ----E2AnalysisTFTcrit,include=FALSE,cache=TRUE,echo=FALSE---------------

verbTFT <- subset(verb, verb$TFT>0)
verbRRT <- subset(verb, verb$RRT>0)
verbFPRT <- subset(verb, verb$FPRT>0)

#summary(m_LKrepE1crit<-lmer(FPRT~dat+adj+int+(1+dat+adj+int||subject) + (1+dat+adj+int||itemid),verbFPRT))#summary(m_LKrepE1crit<-lmer(FPRT~dat+adj+int+(1+dat+adj+int||subject) + (1+dat+adj+int||itemid),verbFPRT))

#head(verb)

stanDat<-createStanDat(d=verbTFT,rt=verbTFT$TFT,
              form=as.formula("~1+dat+adj+int"))
#str(stanDatE2)

E2_3 <- stan(file = "StanModels/maxModel.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)


# check
#m6 <- lmer(log(TFT)~dat+adj+int+(1+dat+adj+int||subject)+(1+dat+adj+int||itemid),verbTFT)
#summary(m6)

E2_3_res<-stan_results(m=E2_3,params=pars)

## ----E2TFTcritraw,echo=FALSE---------------------------------------------
#mE2TFTcritraw <- lmer(TFT~dat+adj+int+(1+dat+adj+int||subject)+(1+dat+adj+int||itemid),verbTFT)

## ----E2AnalysisTFTpostcrit,cache=TRUE,echo=FALSE,include=FALSE-----------

verb1TFT <- subset(verb1, verb1$TFT>0)
verb1FPRT <- subset(verb1, verb1$FPRT>0)
verb1RRT <- subset(verb1, verb1$RRT>0)

#head(verb)
#summary(m_LKrepE1postcrit<-lmer(RRT~dat+adj+int+(1+dat+adj+int||subject) + (1+dat+adj+int||itemid),verb1RRT))
#summary(m_LKrepE1postcrit<-lmer(FPRT~dat+adj+int+(1+dat+adj+int||subject) + (1+dat+adj+int||itemid),verb1FPRT))
stanDat<-createStanDat(d=verb1TFT,rt=verb1TFT$TFT,
              form=as.formula("~1+dat+adj+int"))

E2_3post <- stan(file = "StanModels/maxModel.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)


#print(E2_3post, pars=c("beta","sigma_e","sigma_u","sigma_w"), probs=c(.025,.5,.975))

E2_3post_res<-stan_results(m=E2_3post,params=pars)

## ----E4AccData,cache=FALSE,echo=FALSE,include=FALSE----------------------
dat<-read.table("../data/E4ETlevykellerExp2.txt",header=T)

# subset Accuracy data
datAcc<-subset(dat,roi==1)
### quesitons only follow 50% of items: remove -2 (RESPONSE ACC) or 0 (answer) (i.e., no questions) and NAs


# subset data filler items (exclude LKrep and practice)
datFillerAcc<-subset(datAcc,condition=="f" & answer!=0)
datFillerAcc<-na.omit(datFillerAcc)

# subset data experimental items (exclude fillers and practice)
datAcc<-subset(datAcc,condition!="f" & condition!="p" & answer!=0)
datAcc<-na.omit(datAcc)

# sanity checks
#
#xtabs(~subject+condition,datAcc)
#length(unique(datAcc$itemid))


# "acc" column: 
datAcc <- plyr::rename(datAcc, c(RESPONSE_ACCURACY="acc"))
datAcc$acc<-as.numeric(as.character(datAcc$acc))
#str(datAcc)
datFillerAcc <- plyr::rename(datFillerAcc, c(RESPONSE_ACCURACY="acc"))
datFillerAcc$acc<-as.numeric(as.character(datFillerAcc$acc))

# mean acc exp items
#mean(datAcc$acc)
# mean acc for fillers
#mean(datFillerAcc$acc) 

# compute means accuracies (condition)
#mean_acc<-round(tapply(datAcc$acc, datAcc$condition, mean),2)
#print(mean_acc)

# contrast coding
datAcc$dat<-ifelse(datAcc$cond%in%c("a","b"),1/2,-1/2)
datAcc$adj<-ifelse(datAcc$cond%in%c("b","d"),-1/2,1/2)
datAcc$int<-ifelse(datAcc$cond%in%c("b","c"),-1/2,1/2)

## ----E4AnalysisAccuracy,cache=TRUE,echo=FALSE,include=FALSE--------------

subj <- as.numeric(as.factor(datAcc$subject))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(datAcc$itemid))
N_items <- length(unique(item))


X <- unname(model.matrix(~ 1 + dat + adj + int, datAcc))  
Z_u <- unname(model.matrix(~ 1, datAcc)) 
Z_w <- unname(model.matrix(~ 1, datAcc))

attr(X, which="assign") <- NULL              

# 2. Make Stan data (list)
stanDatAcc <- list(N = nrow(X), 
                P = ncol(X), 
                n_u = ncol(X),
                n_w = ncol(X),
                X = X,       
                Z_u = X,    
                Z_w = X,      
                J = N_subj, 
                K = N_items,
                acc = datAcc$acc,                  
                subj = as.integer(subj),  
                item = as.integer(item))  


# 3. Fit the model.

mE4Acc <- stan(file = "StanModels/responseaccmaximal.stan", 
                  data = stanDatAcc,
                  iter = 2000, 
                  chains = 4)


#print(mE4Acc, pars=c("beta","sigma_u","sigma_w"), probs=c(.025,.5,.975))

#Traceplot
#traceplot(mE4Acc, pars=c("beta","sigma_u","sigma_w"), inc_warmup=FALSE)

ETLK2meanacc<-100*mean(with(datAcc,tapply(acc,condition,mean)),na.rm=TRUE)

## ----E4RTdata,include=FALSE,echo=FALSE-----------------------------------

# E4: ("identical replication of LK13 E2") Analysis
## Method: ET (EyeLink 1000, SR Research)
## Item 14 excluded from lists 1, 5, 10, 11, 16, 18, 26 due to spelling error (item with error was analyzed in original Levy&Keller for all subjects). None of the other 5 replications are affected.
## few missing data points due to skipping of trials (i.e. participant looked at fixation trigger immediately without reading sentence).
## N=28, items=24

dat<-read.table("../data/E4ETlevykellerExp2.txt",header=T)

# subset data
dat<-subset(dat,condition!="f" & condition!="p")

# contrast coding 
dat$dat<-ifelse(dat$condition%in%c("a","b"),1/2,-1/2)
dat$adj<-ifelse(dat$condition%in%c("b","d"),-1/2,1/2)
dat$int<-ifelse(dat$condition%in%c("b","c"),-1/2,1/2)


dat$roi<-factor(dat$roi)
#unique(dat$roi)
#unique(dat$subject)

#rois:
#precritical region (acc NP):
#a = 23 (was rois 23+24)
#b = 23 (was rois 23+24)
#c = 23 (was rois 23+24)
#d = 23 (was rois 23+24)

#critical region (verb):
#a = 25 (was rois 25+26)
#b = 25 (was rois 25+26)
#c = 25 (was rois 25+26)
#d = 25 (was rois 25+26)

#postcritical region (verb1):
#a = 27 (was rois 27+28) * for all conditions (except for items 15,17,18,20; only roi 27, only 1 word)
#b = 27 (was rois 27+28)
#c = 27 (was rois 27+28)
#d = 27 (was rois 27+28)

region<-ifelse(dat$roi==23,"npacc",
               ifelse(dat$roi==25,"verb",
                      ifelse(dat$roi==27,"verb1","noncritical")))

dat$region<-region
dat<-subset(dat,region!="noncritical")
dat$region<-factor(dat$region)

dat$subj<-dat$subject
dat$item<-dat$itemid

# subset critical and postcritical region
verb <-subset(dat,region=="verb")
verb1 <-subset(dat,region=="verb1")

## ----AnalysisE4TFTcrit,cache=TRUE,echo=FALSE,include=FALSE---------------
verbTFT <- subset(verb, verb$TFT>0)
verbRRT <- subset(verb, verb$RRT>0)
verb1RRT <- subset(verb1, verb1$RRT>0)
verb1TFT <- subset(verb, verb1$TFT>0)

#summary(m_LKE2repRRT<-lmer(RRT~dat+adj+int+(1+dat+adj+int||subject)+(1|item),verb1RRT))
# 

## regression prob: does not converge
verb$reg<-ifelse(verb$FPRT==verb$RPD,0,1)
verb1$reg<-ifelse(verb1$FPRT==verb1$RPD,0,1)
#with(verb1,tapply(reg,condition,mean))
#summary(m_LKE2repregp<-glmer(reg~dat+adj+int+(1+dat+adj+int|subj)+(1|item),family=binomial(),verb1))
#summary(m_LKE2repskip<-glmer(FFP~dat+adj+int+(1|subj)+(1|item),family=binomial(),verb1))

#head(verb)

stanDat<-createStanDat(d=verbTFT,rt=verbTFT$TFT,
              form=as.formula("~1+dat+adj+int"))

E4_3 <- stan(file = "StanModels/maxModel.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)

E4_3_res<-stan_results(m=E4_3,params=pars)

## ----E4ETraw,echo=FALSE--------------------------------------------------
#m_E4ETTFTrawcrit<-lmer(TFT~dat+adj+int+(1+dat+adj+int||subj)+(1|item),verbTFT)
#m_E4ETTFTrawpostcrit<-lmer(TFT~dat+adj+int+(1+dat+adj+int||subj)+(1|item),verbTFT)

## ----AnalysisE4TFTpostcrit,cache=TRUE,echo=FALSE,include=FALSE-----------
verb1TFT <- subset(verb1, verb1$TFT>0)

stanDat<-createStanDat(d=verb1TFT,rt=verb1TFT$TFT,
              form=as.formula("~1+dat+adj+int"))

E4_3post <- stan(file = "StanModels/maxModel.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)

E4_3post_res<-stan_results(m=E4_3post,params=pars)

## ----figurese1,echo=FALSE,include=FALSE,warning=FALSE,message=FALSE------
cond<-factor(c("Dat","Adj","DatxAdj"),levels=c("Dat","Adj","DatxAdj"))
cond<-rep(cond,3)
expt<-factor(c(rep("LK Expt 1",3),
        rep("Expt 1, SPR",3),rep("Expt 2, ET",3)),levels=c("LK Expt 1","Expt 1, SPR","Expt 2, ET"))
  
rownames(LKE1_res)<-NULL
rownames(SPRE1_res)<-NULL
rownames(E2_3_res)<-NULL

LKE1all<-data.frame(rbind(LKE1_res,SPRE1_res,E2_3_res))
LKE1all<-data.frame(cond=cond,expt=expt,LKE1all)

p1<-plotresults(LKE1all,xaxisblank=TRUE,
                maintitle="Replications of LK Expt 1 \n (critical region)",removelegend=TRUE)

rownames(LKE1post_res)<-NULL
rownames(SPRE1post_res)<-NULL
rownames(E2_3post_res)<-NULL

LKE1postall<-data.frame(rbind(LKE1post_res,SPRE1post_res,E2_3post_res))
LKE1postall<-data.frame(cond=cond,expt=expt,LKE1postall)

p2<-plotresults(LKE1postall,maintitle="Replications of LK Expt 1 \n (post-critical region)",ylabel="",xaxisblank=TRUE)
#multiplot(p1,p2,cols=1)

## ----figurese2,echo=FALSE,warning=FALSE,message=FALSE,include=TRUE,fig.width=7,fig.height=7----
rownames(LKE2_res)<-NULL
rownames(SPRE2_res)<-NULL
rownames(E4_3_res)<-NULL

LKE2all<-data.frame(rbind(LKE2_res,SPRE2_res,E4_3_res))

expt<-factor(c(rep("LK Expt 2",3),
        rep("Expt 3, SPR",3),rep("Expt 4, ET",3)),levels=c("LK Expt 2","Expt 3, SPR","Expt 4, ET"))


LKE2all<-data.frame(cond=cond,expt=expt,LKE2all)

p1_2<-plotresults(LKE2all,maintitle="Replications of LK Expt 2 \n (critical region)",cols=c("darkgray", "black","black"),removelegend=TRUE)

rownames(LKE2post_res)<-NULL
rownames(SPRE2post_res)<-NULL
rownames(E4_3post_res)<-NULL

LKE2postall<-data.frame(rbind(LKE2post_res,SPRE2post_res,E4_3post_res))
LKE2postall<-data.frame(cond=cond,expt=expt,LKE2postall)

p2_2<-plotresults(LKE2postall,maintitle="Replications of LK Expt 2 \n (post-critical region)",cols=c("darkgray", "black","black"),ylabel="")
multiplot(p1,p1_2,p2,p2_2,cols=2)

## ----CombinedOriginalLK1LK2cd,cache=FALSE,include=FALSE,echo=FALSE,warning=FALSE----

# Prepare data (which was subset earlier as dataLK1cd and dataLK2cd) to merge as one data frome.
#head(dataLK1cd)
#head(dataLK2cd)

## LK EXPERIMENT 1

# relabel LK1 c, d as conditions a and b
dataLK1cd$condition[dataLK1cd$condition=="c"] <- "a"
dataLK1cd$condition[dataLK1cd$condition=="d"] <- "b"

dataLK1cd$subj<-as.factor(dataLK1cd$subj)
dataLK1cd$exp<-"LK1"
dataLK1cd$exp<-as.factor(dataLK1cd$exp)
#str(dataLK1cd)

# add new subject column that makes subj "1" "1LK1"
dataLK1cd$subject<-with(dataLK1cd, paste(subj,exp, sep=""))
dataLK1cd$subject<-as.factor(dataLK1cd$subject)
dataLK1cd<-dataLK1cd[,c(2,5,6,7,8,9,10,11,12,13,14,17)]
#head(dataLK1cd)

## LK EXPERIMENT 2

# recode such that we have the same critcal/postcritical regions for LK1 and LK2 in my merged data frame
# ==> for LK2, delete region 7

dataLK2cd <- dataLK2cd[,-11]

# ==> for LK2 rename region 8: now region 7
# ==> for LK2 rename region 9: now region 8

names(dataLK2cd)[names(dataLK2cd)=="region8"] <- "region7"
names(dataLK2cd)[names(dataLK2cd)=="region9"] <- "region8"
names(dataLK2cd)[names(dataLK2cd)=="region10"] <- "region9"


dataLK2cd$subj<-as.factor(dataLK2cd$subj)
dataLK2cd$exp<-"LK2"
dataLK2cd$exp<-as.factor(dataLK2cd$exp)
#str(dataLK2cd)

# add new subject column that makes subj "1" "1LK2"
dataLK2cd$subject<-with(dataLK2cd, paste(subj,exp, sep=""))
dataLK2cd$subject<-as.factor(dataLK2cd$subject)
dataLK2cd<-dataLK2cd[,c(2,5,6,7,8,9,10,11,12,13,14,17)]


#head(dataLK1cd)
#head(dataLK2cd)
#str(dataLK1cd)

## bind LK1 (c,d) and LK2 (c,d) together as one data frame:

dataLK1LK2cd<-rbind(dataLK1cd,dataLK2cd)
#head(dataLK1LK2cd)
#str(dataLK1LK2cd)

## contrast coding
dataLK1LK2cd$load<-ifelse(dataLK1LK2cd$condition%in%c("a","b"),-1/2,1/2)
dataLK1LK2cd$dist<-ifelse(dataLK1LK2cd$condition%in%c("a","c"),-1/2,1/2)
dataLK1LK2cd$int<-ifelse(dataLK1LK2cd$condition%in%c("a","d"),-1/2,1/2)


                         ## ME LOAD ## ME DIST  ## INT
# a DAT-MC;     PP-SC      -0.5        -0.5      -0.5   (originally E1 LK13 cond c)
# b DAT-MC;     PP-MC      -0.5         0.5       0.5   (orininally E1 LK13 cond d)
# c DAT-MC emb; PP-SC emb   0.5        -0.5       0.5   (originally E2 LK13 cond c)
# d DAT-MC emb; PP-MC emb   0.5         0.5      -0.5   (originally E2 LK13 cond d)

# ME load: positive coef = longer RTs in higher memory load conditions (MC embedding in RC); negative coefficient = shorter RTs (speedup) in high memory load conditions
# ME dist: positive coef = longer RTs when distance btw. subject and verb increased, i.e. when both ADJ and DAT in MC independent of whether MC is embedded in RC or not. 

# Region 7 (merged data frame) = critical region, 
# Region 8 = postcritical region 



## ----AnalysisLK1LK2critical,cache=TRUE,echo=FALSE,include=FALSE,warning=FALSE,message=FALSE----

# original LK1 and LK2 merged (conditions c, d)

# subset data
dataLK1LK2cdCRIT<-subset(dataLK1LK2cd, dataLK1LK2cd$region7>0)
#head(dataLK1LK2cdCRIT)
## load is a between subjects factor:
#xtabs(~subj+load,dataLK1LK2cdCRIT)
## but all predictors are within items
##xtabs(~item+load,dataLK1LK2cdCRIT)

subj <- as.numeric(as.factor(dataLK1LK2cdCRIT$subject))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(dataLK1LK2cdCRIT$item))
N_items <- length(unique(item))

X <- unname(model.matrix(~ 1 + load + dist + int, dataLK1LK2cdCRIT))
Z_u <- unname(model.matrix(~ 1 + dist + int, dataLK1LK2cdCRIT))
Z_w <- unname(model.matrix(~ 1 + load + dist + int, dataLK1LK2cdCRIT))

attr(X, which="assign") <- NULL

# 2. Make Stan data (list)
stanDat <- list(N = nrow(X),          
                P = ncol(X),               
                n_u = ncol(Z_u),              
                n_w = ncol(Z_w),       
                X = X,                      
                Z_u = Z_u,              
                Z_w = Z_w,                   
                J = N_subj,                
                K = N_items,
                rt = dataLK1LK2cdCRIT$region7,                   
                subj = as.integer(subj),
                item = as.integer(item))
#str(stanDat)        

# 3. Fit the model.

LKmerged <- stan(file = "StanModels/maxModelmerged.stan", 
                    data = stanDat,
                    iter = 2000, 
                    chains = 4,
                    control = list(adapt_delta = 0.99)) 

# added adapt delta as got divergent transitions warning once, and low neff/high Rhat another time. 
# sigmau1 and sigmau2 still low neffs.
# increase iteration to 4000

#print(LKmerged, pars=c("beta","sigma_e","sigma_u","sigma_w"), probs=c(.025,.5,.975))
#print(LKmerged, pars=c("Load","Dist","LoadxDist"), probs=c(.025,.5,.975))

# Traceplot
#traceplot(m_LKmerged, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)
# sigmau1 and sigmau2 ?

pars<-c("Load","Dist","LoadxDist")
# extract the model estimates
LKmerged_res<-stan_results(m=LKmerged,params=pars)

## ----AnalysisLK1LK2postcritical,cache=TRUE,echo=FALSE,include=FALSE,message=FALSE----

# subset data
dataLK1LK2cdPOST<-subset(dataLK1LK2cd, dataLK1LK2cd$region8>0)
#head(dataLK1LK2cdPOST)

subj <- as.numeric(as.factor(dataLK1LK2cdPOST$subject))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(dataLK1LK2cdPOST$item))
N_items <- length(unique(item))


X <- unname(model.matrix(~ 1 + load + dist + int, dataLK1LK2cdPOST))
Z_u <- unname(model.matrix(~ 1 + dist + int, dataLK1LK2cdPOST))
Z_w <- unname(model.matrix(~ 1 + load + dist + int, dataLK1LK2cdPOST))

attr(X, which="assign") <- NULL

# 2. Make Stan data (list)
stanDat <- list(N = nrow(X),               
                P = ncol(X),               
                n_u = ncol(Z_u),           
                n_w = ncol(Z_w),             
                X = X,                      
                Z_u = Z_u,               
                Z_w = Z_w,                   
                J = N_subj,                 
                K = N_items,
                rt = dataLK1LK2cdPOST$region8,                   
                subj = as.integer(subj),
                item = as.integer(item))

# 3. Fit the model.

### 
LKmergedpost <- stan(file = "StanModels/maxModelmerged.stan", 
                    data = stanDat,
                    iter = 2000, 
                    chains = 4,
                    control = list(adapt_delta = 0.99))


#print(LKmergedpost, pars=c("beta","sigma_e","sigma_u","sigma_w"), probs=c(.025,.5,.975))

# Traceplot
#traceplot(m_LKmergedpost, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)

# extract the model estimates
LKmergedpost_res<-stan_results(m=LKmergedpost,params=pars)

## ----E5ResponseAccdata,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE----
dat<-read.table("../data/E5SPRlevykellerExp12.txt",header=T)

dat$item<-factor(dat$item)
dat$subj<-factor(dat$subj)

# subset acc data
datAcc<-subset(dat, dat$roi=="?")

# adjust column names:
datAcc <- plyr::rename(datAcc, c(word="response"))  # plyr

# "acc" column: 
datAcc <- plyr::rename(datAcc, c(region="acc"))
#str(datAcc)

# convert acc column to integer
datAcc$acc <- as.numeric(as.character(datAcc$acc))
#head(datAcc)

# subset acc data (filler items)
#datFillerAcc<-subset(datAcc, expt=="filler")
# mean acc for fillers
#mean(datFillerAcc$acc)

# subset acc data (exp items)
datAcc<-subset(datAcc, expt=="LKrep")
# mean accuracies exp items
#mean(datAcc$acc)

# compute means accuracies (condition)
#mean_acc<-round(tapply(datAcc$acc, datAcc$cond, mean),2)
#print(mean_acc)

# contrast coding
datAcc$load<-ifelse(datAcc$cond%in%c("a","b"),-1/2,1/2)
datAcc$dist<-ifelse(datAcc$cond%in%c("a","c"),-1/2,1/2)
datAcc$int<-ifelse(datAcc$cond%in%c("a","d"),-1/2,1/2)

# "combined exp, conds c, d of LK1 and LK2

#              load    dist    int
# a            -0.5   -0.5    -0.5
# b            -0.5   +0.5    +0.5 
# c            +0.5   -0.5    +0.5
# d            +0.5   +0.5    -0.5

E5meanacc<-round(mean(100*with(datAcc,tapply(acc,cond,mean)),na.rm=TRUE))

## ----AnalysisResponseAccE5SPR, include=FALSE,results='hide',cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE----

subj <- as.numeric(as.factor(datAcc$subj))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(datAcc$item))
N_items <- length(unique(item))


X <- unname(model.matrix(~ 1 + load + dist + int, datAcc))  
Z_u <- unname(model.matrix(~ 1, datAcc)) 
Z_w <- unname(model.matrix(~ 1, datAcc))

attr(X, which="assign") <- NULL              

# 2. Make Stan data (list)
stanDatAcc <- list(N = nrow(X), 
                P = ncol(X), 
                n_u = ncol(X),
                n_w = ncol(X),
                X = X,       
                Z_u = X,    
                Z_w = X,      
                J = N_subj, 
                K = N_items,
                acc = datAcc$acc,                  
                subj = as.integer(subj),  
                item = as.integer(item))  


# 3. Fit the model.

mE5Acc <- stan(file = "StanModels/responseaccmaximal.stan", 
                  data = stanDatAcc,
                  iter = 2000, 
                  chains = 4)


#print(mE5Acc, pars=c("beta","sigma_u","sigma_w"), probs=c(.025,.5,.975))

#Traceplot
#traceplot(mE5Acc, pars=c("beta","sigma_u","sigma_w"), inc_warmup=FALSE)


## ----AnalysisE5SPRcrit,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE----

# E5: LK13 combined ("replication of LK13 E1 and E2 conds c, d") Analysis
## Method: SPR (Linger)
## no participants or trials excluded
## N=28, items=24

dat<-read.table("../data/E5SPRlevykellerExp12.txt",header=T)

E5spr<-subset(dat,roi!="?" & expt=="LKrep")

#round(with(subset(E5spr,region=="verb1"),tapply(rt,cond,mean)))
#round(with(subset(E5spr,region=="verb2"),tapply(rt,cond,mean)))


# ROIs: crit: "verb" (a,b = "versteckt, c,d = "versteckt hat"), 
# postcrit: "verb1" (a,b = "und damit", c,d = "den Besuch")
# postcrit for a,b not presented together in Linger, therefore, merge here (verb1+verb2)

#datab<-subset(E5spr,cond!="c" & cond!="d")

#E5sprCRIT<-subset(datab,region=="verb")
#E5sprPOST1<-subset(datab,region=="verb1")

#E5sprPOST2<-subset(datab,region=="verb2")

#E5sprCRIT<-E5sprCRIT[,c(1,3,4,7,8)]
#E5sprPOST1<-E5sprPOST1[,c(1,3,4,7,8)]
#E5sprPOST2<-E5sprPOST2[,c(1,3,4,7,8)]

#E5sprPOST1$rt<-E5sprPOST1$rt + E5sprPOST2$rt

#datab<-rbind(E5sprCRIT,E5sprPOST1)

#datcd<-subset(E5spr,cond!="a" & cond!="b")
#datcd<-subset(datcd,region==c("verb","verb1"))   # "verb1 region",e.g. "den_Konkurs" was shown as one roi in Linger
#datcd<-datcd[,c(1,3,4,7,8)]

#E5spr<-rbind(datab,datcd)

#boxplot(log(rt)~cond,subset(E5spr,region=="verb1"))


# contrast coding
E5spr$load<-ifelse(E5spr$cond%in%c("a","b"),-1/2,1/2)
E5spr$dist<-ifelse(E5spr$cond%in%c("a","c"),-1/2,1/2)
E5spr$int<-ifelse(E5spr$cond%in%c("a","d"),-1/2,1/2)

                         ## ME LOAD ## ME DIST  ## INT
# a DAT-MC;     PP-SC      -0.5        -0.5      -0.5   (originally E1 LK13 cond c)
# b DAT-MC;     PP-MC      -0.5         0.5       0.5   (orininally E1 LK13 cond d)
# c DAT-MC emb; PP-SC emb   0.5        -0.5       0.5   (originally E2 LK13 cond c)
# d DAT-MC emb; PP-MC emb   0.5         0.5      -0.5   (originally E2 LK13 cond d)



# subset crit, postcrit
verb<-subset(E5spr,region=="verb")
verb1<-subset(E5spr,region=="verb1")



## ----AnalysisE5RTcrit,cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE----

subj <- as.numeric(as.factor(verb$subj))
#unique(verb$subj)
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(verb$item))
#unique(verb$item)
N_items <- length(unique(item))

X <- unname(model.matrix(~ 1 + load + dist + int, verb))
Z_u <- unname(model.matrix(~ 1, verb))
Z_w <- unname(model.matrix(~ 1, verb))

attr(X, which="assign") <- NULL

# 2. Make Stan data (list)
stanDat <- list(N = nrow(X),           
                P = ncol(X),               
                n_u = ncol(X),            
                n_w = ncol(X),              
                X = X,                      
                Z_u = X,                  
                Z_w = X,                  
                J = N_subj,
                K = N_items, 
                rt = verb$rt,                   
                subj = as.integer(subj),   
                item = as.integer(item))   

# 3. Fit the model.

E5 <- stan(file = "StanModels/maxModelmerged2.stan", 
                    data = stanDat,
                    iter = 2000, 
                    chains = 4)


#print(E5, pars=c("beta","sigma_e","sigma_u","sigma_w"), probs=c(.025,.5,.975))

# Traceplot
#traceplot(m_E5, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)

# extract the model estimates
E5_res<-stan_results(m=E5,params=c("Load","Dist","LoadxDist"))

## ----AnalysisE5RTpostcrit,cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE----

subj <- as.numeric(as.factor(verb1$subj))
#unique(verb1$subj)
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(verb1$item))
#unique(verb1$item)
N_items <- length(unique(item))

X <- unname(model.matrix(~ 1 + load + dist + int, verb1))
Z_u <- unname(model.matrix(~ 1, verb1))
Z_w <- unname(model.matrix(~ 1, verb1))

attr(X, which="assign") <- NULL

# 2. Make Stan data (list)
stanDat <- list(N = nrow(X),           
                P = ncol(X),               
                n_u = ncol(X),            
                n_w = ncol(X),              
                X = X,                      
                Z_u = X,                  
                Z_w = X,                  
                J = N_subj,
                K = N_items, 
                rt = verb1$rt,                   
                subj = as.integer(subj),   
                item = as.integer(item))   

# 3. Fit the model.

E5post <- stan(file = "StanModels/maxModelmerged2.stan", 
                    data = stanDat,
                    iter = 2000, 
                    chains = 4)


#print(E5post, pars=c("beta","sigma_e","sigma_u","sigma_w"), probs=c(.025,.5,.975))

# Traceplot
#traceplot(m_E5post, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)

# extract the model estimates
E5post_res<-stan_results(m=E5post,params=c("Load","Dist","LoadxDist"))

## ----sprmergedraw,echo=FALSE---------------------------------------------
#msprmergedcritraw<-lmer(rt~load + dist + int+(1+load + dist + int||subj)+(1+load + dist + int||item),verb)
#msprmergedpostcritraw<-lmer(rt~load + dist + int+(1+load + dist + int||subj)+(1+load + dist + int||item),verb1)

## ----E6ResponseAccdata,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE----
dat<-read.table("../data/E6ETlevykellerExp12.txt",header=T)

# subset Accuracy data
datAcc<-subset(dat,roi==1)
#head(datAcc)
#str(datAcc)
 
### quesitons only follow 50% of items: remove -2 (RESPONSE ACC) or 0 (answer) (i.e., no questions) and NAs

# subset data filler items (exclude LKrep and practice)
datFillerAcc<-subset(datAcc,condition=="f" & answer!=0)
datFillerAcc<-na.omit(datFillerAcc)

# subset data experimental items (exclude fillers and practice)
datAcc<-subset(datAcc,condition!="f" & condition!="p" & answer!=0)
datAcc<-na.omit(datAcc)

# sanity checks
#xtabs(~subject+condition,datAcc)
#length(unique(datAcc$itemid))


# "acc" column: 
datAcc <- plyr::rename(datAcc, c(RESPONSE_ACCURACY="acc"))
datFillerAcc <- plyr::rename(datFillerAcc, c(RESPONSE_ACCURACY="acc"))
# mean acc exp items
#mean(datAcc$acc) # [1] 0.6448598
# mean acc for fillers
#mean(datFillerAcc$acc) # [1] 0.898773

# compute means accuracies (condition)
#mean_acc<-round(tapply(datAcc$acc, datAcc$condition, mean),2)
#print(mean_acc)

# contrast coding: 
datAcc$load<-ifelse(datAcc$condition%in%c("a","b"),-1/2,1/2)
datAcc$dist<-ifelse(datAcc$condition%in%c("a","c"),-1/2,1/2)
datAcc$int<-ifelse(datAcc$condition%in%c("a","d"),-1/2,1/2)

E6meanacc<-round(mean(100*with(datAcc,tapply(acc,condition,mean)),na.rm=TRUE))

## ----E6AnalysisAccuracy,cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE----

subj <- as.numeric(as.factor(datAcc$subject))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(datAcc$itemid))
N_items <- length(unique(item))


X <- unname(model.matrix(~ 1 + load + dist + int, datAcc))  
Z_u <- unname(model.matrix(~ 1, datAcc)) 
Z_w <- unname(model.matrix(~ 1, datAcc))

attr(X, which="assign") <- NULL              

# 2. Make Stan data (list)
stanDatAcc <- list(N = nrow(X), 
                P = ncol(X), 
                n_u = ncol(X),
                n_w = ncol(X),
                X = X,       
                Z_u = X,    
                Z_w = X,      
                J = N_subj, 
                K = N_items,
                acc = datAcc$acc,                  
                subj = as.integer(subj),  
                item = as.integer(item))  


# 3. Fit the model.

E6Acc <- stan(file = "StanModels/responseaccmaximal.stan", 
                  data = stanDatAcc,
                  iter = 2000, 
                  chains = 4)


#print(mE6Acc, pars=c("beta","sigma_u","sigma_w"), probs=c(.025,.5,.975))

#Traceplot
#traceplot(mE6Acc, pars=c("beta","sigma_u","sigma_w"), inc_warmup=FALSE)


## ----E6RTData,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE---------


dat<-read.table("../data/E6ETlevykellerExp12.txt",header=T)

# contrast coding (same as for E5 SPR merged)
dat$load<-ifelse(dat$condition%in%c("a","b"),-1/2,1/2)
dat$dist<-ifelse(dat$condition%in%c("a","c"),-1/2,1/2)
dat$int<-ifelse(dat$condition%in%c("a","d"),-1/2,1/2)

#str(dat)
dat$roi<-factor(dat$roi)


## rois
#critical region (verb):
#a = 24 
#b = 24 
#c = 25 (merged roi 25 + 26)
#d = 25 (merged roi 25 + 26)

#precritical region (acc NP)
#a = 22 (merged 22 + 23)
#b = 22 (merged 22 + 23)
#c = 23 (merged 23 + 24)
#d = 23 (merged 23 + 24)

#postcritical region (varies: a/b = "und so/damit", c/d= "einen Verlust"/"Recht")
#a = 25 (merged 25 + 26)
#b = 25 (merged 25 + 26)
#c = 26 (merged 27 + 28) * expect for items 15,17,18,20 only region 27, now 26, (as only e.g. "Recht")
#d = 26 (merged 27 + 28) * same as above


# rois in ET (differ for conds a, b vs c,d)
region<-ifelse(dat$condition%in%c("a","b") & dat$roi==22,"npacc",
               ifelse(dat$condition%in%c("c","d") & dat$roi==23,"npacc",
                ifelse(dat$condition%in%c("a","b") & dat$roi==24,"verb",
                    ifelse(dat$condition%in%c("c","d") & dat$roi==25,"verb",   
                              ifelse(dat$condition%in%c("a","b") & dat$roi==25,"verb1",
                    ifelse(dat$condition%in%c("c","d") & dat$roi==26,"verb1","noncritical"))))))

dat$region<-region
dat<-subset(dat,region!="noncritical")
dat$region<-factor(dat$region)

# check data:
#dim(dat)
#str(dat)
#xtabs(~subject+region,dat)
#xtabs(~subject+itemid,dat)
#xtabs(~itemid+region,dat)

# subset critical region 
verb <-(subset(dat,region=="verb"))
verb1<-subset(dat,region=="verb1")

## ----AnalysisE6TFTcrit,cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE----

verbTFT <- subset(verb,verb$TFT>0)

subj <- as.numeric(as.factor(verbTFT$subject))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(verbTFT$itemid))
N_items <- length(unique(item))

X <- unname(model.matrix(~ 1 + load + dist + int, verbTFT))  
Z_u <- unname(model.matrix(~ 1, verbTFT)) 
Z_w <- unname(model.matrix(~ 1, verbTFT))

attr(X, which="assign") <- NULL   

# 2. Make Stan data.
stanDat <- list(N = nrow(X),    
                P = ncol(X), 
                n_u = ncol(X),  
                n_w = ncol(X),    
                X = X,           
                Z_u = X,         
                Z_w = X,         
                J = N_subj,
                K = N_items, 
                rt = verbTFT$TFT,                                                                      
                subj = as.integer(subj), 
                item = as.integer(item))

# 3. Fit the model.

E6_3 <- stan(file = "StanModels/maxModelmerged2.stan", 
                  data = stanDat,
                  iter = 2000, 
                 chains = 4)

# print parameters
#print(E6_3, pars=c("beta","sigma_e","sigma_u","sigma_w"), probs=c(.025,.5,.975))

# Traceplot
#traceplot(m_E6_3, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)

# extract the model estimates
E6_3_res<-stan_results(m=E6_3,params=c("Load","Dist","LoadxDist"))

## ----E6etcritraw,echo=FALSE----------------------------------------------
#mE6etcritraw<-lmer(TFT~load + dist + int+(load + dist + int||subject)+(load + dist + int||itemid),verbTFT)

## ----AnalysisE6TFTpostcrit,cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE----

verb1TFT <- subset(verb1,verb1$TFT>0)

subj <- as.numeric(as.factor(verb1TFT$subject))
#unique(verb1TFT$subject)
#unique(subj)
N_subj <- length(unique(subj))

#unique(verb1TFT$itemid)
item <- as.numeric(as.factor(verb1TFT$itemid))
#unique(item)
N_items <- length(unique(item))

X <- unname(model.matrix(~ 1 + load + dist + int, verb1TFT))  
Z_u <- unname(model.matrix(~ 1, verb1TFT)) 
Z_w <- unname(model.matrix(~ 1, verb1TFT))

attr(X, which="assign") <- NULL   

# 2. Make Stan data.
stanDat <- list(N = nrow(X),    
                P = ncol(X), 
                n_u = ncol(X),  
                n_w = ncol(X),    
                X = X,           
                Z_u = X,         
                Z_w = X,         
                J = N_subj,
                K = N_items, 
                rt = verb1TFT$TFT,                                                                     
                subj = as.integer(subj), 
                item = as.integer(item))

# 3. Fit the model.

E6_3post <- stan(file = "StanModels/maxModelmerged2.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)

# print parameters
#print(m_E6_3post, pars=c("beta","sigma_e","sigma_u","sigma_w"), probs=c(.025,.5,.975))

# Traceplot
#traceplot(m_E6_3post, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)

# extract the model estimates
E6_3post_res<-stan_results(m=E6_3post,params=c("Load","Dist","LoadxDist"))

## ----E6etpostcritraw,echo=FALSE------------------------------------------
#mE6etpostcritraw<-lmer(TFT~load + dist + int+(load + dist + int||subject)+(load + dist + int||itemid),verb1TFT)

## ----E7AccData,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE--------

dat<-read.table("../data/data_LK13rep100subj.txt",header=T)

# subset Accuracy data
datAcc<-subset(dat,roi==1)
#head(datAcc)
#str(datAcc)
 
### quesitons only follow 50% of items: remove -2 (RESPONSE ACC) or 0 (answer) (i.e., no questions) and NAs

# subset data filler items (exclude LKrep and practice)
datFillerAcc<-subset(datAcc,condition=="f" & answer!=0)
datFillerAcc<-na.omit(datFillerAcc)

# subset data experimental items (exclude fillers and practice)
datAcc<-subset(datAcc,condition!="f" & condition!="p" & answer!=0)
datAcc<-na.omit(datAcc)

# sanity checks
#xtabs(~subject+condition,datAcc)
#length(unique(datAcc$itemid))


# "acc" column: 
datAcc <- plyr::rename(datAcc, c(RESPONSE_ACCURACY="acc"))
datFillerAcc <- plyr::rename(datFillerAcc, c(RESPONSE_ACCURACY="acc"))
# mean acc exp items
#mean(datAcc$acc) # [1] 0.6448598
# mean acc for fillers
#mean(datFillerAcc$acc) # [1] 0.898773

# compute means accuracies (condition)
#mean_acc<-round(tapply(datAcc$acc, datAcc$condition, mean),2)
#print(mean_acc)

# contrast coding: 
datAcc$load<-ifelse(datAcc$condition%in%c("a","b"),-1/2,1/2)
datAcc$dist<-ifelse(datAcc$condition%in%c("a","c"),-1/2,1/2)
datAcc$int<-ifelse(datAcc$condition%in%c("a","d"),-1/2,1/2)

E7meanacc<-round(mean(100*with(datAcc,tapply(acc,condition,mean)),na.rm=TRUE))


## ----E7AnalysisAcc,cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE-----

subj <- as.numeric(as.factor(datAcc$subject))
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(datAcc$itemid))
N_items <- length(unique(item))


X <- unname(model.matrix(~ 1 + load + dist + int, datAcc))  
Z_u <- unname(model.matrix(~ 1, datAcc)) 
Z_w <- unname(model.matrix(~ 1, datAcc))

attr(X, which="assign") <- NULL              

# 2. Make Stan data (list)
stanDatAcc <- list(N = nrow(X), 
                P = ncol(X), 
                n_u = ncol(X),
                n_w = ncol(X),
                X = X,       
                Z_u = X,    
                Z_w = X,      
                J = N_subj, 
                K = N_items,
                acc = datAcc$acc,                  
                subj = as.integer(subj),  
                item = as.integer(item))  


# 3. Fit the model.

mE7Acc <- stan(file = "StanModels/responseaccmaximal.stan", 
                  data = stanDatAcc,
                  iter = 2000, 
                  chains = 4)


#print(mE7Acc, pars=c("beta","sigma_u","sigma_w"), probs=c(.025,.5,.975))

#Traceplot
#traceplot(mE7Acc, pars=c("beta","sigma_u","sigma_w"), inc_warmup=FALSE)


## ----E7rtData,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE---------

dat<-read.table("../data/data_LK13rep100subj.txt",header=T)

dat<-subset(dat,condition!="f" & condition!="p")

# contrast coding (same as for E5 SPR merged)
dat$load<-ifelse(dat$condition%in%c("a","b"),-1/2,1/2)
dat$dist<-ifelse(dat$condition%in%c("a","c"),-1/2,1/2)
dat$int<-ifelse(dat$condition%in%c("a","d"),-1/2,1/2)

## nested contrasts
dat$lodist<-ifelse(dat$condition=="a",-1/2,
                   ifelse(dat$condition=="b",1/2,0))
dat$hidist<-ifelse(dat$condition=="c",-1/2,
                   ifelse(dat$condition=="d",1/2,0))

#xtabs(~condition+hidist,dat)


#str(dat)
#head(dat)
#dat$roi<-factor(dat$roi)

# rois in ET (differ for conds a, b vs c,d)
region<-ifelse(dat$condition%in%c("a","b") & dat$roi==22,"npacc",
               ifelse(dat$condition%in%c("c","d") & dat$roi==23,"npacc",
                ifelse(dat$condition%in%c("a","b") & dat$roi==24,"verb",
                    ifelse(dat$condition%in%c("c","d") & dat$roi==25,"verb",   
                              ifelse(dat$condition%in%c("a","b") & dat$roi==25,"verb1",
                    ifelse(dat$condition%in%c("c","d") & dat$roi==26,"verb1","noncritical"))))))

#length(region)  

dat$region<-region
dat<-subset(dat,region!="noncritical")
dat$region<-factor(dat$region)


# check data:
#xtabs(~subject+region,dat)    # subj 41 (bad cali) saw rois x 19 (instead of 24)
#xtabs(~subject+itemid,dat)
#xtabs(~itemid+region,dat)

#test1<-subset(dat,dat$itemid==6)

#xtabs(~subject+itemid,test1)
#xtabs(~subject+condition,test1)
#xtabs(~itemid+condition,test1)

#subset critical region 
verb <-subset(dat,region=="verb")

## quick sanity check to see if the conventional analysis would have found anything that TRT could not tell you, answer is no:
#summary(mFPRT<-lmer(log(FPRT)~load+dist+int+(1|subject)+(1|itemid),subset(verb,FPRT>0)))
#summary(mRPD<-lmer(log(RPD)~load+dist+int+(1|subject)+(1|itemid),subset(verb,RPD>0)))
#summary(mRRT<-lmer(log(RRT)~load+dist+int+(1|subject)+(1|itemid),subset(verb,RRT>0)))
## summary(mTFT<-lmer(log(TFT)~load+dist+int+(1|subject)+(1|itemid),subset(verb,TFT>0)))

#summary(mTFT<-lmer(log(TFT)~load+lodist+hidist+(1|subject)+(1|itemid),subset(verb,TFT>0)))

#summary(m)
#subset postcritical region 
verb1 <-subset(dat,region=="verb1")

#xtabs(~subject+condition,verb)
#xtabs(~subject+condition,verb1)

## ----AnalysisE7TFTcrit,include=FALSE,results='hide',cache=TRUE,echo=FALSE,warning=FALSE----

verbTFT <- subset(verb,verb$TFT>0)
#head(verb)

subj <- as.numeric(as.factor(verbTFT$subject))
#unique(verbTFT$subject)
#unique(subj)
N_subj <- length(unique(subj))

#unique(verbTFT$itemid)
item <- as.numeric(as.factor(verbTFT$itemid))
#unique(item)
N_items <- length(unique(item))

X <- unname(model.matrix(~ 1 + load + dist + int, verbTFT))  
attr(X, which="assign") <- NULL   

# 2. Make Stan data.
stanDat <- list(N = nrow(X),    
                P = ncol(X), 
                n_u = ncol(X),  
                n_w = ncol(X),    
                X = X,           
                Z_u = X,         
                Z_w = X,         
                J = N_subj,
                K = N_items, 
                rt = verbTFT$TFT,                                                                      
                subj = as.integer(subj), 
                item = as.integer(item))


# 3. Fit the model.

E7_3 <- stan(file = "StanModels/maxModelmerged2.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)

# print parameters
#print(E7_3, pars=c("beta","sigma_e","sigma_u","sigma_w"), digits=4, probs=c(.025,.5,.975))

# check
#m13 <- lmer(log(TFT)~load+dist+int+(1+load+dist+int||subject)+(1+load+dist+int||itemid),verbTFT)
#summary(m13)

# Traceplot
#traceplot(m_E7_3, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)

# extract the model estimates
E7_3_res<-stan_results(m=E7_3,params=c("Load","Dist","LoadxDist"))

#round(with(verbTFT,tapply(TFT,condition,mean)))

## ----E7etcritraw,echo=FALSE----------------------------------------------
#mE7etcritraw<-lmer(TFT~load + dist + int+(load + dist + int||subject)+(load + dist + int||itemid),verbTFT)

## ----E7etpostcritraw,echo=FALSE------------------------------------------
#mE7etpostcritraw<-lmer(TFT~load + dist + int+(load + dist + int||subject)+(load + dist + int||itemid),verb1TFT)

## ----AnalysisE7TFTcritsmallsample,echo=FALSE,eval=FALSE------------------
## nsim<-50
## estimates<-matrix(rep(NA,nsim*3*3),ncol=9)
## subjID<-unique(verbTFT$subj)
## 
## for(i in 1:nsim){
## small_sample<-sample(subjID,28)
## dsmall<-subset(verbTFT,subject%in%c(small_sample))
## stanDatsmall<-createStanDat(d=dsmall,rt=dsmall$TFT,
##               form=as.formula("~1+load+dist+int"))
## #str(standatsmall)
## E7_3small <- stan(file = "StanModels/maxModelmerged2.stan",
##                   data = stanDatsmall,
##                   iter = 2000,
##                   chains = 4)
## estimates[i,]<-as.vector(stan_results(m=E7_3small,params=c("Load","Dist","LoadxDist")))
## }
## save(estimates,file="../data/smallsamplesestimates.Rda")

## ----AnalysisE7TFTcritnested,include=FALSE,results='hide',cache=TRUE,echo=FALSE,warning=FALSE----
## nested:
X <- unname(model.matrix(~ 1 + load + lodist + hidist, verbTFT))  
attr(X, which="assign") <- NULL   
#head(X)
Z_u <- Z_w <- X
#  unname(model.matrix(~ 1 , verbTFT)) 
#attr(Z_u, which="assign") <- NULL   
#attr(Z_w, which="assign") <- NULL   



stanDat <- list(N = nrow(X),    
                P = ncol(X), 
                n_u = ncol(Z_u),  
                n_w = ncol(Z_w),    
                X = X,           
                Z_u = Z_u,         
                Z_w = Z_w,         
                J = N_subj,
                K = N_items, 
                rt = verbTFT$TFT,               
                subj = as.integer(subj), 
                item = as.integer(item))

#str(stanDat)

E7_3nested <- stan(file = "StanModels/maxModelmerged2nested.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)

E7_3nested_res<-stan_results(m=E7_3nested,params=c("Load","LoDist","HiDist"))

## ----AnalysisE7TFTpostcrit,cache=TRUE,echo=FALSE,warning=FALSE,message=FALSE----
verb1TFT <- subset(verb1,verb1$TFT>0)
#head(verb)
#round(with(verb1TFT,tapply(TFT,condition,mean)))


subj <- as.numeric(as.factor(verb1TFT$subject))
#unique(verb1TFT$subject)
#unique(subj)
N_subj <- length(unique(subj))

item <- as.numeric(as.factor(verb1TFT$itemid))
#unique(item)
N_items <- length(unique(item))

X <- unname(model.matrix(~ 1 + load + dist + int, verb1TFT))  

attr(X, which="assign") <- NULL   

# 2. Make Stan data.
stanDat <- list(N = nrow(X),    
                P = ncol(X), 
                n_u = ncol(X),  
                n_w = ncol(X),    
                X = X,           
                Z_u = X,         
                Z_w = X,         
                J = N_subj,
                K = N_items, 
                rt = verb1TFT$TFT,               
                subj = as.integer(subj), 
                item = as.integer(item))


# 3. Fit the model.
#library(rstan)
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

E7_3post <- stan(file = "StanModels/maxModelmerged2.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)

# print parameters
#print(E7_3post, pars=c("beta","sigma_e","sigma_u","sigma_w"), digits=4, probs=c(.025,.5,.975))

# check
#m13 <- lmer(log(TFT)~load+dist+int+(1+load+dist+int||subject)+(1+load+dist+int||itemid),verb1TFT)
#m13raw <- lmer(TFT~load+dist+int+(1+load+dist+int||subject)+(1+load+dist+int||itemid),verb1TFT)
#summary(m13raw)

# Traceplot
#traceplot(m_E7_3post, pars=c("beta","sigma_e","sigma_u","sigma_w"), inc_warmup=FALSE)

# extract the model estimates
E7_3post_res<-stan_results(m=E7_3post,params=c("Load","Dist","LoadxDist"))


## ----AnalysisE7TFTpostcritnested,include=FALSE,results='hide',cache=TRUE,echo=FALSE,warning=FALSE----
## nested:
X <- unname(model.matrix(~ 1 + load + lodist + hidist, verb1TFT))  
attr(X, which="assign") <- NULL   

Z_u <- Z_w <- X

stanDat <- list(N = nrow(X),    
                P = ncol(X), 
                n_u = ncol(Z_u),  
                n_w = ncol(Z_w),    
                X = X,           
                Z_u = Z_u,         
                Z_w = Z_w,         
                J = N_subj,
                K = N_items, 
                rt = verb1TFT$TFT,                                subj = as.integer(subj), 
                item = as.integer(item))

E7_3postnested <- stan(file = "StanModels/maxModelmerged2nested.stan", 
                  data = stanDat,
                  iter = 2000, 
                  chains = 4)

E7_3postnested_res<-stan_results(m=E7_3postnested,params=c("Load","LoDist","HiDist"))

## ----figurese3,echo=FALSE,include=TRUE,warning=FALSE,message=FALSE-------
cond<-factor(c("Load","Dist","LoadxDist"),levels=c("Load","Dist","LoadxDist"))
cond<-rep(cond,4)
expt<-factor(c(rep("Original LK data",3),
        rep("E5, SPR",3),rep("E6, ET",3),rep("E7, ET (n=100)",3)),levels=c("Original LK data","E5, SPR","E6, ET","E7, ET (n=100)"))
  

rownames(LKmerged_res)<-NULL
rownames(E5_res)<-NULL
rownames(E6_3_res)<-NULL
rownames(E7_3_res)<-NULL

LKmergedall<-data.frame(rbind(LKmerged_res,E5_res,E6_3_res,E7_3_res))
LKmergedall<-data.frame(cond=cond,expt=expt,LKmergedall)


p1<-plotresults(LKmergedall,maintitle="Load vs. Distance (critical region)",,cols=c("darkgray", "black","black","black"),legendposition=c(0.8,0.75),lowery=-150,uppery=400)

rownames(LKmergedpost_res)<-NULL
rownames(E5post_res)<-NULL
rownames(E6_3post_res)<-NULL
rownames(E7_3post_res)<-NULL

LKmergedpostall<-data.frame(rbind(LKmergedpost_res,E5post_res,
                                  E6_3post_res,E7_3post_res))
LKmergedpostall<-data.frame(cond=cond,expt=expt,LKmergedpostall)
p2<-plotresults(LKmergedpostall,maintitle="Load vs. Distance (post-critical region)",cols=c("darkgray", "black","black","black"),removelegend=TRUE,lowery=-150,uppery=400)
multiplot(p1,p2,cols=1)

## ----loadsmallsampleestimates,echo=FALSE---------------------------------
load(file="../data/smallsamplesestimates.Rda")

## ----E7smallsampleplot,echo=FALSE,include=TRUE,fig.width=8,fig.height=6----
est_Load<-estimates[,c(1,4,7)]
#dim(est_Load)
est_Dist<-estimates[,c(1,4,7)+1]
est_LoadxDist<-estimates[,c(1,4,7)+2]

est<-data.frame(rbind(est_Load,est_Dist,est_LoadxDist))
colnames(est)<-c("mean","lower","upper")
est$significance<-ifelse(sign(est$lower)==sign(est$upper),"p<0.05","p>0.05")
est$predictor<-factor(rep(c("Load","Dist","LoadxDist"),each=50))
est$sample<-rep(1:50,3)

estLoad<-subset(est,predictor=="Load")
estDist<-subset(est,predictor=="Dist")
estLoadxDist<-subset(est,predictor=="LoadxDist")

#summary(abs(estDist$mean/E7_3_res[2,1]))

pd<-position_dodge(0.6)

plot_Load<-ggplot(estLoad, aes(x=sample, 
                             y=mean)) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=.25, size=.5, position=pd) +
  labs(title="Effect of Load: \n Estimates from repeated samples (n=28)") +
  xlab("sample")+
  ylab("estimate (ms)")+
  geom_hline(yintercept=E7_3_res[1,1])+
  geom_hline(yintercept=E7_3_res[1,2],linetype=2)+
  geom_hline(yintercept=E7_3_res[1,3],linetype=2)+
  geom_point(position=pd, size=2)+
  theme_bw()+magnifytext()+coord_flip()

estDist<-estDist[order(estDist$mean), ]
## relabel samples in increasing order:
estDist$sample<-1:50

plot_Dist<-ggplot(estDist, aes(x=sample, 
                             y=mean,shape=significance,ymin=lower, ymax=upper)) +
  #geom_errorbar(aes(ymin=lower, ymax=upper),
  #              width=.25, size=.5, position=pd) +
  geom_pointrange()+
  geom_point(size=2.5)+
  scale_shape_manual(values=c(1, 17))+
  labs(title="Effect of Distance: \n Estimates from repeated samples (n=28)") +
  xlab("Sample id")+
  ylab("Estimates (ms)")+
  geom_hline(yintercept=E7_3_res[2,1])+
  xlim(1,50)+
  geom_hline(yintercept=E7_3_res[2,2],linetype=2)+
  geom_hline(yintercept=E7_3_res[2,3],linetype=2)+
  geom_point(position=pd, size=2)+
  theme_bw()+magnifytext()+geom_hline(yintercept=0, linetype="dotted")
#+coord_flip()

plot_LoadxDist<-ggplot(estLoadxDist, aes(x=sample, 
                             y=mean)) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=.25, size=.5, position=pd) +
  labs(title="Effect of Load x Distance: \n Estimates from repeated samples (n=28)") +
  xlab("")+
  ylab("Estimates (ms)")+
  geom_hline(yintercept=E7_3_res[3,1])+
  geom_hline(yintercept=E7_3_res[3,2],linetype=2)+
  geom_hline(yintercept=E7_3_res[3,3],linetype=2)+
  geom_point(position=pd, size=2)+
  theme_bw()+magnifytext()+coord_flip()

plot_Dist

## ----echo=FALSE,eval=FALSE-----------------------------------------------
## test <- rbind(mtcars[1:10, 1:5],
##               colSums(mtcars[1:10, 1:5]))
## rownames(test)[11] <- "Sum"
## 
## test[11, ] <- paste0("BOLD", test[11, ])
## 
## bold.somerows <-
##         function(x) gsub('BOLD(.*)',paste('\\\\textbf{\\1','}'),x)
## 
## test.xt <- xtable(test, label="table", caption='test')
## #align(test.xt) <- "|l|l|l|r|r|r|"#
## #print(test.xt, type="latex",tabular.environment='tabular', include.rownames = FALSE, floating=TRUE, sanitize.text.function = bold.somerows)

## ----agrmtattrn,echo=FALSE,warning=FALSE,message=FALSE,eval=FALSE--------
## 
## ## Dillon et al 2013 Expt 1
## DillonE1<-read.table("../data/DillonE1.txt",header=T)
## ## Lago et al 2015 data (all expts):
## Lago<-read.csv("../data/Lago.csv",header=T)
## ##Wagers et al 2009 data (all expts):
## load("../data/Wagers.Rdata")
## load("../data/Tucker.RData")
## 
## ## Dillon E1:
## DillonE1$cond<-factor(DillonE1$cond)
## DillonE1Mism<-subset(DillonE1,fixationtype=="tt" & cond%in%c(3,4) & value!="NA")
## DillonE1Mism$cond<-factor(DillonE1Mism$cond)
## DillonE1Mism$int<-ifelse(DillonE1Mism$cond==3,"low","high")
## DillonE1Mism$x<-ifelse(DillonE1Mism$cond==3,-1,1)
## dillonE1<-DillonE1Mism[,c(1,3,4,14,15)]
## dillonE1$expt<-factor("dillonE1")
## colnames(dillonE1)[3]<-"rt"
## 
## nsubj_dillonE1<-length(unique(dillonE1$subj))
## 
## ##Lago:
## dat<-Lago
## ## critical region: not used because published paper found
## ## significant effects in postcrit region only
## e1<-subset(dat,Experiment=="Experiment1" & Region=="06v1")
## e2<-subset(dat,Experiment=="Experiment2" & Region=="06aux")
## e3a<-subset(dat,Experiment=="Experiment3A" & Region=="06aux")
## e3b<-subset(dat,Experiment=="Experiment3B" & Region=="aux")
## 
## nsubj_lagoe1<-length(unique(e1$Subject))
## nsubj_lagoe2<-length(unique(e2$Subject))
## nsubj_lagoe3a<-length(unique(e3a$Subject))
## nsubj_lagoe3b<-length(unique(e3b$Subject))
## 
## 
## ## postcritical region:
## poste1<-subset(dat,Experiment=="Experiment1" & Region=="07prep")
## poste2<-subset(dat,Experiment=="Experiment2" & Region=="07adv")
## poste3a<-subset(dat,Experiment=="Experiment3A" & Region=="07a")
## poste3b<-subset(dat,Experiment=="Experiment3B" & Region=="a")
## 
## ##e1: a,b
## #-(a) Ungram , singular attractor (interference condition)
## #La *nota* que la chica escribieron en la clase alegró a su amiga
## #The note that the girl wrotepl during class cheered her friend up
## #-(b) Ungram , plural attractor (baseline condition)
## #Las *notas* que la chica escribieron en la clase alegraron a su amiga
## #The notes that the girl wrotepl during class cheered her friend up
## poste1<-subset(poste1,Condition%in%c("a","b"))
## poste1$Condition<-factor(poste1$Condition)
## poste1$x<-ifelse(poste1$Condition=="a",-1,1)
## poste1$int<-ifelse(poste1$Condition=="a","low","high")
## poste1<-poste1[,c(1,3,8,15,14)]
## poste1$expt<-factor("lagoE1")
## lagoE1<-poste1
## colnames(lagoE1)<-c("subj","item","rt","int","x","expt")
## 
## ## e2: c,d
## poste2<-subset(poste2,Condition%in%c("c","d"))
## poste2$Condition<-factor(poste2$Condition)
## poste2$x<-ifelse(poste2$Condition=="c",-1,1)
## poste2$int<-ifelse(poste2$Condition=="c","low","high")
## #head(poste2)
## poste2<-poste2[,c(1,3,8,15,14)]
## poste2$expt<-factor("lagoE2")
## lagoE2<-poste2
## colnames(lagoE2)<-c("subj","item","rt","int","x","expt")
## 
## ## e3a: e,f
## poste3a<-subset(poste3a,Condition%in%c("e","f"))
## poste3a$Condition<-factor(poste3a$Condition)
## #-(e) Ungram, singular attractor (interference condition)
## #La *nota* que la chica van a escribir en la clase alegrará a su amiga
## #The note that the girl are going to write during class will cheer her friend up
## #-(f) Ungram, plural attractor (baseline condition)
## #Las *notas* que la chica van a escribir en la clase alegrarán a su amiga
## #The notes that the girl are going to write during class will cheer her friend up
## #boxplot(RT~Condition,poste3a)
## poste3a$x<-ifelse(poste3a$Condition=="e",-1,1)
## poste3a$int<-ifelse(poste3a$Condition=="e","low","high")
## poste3a<-poste3a[,c(1,3,8,15,14)]
## poste3a$expt<-factor("lagoE3a")
## lagoE3a<-poste3a
## colnames(lagoE3a)<-c("subj","item","rt","int","x","expt")
## 
## ## e3b: e,f
## poste3b<-subset(poste3b,Condition%in%c("e","f"))
## poste3b$Condition<-factor(poste3b$Condition)
## #-(e) Ungram, singular attractor (baseline condition)
## #The player that the coach were always praising very enthusiastically decided to     leave the team
## #-(f) Ungram, plural attractor (interference condition)
## #The players that the coach were always praising very enthusiastically decided to     leave the team
## poste3b$x<-ifelse(poste3b$Condition=="e",-1,1)
## poste3b$int<-ifelse(poste3b$Condition=="e","low","high")
## poste3b<-poste3b[,c(1,3,8,15,14)]
## poste3b$expt<-factor("lagoE3b")
## lagoE3b<-poste3b
## colnames(lagoE3b)<-c("subj","item","rt","int","x","expt")
## 
## ## Wagers:
## E2postcrit<-subset(Experiment2,Region==7)
## nsubj_wagerse2<-length(unique(E2postcrit$Subj))
## #E2$intr.au<-ifelse(E2$rchead=="pl" & E2$gramm=="ungram",1/2,
## #                   ifelse(E2$rchead=="sg" & E2$gramm=="ungram",-1/2,
## #                          0))
## ## d (sing),h (plu)
## #unique(subset(E2postcrit,gramm=="ungram")$Condition)
## E2postcrit<-subset(E2postcrit,Condition%in%c("d","h"))
## E2postcrit$Condition<-factor(E2postcrit$Condition)
## E2postcrit$x<-ifelse(E2postcrit$Condition=="d",-1,1)
## E2postcrit$int<-ifelse(E2postcrit$Condition=="d","low","high")
## #colnames(E2postcrit)
## E2postcrit<-E2postcrit[,c(4,3,8,13,12)]
## E2postcrit$expt<-factor("wagersE2")
## wagersE2<-E2postcrit
## colnames(wagersE2)<-c("subj","item","rt","int","x","expt")
## 
## ## E3
## E3postcrit<-subset(Experiment3,Region==7)
## nsubj_wagerse3<-length(unique(E3postcrit$Subj))
## #E3crit$intr.au.pl<-ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="sg" &
## #                            E3crit$rchead=="pl",1/2,
## #                         ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="sg" &
## #                                   E3crit$rchead=="sg",-1/2,0))
## 
## #E3crit$intr.au.sg<-ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="pl" &
## #                            E3crit$rchead=="sg",1/2,
## #                          ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="pl" &
## #                                   E3crit$rchead=="pl",-1/2,0))
## 
## E3postcrit_pl<-subset(E3postcrit,gramm=="ungram" & rcsubj=="sg")
## E3postcrit_pl$Condition<-factor(E3postcrit_pl$Condition)
## E3postcrit_sg<-subset(E3postcrit,gramm=="ungram" & rcsubj=="pl")
## E3postcrit_sg$Condition<-factor(E3postcrit_sg$Condition)
## 
## #unique(E3postcrit_pl$Condition) ## b,f
## #unique(E3postcrit_sg$Condition) ## c,g
## 
## #head(subset(E3postcrit_sg,rchead=="sg"))
## #head(E3postcrit_pl)
## 
## ## plural:
## E3postcrit_pl$x<-ifelse(E3postcrit_pl$Condition=="b",-1,1)
## E3postcrit_pl$int<-ifelse(E3postcrit_pl$Condition=="b","low","high")
## E3postcrit_pl<-E3postcrit_pl[,c(4,3,8,15,14)]
## E3postcrit_pl$expt<-factor("wagersE3pl")
## colnames(E3postcrit_pl)<-c("subj","item","rt","int","x","expt")
## wagersE3pl<-E3postcrit_pl
## 
## ## singular:
## E3postcrit_sg$x<-ifelse(E3postcrit_sg$Condition=="c",-1,1)
## E3postcrit_sg$int<-ifelse(E3postcrit_sg$Condition=="c","low","high")
## E3postcrit_sg<-E3postcrit_sg[,c(4,3,8,15,14)]
## E3postcrit_sg$expt<-factor("wagersE3sg")
## colnames(E3postcrit_sg)<-c("subj","item","rt","int","x","expt")
## wagersE3sg<-E3postcrit_sg
## 
## ## E4
## E4postcrit<-subset(Experiment4,Region==8) ##
## nsubj_wagerse4<-length(unique(E4postcrit$Subj))
## #head(subset(Experiment4,Condition=="c"),n=10)
## #postcritical region
## #E4postcrit$intr.au<-ifelse(E4postcrit$gramm=="ungram" & E4postcrit$match=="match",-1/2,
## #                           ifelse(E4postcrit$gramm=="ungram" & E4postcrit$match=="mismatch",1/2,0))
## E4postcrit<-subset(E4postcrit,gramm=="ungram")
## E4postcrit$Condition<-factor(E4postcrit$Condition)
## E4postcrit$x<-ifelse(E4postcrit$Condition=="c",-1,1)
## E4postcrit$int<-ifelse(E4postcrit$Condition=="c","low","high")
## E4postcrit<-E4postcrit[,c(4,3,8,13,12)]
## E4postcrit$expt<-factor("wagersE4")
## colnames(E4postcrit)<-c("subj","item","rt","int","x","expt")
## wagersE4<-E4postcrit
## 
## # E5
## E5postcrit<-subset(Experiment5,Region==8) ##postcritical region
## nsubj_wagerse5<-length(unique(E5postcrit$Subj))
## E5postcrit<-subset(E5postcrit,gramm=="ungram")
## E5postcrit$Condition<-factor(E5postcrit$Condition)
## ## c,d
## E5postcrit$x<-ifelse(E5postcrit$Condition=="c",-1,1)
## E5postcrit$int<-ifelse(E5postcrit$Condition=="c","low","high")
## E5postcrit<-E5postcrit[,c(4,3,8,13,12)]
## colnames(E5postcrit)<-c("subj","item","rt","int","x")
## E5postcrit$expt<-factor("wagersE5")
## wagersE5<-E5postcrit
## 
## #head(wagersE5)
## 
## dat<-rbind(dillonE1,wagersE2,
##            lagoE1,lagoE2,
##            lagoE3a,lagoE3b,
##            wagersE2,
##            wagersE3pl,wagersE3sg,
##            wagersE4,wagersE5)
## dat$subj<-factor(paste(dat$expt,dat$subj,sep=""))
## dat$item<-factor(paste(dat$expt,dat$item,sep=""))
## with(dat,tapply(subj,expt,function(x)length(unique(x))))
## 

## ----agrmtattrn2,echo=FALSE,warning=FALSE,message=FALSE,eval=FALSE-------
## ## Dillon E1:
## stanDat<-createStanDat(d=dillonE1,
##                        rt=dillonE1$rt,
##                        form=as.formula("~ 1 + x"))
## #str(stanDat)
## DillonE1 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##              data = stanDat,
##              iter = 2000,
##              chains = 4)
## ## Int is Interference:
## pars<-c("Int","beta[2]","sigma_u[1]","sigma_u[2]","sigma_w[1]","sigma_w[2]","sigma_e")
## DillonE1_res<-stan_results(DillonE1,params=pars[1])
## 
## stanDat<-createStanDat(d=lagoE1,
##                        rt=lagoE1$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE1 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4)
## LagoE1_res<-stan_results(LagoE1,params=pars[1])
## 
## stanDat<-createStanDat(d=lagoE2,
##                        rt=lagoE2$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE2 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                data = stanDat,
##                iter = 2000,
##                chains = 4)
## LagoE2_res<-stan_results(LagoE2,params=pars[1])
## 
## stanDat<-createStanDat(d=lagoE3a,
##                           rt=lagoE3a$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE3a <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                data = stanDat,
##                iter = 2000,
##                chains = 4)
## 
## LagoE3a_res<-stan_results(LagoE3a,params=pars[1])
## 
## stanDat<-createStanDat(d=lagoE3b,
##                           rt=lagoE3b$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE3b <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                 data = stanDat,
##                 iter = 2000,
##                 chains = 4)
## 
## LagoE3b_res<-stan_results(LagoE3b,params=pars[1])
## 
## stanDat<-createStanDat(d=wagersE2,
##                           rt=wagersE2$rt,
##                        form=as.formula("~ 1 + x"))
## 
## WagersE2 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                 data = stanDat,
##                 iter = 2000,
##                 chains = 4)
## WagersE2_res<-stan_results(WagersE2,params=pars[1])
## 
## stanDat<-createStanDat(d=wagersE3pl,
##                           rt=wagersE3pl$rt,
##                        form=as.formula("~ 1 + x"))
## 
## WagersE3pl <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4)
## WagersE3pl_res<-stan_results(WagersE3pl,params=pars[1])
## 
## stanDat<-createStanDat(d=wagersE3sg,
##                           rt=wagersE3sg$rt,
##                           form=as.formula("~ 1 + x"))
## 
## WagersE3sg <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                    data = stanDat,
##                    iter = 2000,
##                    chains = 4)
## WagersE3sg_res<-stan_results(WagersE3sg,params=pars[1])
## 
## stanDat<-createStanDat(d=wagersE4,
##                           rt=wagersE4$rt,
##                           form=as.formula("~ 1 + x"))
## 
## WagersE4 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                    data = stanDat,
##                    iter = 2000,
##                    chains = 4)
## WagersE4_res<-stan_results(WagersE4,params=pars[1])
## 
## stanDat<-createStanDat(d=wagersE5,
##                        rt=wagersE5$rt,
##                        form=as.formula("~ 1 + x"))
## 
## WagersE5 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4)
## WagersE5_res<-stan_results(WagersE5,params=pars[1])
## 
## ## cunnings and sturt 2018
## CS18E1<-c(-22,-42,-4)
## 
## CS18E2<-c(-19,-40,1)
## 
## posteriors<-data.frame(expt=factor(1:12),rbind(DillonE1_res,LagoE1_res,LagoE2_res,LagoE3a_res,LagoE3b_res,WagersE2_res,WagersE3pl_res,WagersE3sg_res,WagersE4_res,WagersE5_res,CS18E1,CS18E2))
## 
## ## order expts to match posteriors:
## dat$expt<-factor(dat$expt,levels=c("dillonE1","lagoE1","lagoE2","lagoE3a","lagoE3b","wagersE2","wagersE3pl","wagersE3sg","wagersE4","wagersE5","CS18E1","CS18E2"))
## #levels(dat$expt)
## 
## n_subj<-as.vector(with(dat,tapply(subj,expt,function(x)length(unique(x)))))
## n_item<-as.vector(with(dat,tapply(item,expt,function(x)length(unique(x)))))
## n_subj[11:12]<-c(48,48)
## n_item[11:12]<-c(32,32)
## 
## posteriors$n_subj<-n_subj
## posteriors$n_item<-n_item
## ## width of credible interval
## posteriors$width<-posteriors$upper-posteriors$lower
## 
## ## reorder by mean magnitude:
## posteriors<-posteriors[with(posteriors,order(mean)),]
## 
## ## reorder expt by mean's magnitude:
## posteriors$expt <- factor(posteriors$expt,levels=posteriors$expt[order(posteriors$mean)])
## save(posteriors,file="../data/posteriorsTargetMismatch.Rda")
## #xtable(posteriors)

## ----echo=FALSE----------------------------------------------------------
load(file="../data/posteriorsTargetMismatch.Rda")
posteriors<-rbind(posteriors,
      c(NA,-26,-57,-10,NA,NA,NA))


## ----powerz,echo=FALSE,fig.width=4,fig.height=4,include=TRUE-------------
zalpha<-qnorm(0.95)
zp<-seq(0,4,by=0.01)
pval<-1-pnorm(zp)
plot(zp,(1-pnorm(zalpha-zp)),type="l",
     xlab="Unknown z-score of true effect",
     ylab="Power",
     ylim=c(0,1),axes=FALSE)
axis(1, pos=0)
axis(2, pos=0)
#plot(pval,(1-pnorm(zalpha-zp)),type="l",
#     xlab="p-value",ylab="Power")

## ----echo=FALSE----------------------------------------------------------
n<-36
d<-1/10
s<-1
pow<-power.t.test(delta=d,sd=s,n=n,alternative = "one.sided",type="one.sample",strict=TRUE)$power
#pow<-1-pnorm(zalpha-.1/(1/sqrt(36)))

## ----fakedatasim,echo=FALSE,include=TRUE---------------------------------
library(MASS)
## assumes that no. of subjects and no. of items is divisible by 4.
gen_fake_lnorm2x2<-function(nitem=24,
                         nsubj=28,
                         beta=NULL,
                         Sigma_u=NULL, # subject vcov matrix
                         Sigma_w=NULL, # item vcov matrix
                         sigma_e=NULL){
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
  
  ## subject random effects:
  u<-mvrnorm(n=length(unique(fakedat$subj)),
             mu=c(0,0,0,0),Sigma=Sigma_u)
  
  ## item random effects
  w<-mvrnorm(n=length(unique(fakedat$item)),
             mu=c(0,0,0,0),Sigma=Sigma_w)

  ## generate data row by row:  
  N<-dim(fakedat)[1]
  rt<-rep(NA,N)
  for(i in 1:N){
    rt[i] <- rlnorm(1,beta[1] + 
                      u[fakedat[i,]$subj,1] +
                      w[fakedat[i,]$item,1] + 
                      (beta[2]+u[fakedat[i,]$subj,2]+
                         w[fakedat[i,]$item,2])*fakedat$c1[i]+
                      (beta[3]+u[fakedat[i,]$subj,3]+
                         w[fakedat[i,]$item,3])*fakedat$c2[i]+
                      (beta[4]+u[fakedat[i,]$subj,4]+
                         w[fakedat[i,]$item,4])*fakedat$c3[i],
                   sigma_e) 
  }   
  fakedat$rt<-rt
  fakedat$subj<-factor(fakedat$subj)
  fakedat$item<-factor(fakedat$item)
  fakedat}

## ----reanalyzeE1,warning=FALSE,echo=FALSE,include=TRUE-------------------
reading_time <- read.table('../data/OrigLevyKellerData/prediction_experiment_data/experiment1/lmr/results/exp1_tt_r.res', header=TRUE)

condition<-ifelse(reading_time$dat=="sub" & reading_time$adj=="sub","a",
                  ifelse(reading_time$dat=="sub" & reading_time$adj=="main","b",
                         ifelse(reading_time$dat=="main" & reading_time$adj=="sub","c", 
                                ifelse(reading_time$dat=="main" & reading_time$adj=="main","d","NA"))))

reading_time$condition<-factor(condition)

# contrast coding: 
reading_time$dat<-ifelse(reading_time$condition%in%c("a","b"),1/2,-1/2)
reading_time$adj<-ifelse(reading_time$condition%in%c("b","d"),-1/2,1/2)
reading_time$int<-ifelse(reading_time$condition%in%c("b","c"),-1/2,1/2)

                 ## ME DAT ## ME PP-ADJ ## INT
# a DAT-SC; PP-SC    0.5         0.5       0.5    
# b DAT-SC; PP-MC    0.5        -0.5      -0.5
# c DAT-MC; PP-SC   -0.5         0.5      -0.5
# d DAT-MC; PP-MC   -0.5        -0.5       0.5

# remove zeros
reading_time_nozeros <- reading_time[reading_time$region7 != 0,]

library(lme4)
m<-lmer(log(region7) ~ dat+adj+int + 
          (dat+adj+int|subj) + 
          (dat+adj+int|item), 
        data=reading_time_nozeros)

## ----setparameters,echo=FALSE,include=TRUE-------------------------------
## set true parameter values:
beta<-round(summary(m)$coefficients[,1],4)
sigma_e<-round(attr(VarCorr(m),"sc"),4)
subj_ranefsd<-round(attr(VarCorr(m)$subj,"stddev"),4)
subj_ranefcorr<-round(attr(VarCorr(m)$subj,"corr"),1)
## choose some intermediate values for correlations:
corr_matrix<-(diag(4) + matrix(rep(1,16),ncol=4))/2

## assemble variance matrix for subjects:
Sigma_u<-SIN::sdcor2cov(stddev=subj_ranefsd,corr=corr_matrix)

item_ranefsd<-round(attr(VarCorr(m)$item,"stddev"),4)

## assemble variance matrix for items:
Sigma_w<-SIN::sdcor2cov(stddev=item_ranefsd,corr=corr_matrix)

## ----simulatedata,echo=FALSE,include=FALSE,cache=TRUE--------------------
set.seed(4321)
nsim<-100
## effect size ranging from 30 to 80 ms:
(beta2<-c(0.056,0.095,0.155))
estc1<-tvalsc1<-tvalsc2<-tvalsc3<-matrix(rep(NA,nsim*length(beta2)),ncol=nsim)
failed<-matrix(rep(0,nsim*length(beta2)),ncol=nsim)
for(j in 1:length(beta2)){
for(i in 1:nsim){
  beta[2]<-beta2[j]
  dat<-gen_fake_lnorm2x2(nitem=24,
                         nsubj=28,
                       beta=beta,
                       Sigma_u=Sigma_u,
                       Sigma_w=Sigma_w,
                      sigma_e=sigma_e)

## no correlations estimated to avoid convergence problems: 
## analysis done after log-transforming:  
m<-lmer(log(rt) ~ c1+c2+c3 + (c1+c2+c3||subj) + 
          (c1+c2+c3||item), data=dat)
estc1[j,i]<-round(summary(m)$coefficients[2,1],4)
## ignore failed trials
if(any( grepl("failed to converge", m@optinfo$conv$lme4$messages) )){
  failed[j,i]<-1
} else{
tvalsc1[j,i]<-summary(m)$coefficients[2,3]
tvalsc2[j,i]<-summary(m)$coefficients[3,3]
tvalsc3[j,i]<-summary(m)$coefficients[4,3]
}}}
## proportion of convergence failures:
rowMeans(failed)

## ----powercalculation,eval=TRUE,echo=FALSE,include=FALSE-----------------
pow<-rep(NA,length(beta2))
for(k in 1:length(beta2)){
  pow[k]<-mean(abs(tvalsc1[k,])>2,na.rm=TRUE)
}

## ----figagrmtattrn,echo=FALSE,fig.width=6,fig.height=4,include=TRUE,warning=FALSE----
posteriors$data_model<-c(rep("data",12),"model")

posteriors$expt<-factor(posteriors$expt,levels=c(levels(posteriors$expt),"model"))

pd<-position_dodge(0.6)

plot_targetmismatch<-ggplot(posteriors, aes(x=expt,
                                            y=mean, 
                             group=expt,shape=data_model)) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=.25, size=.5, position=pd) +
  scale_x_discrete(labels=levels(posteriors$expt[1:12]))+
  labs(title="Subject-verb dependencies") +
  #scale_x_continuous(limits = c(0, 11))+
  xlab("Experiment id / LV05 model")+
  ylab("Estimates (ms)")+
    theme_bw()+
  scale_shape_discrete(labels=c("data",
                                "model"))+
  theme(legend.position=c(0.8,0.2)) +
  #theme(legend.title=element_text("data / model"))+
  #guide_legend(title="data / model")+
  geom_hline(yintercept=0)+
  geom_point(position=pd, size=2)
#+coord_flip()

num_subj<-as.vector(posteriors$n_subj)
num_item<-as.vector(posteriors$n_item)

annsze<-3
plot_targetmismatch+annotate("text", x = 1:13, y = 25, label = num_subj,size=annsze)+
  annotate("text", x = 1:13, y = 15, label = num_item,size=annsze)+annotate("text", x=1,y= 30, label="subj",size=annsze)+annotate("text", x=1,y= 20, label="item",size=annsze)+    magnifytext()

## ----rawrtanalyses,eval=FALSE,echo=FALSE,include=FALSE-------------------
## summary(SPRLK1critm1)
## summary(SPRLK1postcritm1post)
## summary(m2)
## summary(m6)
## summary(m_E4ETTFTrawcrit)
## summary(m_E4ETTFTrawpostcrit)
## summary(mE6etcritraw)
## summary(mE6etpostcritraw)
## summary(msprmergedcritraw)
## summary(msprmergedpostcritraw)

## ----wordlenfreq,cache=TRUE,echo=FALSE,eval=TRUE,include=TRUE,warning=FALSE,message=FALSE----
# Sanity check on filler items: Is there a frequency effect and word length effect as expected?

# Eyetracking data
e2 <- read.table('../data/Exp2fillersFreq.txt', header=TRUE)
#summary(e2)


e2$cfreq<-scale(e2$type_logFreq,scale=FALSE)
e2$clen<-scale(e2$len,scale=FALSE)

#xtabs(~itemid+cfreq,e2)

e2ma<-lmer(FPRT~cfreq+clen+(1+cfreq+clen||subject)+(1+cfreq+clen||itemid),e2)
#summary(e2ma)

e4 <- read.table('../data/Exp4fillersFreq.txt', header=TRUE)

e4$cfreq<-scale(e4$type_logFreq,scale=FALSE)
e4$clen<-scale(e4$len,scale=FALSE)

e4ma<-lmer(FPRT~cfreq+clen+(1+cfreq+clen||subject)+(1+cfreq+clen||itemid),e4)
#summary(e4ma)


e6 <- read.table('../data/Exp6fillersFreq.txt', header=TRUE)
e6$cfreq<-scale(e6$type_logFreq,scale=FALSE)
e6$clen<-scale(e6$len,scale=FALSE)

e6ma<-lmer(FPRT~cfreq+clen+(1+cfreq+clen||subject)+(1+cfreq+clen||itemid),e6)
#summary(e6ma)

e7 <- read.table('../data/Exp7fillersFreq.txt', header=TRUE)

e7$cfreq<-scale(e7$type_logFreq,scale=FALSE)
e7$clen<-scale(e7$len,scale=FALSE)

e7ma<-lmer(FPRT~cfreq+clen+(1+cfreq+clen||subject)+(1+cfreq+clen||itemid),e7)
#summary(e7ma)

est<-rbind(summary(e2ma)$coefficients[2:3,c(1,2,3)],
summary(e4ma)$coefficients[2:3,c(1,2,3)],
summary(e6ma)$coefficients[2:3,c(1,2,3)],
summary(e7ma)$coefficients[2:3,c(1,2,3)])

Expt<-paste(rep("Expt",8),rep(1:4,each=2),sep=" ")

est<-data.frame(Expt=Expt,Predictor=rep(c("Freq","Len"),4),est)

#xtable(est,digits=c(1,1,1,2,2,2))

## ----saveallanalyses,echo=FALSE,eval=FALSE-------------------------------
## ## just for future analyses:
## save.image(file="allanalyses.Rda")

