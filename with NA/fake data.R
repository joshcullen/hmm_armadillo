rm(list=ls(all=TRUE))
set.seed(1)

#parameters
nobs=10000
ngroup=4
mu.SL.true=mu.SL=c(-2,0,1,2)
sd.SL.true=sd.SL=c(1,0.5,0.5,0.5)
mu.TA.true=mu.TA=c(2,0.5,-1,-2)
sd.TA.true=sd.TA=c(0.3,2,0.3,0.3)

#get group assignments
tmp=rmultinom(nobs,size=1,prob=rep(1/ngroup,ngroup))
ind=rep(NA,nobs)
for (i in 1:nobs){
  ind[i]=which(tmp[,i]==1)
}

#get SL
fim=data.frame(log.SL=rep(NA,nobs),
               logit.TA=rep(NA,nobs))
fim$log.SL=rnorm(nobs,mean=mu.SL[ind],sd=sd.SL[ind])

#look at SL
par(mfrow=c(2,2))
rango=range(exp(fim$log.SL))
for (i in 1:ngroup) hist(exp(fim$log.SL[ind==i]),xlim=rango)

#get TA
fim$logit.TA=rnorm(nobs,mean=mu.TA[ind],sd=sd.TA[ind])

#look at TA
tmp=exp(fim$logit.TA)
prob=tmp/(1+tmp)
par(mfrow=c(2,2))
rango=c(0,1)
for (i in 1:ngroup) hist(prob[ind==i],xlim=rango)

#plug in NA's in TA
# ind=rbinom(nobs,size=1,prob=0.002); sum(ind==1)
# fim[ind==1,2]=NA

#export results
setwd('U:\\GIT_models\\hmm_armadillo')
write.csv(fim,'fake data.csv',row.names=F)