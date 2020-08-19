rm(list=ls())
library('viridis')

#get data
setwd('U:\\independent studies\\tsbp\\nina')
dat=read.csv('Armadillo edited.csv',as.is=T)

#get modeling results
setwd('U:\\independent studies\\tsbp\\nina\\hmm results')
#look at convergence
llk=read.csv('llk.csv',as.is=T)
plot(llk$V1,type='l')
seq1=500:nrow(llk)
plot(llk$V1[seq1],type='l')

max.llk=read.csv('max llk.csv',as.is=T)
abline(h=max.llk$x,col='red')
ind.max=which(llk$V1==max.llk$x)

#look at number of groups
theta=read.csv('theta.csv',as.is=T)
round(theta[ind.max,],2)

theta.summary=round(apply(theta[seq1,],2,mean),2); theta.summary
plot(theta.summary,type='h')
plot(unlist(theta[nrow(theta),]),type='h') #max of 5 groups
ngr.max=5

#get parameters
mu.ak=read.csv('mu ak.csv',as.is=T)[seq1,1:ngr.max]
mu.sk=read.csv('mu sk.csv',as.is=T)[seq1,1:ngr.max]
sig2.ak=read.csv('sig2 ak.csv',as.is=T)[seq1,1:ngr.max]
sig2.sk=read.csv('sig2 sk.csv',as.is=T)[seq1,1:ngr.max]

#change order from directed to less directed movement
# seq2=c(4,3,2,5,1)
# mu.ak=mu.ak[,seq2]
# mu.sk=mu.sk[,seq2]
# sig2.ak=sig2.ak[,seq2]
# sig2.sk=sig2.sk[,seq2]

#look at SL
nsim=nrow(mu.sk)
fim=list(); xmax=ymax=-Inf
for (i in 1:ngr.max){
  tmp=rnorm(nsim,mean=mu.sk[,i],sd=sqrt(sig2.sk[,i]))
  tmp1=exp(tmp)
  tmp2=density(tmp1,from=0)
  fim[[i]]=tmp2
  if (ymax<max(tmp2$y)) ymax=max(tmp2$y)
  if (xmax<max(tmp2$x)) xmax=max(tmp2$x)
}

#plot results
setwd('U:\\independent studies\\tsbp\\nina\\hmm results')
xmax=100
png('hmm results SL.png',width=700,height=700)
plot(NA,NA,xlim=c(0,xmax),ylim=c(0,ymax),main='Step length',xlab='',ylab='',
     cex.axis=1.5,cex.main=2)
for (i in 1:ngr.max){
  lines(fim[[i]]$x,fim[[i]]$y,col=i,lwd=2)
}
legend(40,ymax,paste('Groups ',1:5),col=1:5,lty=1,cex=2,lwd=2)
dev.off()

#look at TA
nsim=nrow(mu.ak)
fim=list(); xmax=ymax=-Inf
for (i in 1:ngr.max){
  tmp=rnorm(nsim,mean=mu.ak[,i],sd=sqrt(sig2.ak[,i]))
  tmp1=exp(tmp)
  tmp2=density(pi*tmp1/(1+tmp1),from=0,to=pi)
  fim[[i]]=tmp2
  if (ymax<max(tmp2$y)) ymax=max(tmp2$y)
  if (xmax<max(tmp2$x)) xmax=max(tmp2$x)
}

png('hmm results TA.png',width=700,height=700)
plot(NA,NA,xlim=c(0,xmax),ylim=c(0,ymax),main='Turning angle magnitude',
     cex.axis=1.5,cex.main=2,xlab='',ylab='')
for (i in 1:ngr.max){
  lines(fim[[i]]$x,fim[[i]]$y,col=i,lwd=2)
}
legend(2,ymax,paste('Groups ',1:5),col=1:5,lty=1,cex=2,lwd=2)
dev.off()
#group 4 seems to be dispersal: biggest step lengths and TA closer to 0