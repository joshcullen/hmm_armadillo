rm(list=ls(all=TRUE))
set.seed(3)

#import data
source('gibbs hmm aux functions.R')
source('gibbs hmm main function.R')
dat=read.csv('fake data.csv', as.is=T)
nobs=nrow(dat)

#priors
var.mu=100
sig2.a=0.1; sig2.b=0.1
gamma1=0.1

#initialize parameters
max.group=10

#MCMC stuff
ngibbs=3000
nburn=ngibbs/2

mod=hmm.main.function(dat=dat,var.mu=var.mu,sig2.a=sig2.a,sig2.b=sig2.b,
                      gamma1=gamma1,max.group=max.group,
                      ngibbs=ngibbs,nburn=nburn)
