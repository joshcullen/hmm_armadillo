plot(mod$theta[ngibbs,],type='h')

compare1=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango,col='red')
}

estim.mu.ak=mod$mu.ak[ngibbs,1:4]
estim.mu.sk=mod$mu.sk[ngibbs,1:4]
estim.sd.sk=sqrt(mod$sig2.sk[ngibbs,1:4])
estim.sd.ak=sqrt(mod$sig2.ak[ngibbs,1:4])

ordem=c(1,3,4,2)
par(mfrow=c(2,2))
compare1(true=mu.SL.true,estim=estim.mu.sk[ordem])
compare1(true=mu.TA.true,estim=estim.mu.ak[ordem])
compare1(true=sd.SL.true,estim=estim.sd.sk[ordem])
compare1(true=sd.TA.true,estim=estim.sd.ak[ordem])
