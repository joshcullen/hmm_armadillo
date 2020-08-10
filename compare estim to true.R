plot(theta,type='h')

compare1=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango,col='red')
}

estim.mu.ak=store.mu.ak[ngibbs,1:4]
estim.mu.sk=store.mu.sk[ngibbs,1:4]
estim.sd.sk=sqrt(store.sig2.sk[ngibbs,1:4])
estim.sd.ak=sqrt(store.sig2.ak[ngibbs,1:4])

ordem=c(1,4,3,2)
par(mfrow=c(2,2))
compare1(true=mu.SL.true,estim=estim.mu.sk[ordem])
compare1(true=mu.TA.true,estim=estim.mu.ak[ordem])
compare1(true=sd.SL.true,estim=estim.sd.sk[ordem])
compare1(true=sd.TA.true,estim=estim.sd.ak[ordem])
