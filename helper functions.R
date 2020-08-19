logit = function(p) {
  log(p / (1-p))
}

#------------------------------
transform_data = function(dat) {
  #step = step length
  #angle = turning angle
  
  #adjust any 0 step lengths and any values of 0 or pi for turning angles
  ind.SL<- which(dat$step == 0)
  if (length(ind.SL)) dat[ind.SL,]$step<- 0.01 
  
  ind.TA.0<- which(dat$angle == 0)
  if (length(ind.TA.0)) dat[ind.TA.0,]$angle<- 0.01
  
  ind.TA.pi<- which(abs(round(dat$angle, 3)) == round(pi, 3))
  if (length(ind.TA.pi)) dat[ind.TA.pi,]$angle<- round(pi - 0.01, 3)
  
  
  #log-transform SL (dist)
  dat$log.SL<- log(dat$step)
  
  #logit-transform TA abs(rel.angle) scaled to pi
  dat$logit.TA<- logit(abs(dat$angle / pi))
  
  return(dat)
}
