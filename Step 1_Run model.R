## Run the Bayesian HMM

set.seed(1)

library(tidyverse)
library(lubridate)
library(tictoc)


source('gibbs hmm aux functions.R')
source('gibbs hmm main function.R')
source('helper functions.R')

## Import and prep data
dat<- read.csv("Modified Armadillo Data.csv", header = T, sep = ",")
dat$id<- as.character(dat$id)
dat$date<- as_datetime(dat$date)
dat<- dat %>% 
  rename(step = dist, angle = rel.angle)

#filter out burrow locs by step length
dat<- dat %>% 
  filter(step != 0)

#filter obs for primary time step (300 s)
dat<- bayesmove::round_track_time(dat, "id", "dt", "date", 300, 60, time.zone = "UTC")
dat<- dat[dat$dt == 300,]



#Explore distribution of SL and TA
ggplot(dat, aes(step)) +
  geom_density(fill = "cadetblue") +
  theme_bw()

ggplot(dat, aes(angle)) +
  geom_density(fill = "firebrick") +
  theme_bw()





#Log-transform SL and logit-transform abs. value of pi-scaled TA
dat<- transform_data(dat = dat)

#remove all obs where log.SL or logit.TA are NA
dat<- dat[!is.na(dat$log.SL),]
dat<- dat[!is.na(dat$logit.TA),]

#Explore transformed distribution of SL and TA
ggplot(dat, aes(log.SL)) +
  geom_density(fill = "cadetblue") +
  theme_bw()

ggplot(dat, aes(logit.TA)) +
  geom_density(fill = "firebrick") +
  theme_bw()




## Run model

#priors
var.mu=100
sig2.a=0.1; sig2.b=0.1
gamma1=0.1

#initialize parameters
max.group=7

#MCMC stuff
ngibbs=3000
nburn=ngibbs/2

tic()
mod=hmm.main.function(dat=dat,var.mu=var.mu,sig2.a=sig2.a,sig2.b=sig2.b,
                      gamma1=gamma1,max.group=max.group,
                      ngibbs=ngibbs,nburn=nburn)
toc()
#takes 33 min to run when ngibbs = 3000 and max.group = 7


## Inspect results

#check traceplot of log-likelihood
plot(mod$llk, type='l')

#inpsect boxplot of thetas
boxplot(mod$theta[(nburn + 1):ngibbs,])
colMeans(mod$theta[(nburn + 1):ngibbs,]) %>% cumsum()  #suggests 3 active states

MAP.est<- which.max(mod$llk[(nburn + 1):ngibbs]) + nburn







behav.params<- mod %>% 
  keep(names(.) == names(.)[1:4]) %>% 
  map(., as.data.frame) %>% 
  map(., select, 1:3) %>% 
  map(., slice, MAP.est) %>% 
  bind_rows() %>% 
  t() %>% 
  data.frame() %>% 
  mutate(behav = 1:3)

names(behav.params)[1:4]<- names(mod)[1:4]



plot_data_SL<- 
  pmap_df(behav.params[,-c(1,3)],
          function(mu.sk, sig2.sk, behav) {
            tibble(y = exp(rnorm(1000, mean = mu.sk, sd = sqrt(sig2.sk))),
                   behav = behav)
          })
plot_data_SL$behav<- factor(plot_data_SL$behav)


plot_data_TA<- 
  pmap_df(behav.params[,-c(2,4)],
          function(mu.ak, sig2.ak, behav) {
            tibble(y = plogis(rnorm(1000, mean = mu.ak, sd = sqrt(sig2.ak))) * pi,
                   behav = behav)
          })
plot_data_TA$behav<- factor(plot_data_TA$behav)


# Plot overlapping densities per behavior (step lengths)
ggplot(data = plot_data_SL, aes(color = behav)) +
  geom_density(aes(group = behav, x = y), alpha = 0.4, size = 1) + 
  labs(x = "\nStep Length (m)", y = "Density\n") +
  scale_color_brewer("", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))


# Plot turning angle
ggplot(data = plot_data_TA, aes(color = behav)) +
  geom_density(aes(group = behav, x = y), alpha = 0.4, size = 1) + 
  labs(x = "\nTurning Angle (rad)", y = "Density\n") +
  scale_color_brewer("", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))




# Viz mapped behaviors
dat$behav<- factor(mod$z.k)

ggplot() +
  geom_point(data = dat, aes(x, y, color = behav), size = 2,
             alpha = 0.6) +
  theme_bw() +
  scale_color_brewer("Behavior", palette = "Dark2") +
  facet_wrap(~id, scales = "free")
