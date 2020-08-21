## Run the Bayesian HMM

set.seed(1)

library(tidyverse)
library(lubridate)
library(tictoc)
library(lunar)
library(raster)
library(RColorBrewer)


source('gibbs hmm aux functions.R')
source('gibbs hmm main function.R')
source('helper functions.R')

## Import and prep data
dat<- read.csv("Modified Armadillo Data.csv", header = T, sep = ",")
dat$id<- as.character(dat$id)
dat$date<- as_datetime(dat$date)
dat<- dat %>% 
  rename(step = dist, angle = rel.angle)
dat<- bayesmove::round_track_time(dat, "id", "dt", "date", 300, 60, time.zone = "UTC")
dat$time1<- 1:nrow(dat)  #create index value for merge of results after filtering

#filter out burrow locs by step length
dat2<- dat %>% 
  filter(step != 0)

#filter obs for primary time step (300 s)
dat2<- dat2[dat2$dt == 300,]



#Explore distribution of SL and TA
ggplot(dat2, aes(step)) +
  geom_density(fill = "cadetblue") +
  theme_bw()

ggplot(dat2, aes(angle)) +
  geom_density(fill = "firebrick") +
  theme_bw()





#Log-transform SL and logit-transform abs. value of pi-scaled TA
dat2<- transform_data(dat = dat2)

#remove all obs where log.SL or logit.TA are NA
dat2<- dat2[!is.na(dat2$log.SL),]
dat2<- dat2[!is.na(dat2$logit.TA),]

#Explore transformed distribution of SL and TA
ggplot(dat2, aes(log.SL)) +
  geom_density(fill = "cadetblue") +
  theme_bw()

ggplot(dat2, aes(logit.TA)) +
  geom_density(fill = "firebrick") +
  theme_bw()




## Run model

#priors
var.mu=100
sig2.a=0.1; sig2.b=0.1
gamma1=0.1

#initialize parameters
max.group=2

#MCMC stuff
ngibbs=3000
nburn=ngibbs/2

tic()
mod=hmm.main.function(dat=dat2,var.mu=var.mu,sig2.a=sig2.a,sig2.b=sig2.b,
                      gamma1=gamma1,max.group=max.group,
                      ngibbs=ngibbs,nburn=nburn)
toc()
#takes 29 min to run when ngibbs = 3000 and max.group = 2


## Inspect results

#check traceplot of log-likelihood
plot(mod$llk, type='l')

#inpsect boxplot of thetas
boxplot(mod$theta[(nburn + 1):ngibbs,])
colMeans(mod$theta[(nburn + 1):ngibbs,]) %>% cumsum()  #suggests 2 active states

MAP.est<- which.max(mod$llk[(nburn + 1):ngibbs]) + nburn







behav.params<- mod %>% 
  keep(names(.) == names(.)[1:4]) %>% 
  map(., as.data.frame) %>% 
  map(., select, 1:2) %>% 
  map(., slice, MAP.est) %>% 
  bind_rows() %>% 
  t() %>% 
  data.frame() %>% 
  mutate(behav = 1:2)

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
dat2$behav<- factor(mod$z.k)

# All behavs
ggplot() +
  geom_point(data = dat2, aes(x, y, color = behav), size = 2,
             alpha = 0.6) +
  theme_bw() +
  scale_color_brewer("Behavior", palette = "Dark2") +
  facet_wrap(~id, scales = "free")

# Foraging
ggplot() +
  geom_point(data = dat2 %>% filter(behav == 1), aes(x, y),
             color = RColorBrewer::brewer.pal(3, "Dark2")[1], size = 2,
             alpha = 0.6) +
  theme_bw() +
  scale_color_brewer("Behavior", palette = "Dark2") +
  facet_wrap(~id, scales = "free")

# Transit
ggplot() +
  geom_point(data = dat2 %>% filter(behav == 2), aes(x, y),
             color = RColorBrewer::brewer.pal(3, "Dark2")[2], size = 2,
             alpha = 0.6) +
  theme_bw() +
  scale_color_brewer("Behavior", palette = "Dark2") +
  facet_wrap(~id, scales = "free")




### Merge with original dataset

dat.merge<- left_join(dat, dat2[,c("behav", "time1")], by = "time1")
dat.merge$behav<- ifelse(is.na(dat.merge$behav) & dat.merge$step == 0, 3, dat.merge$behav)
dat.merge$behav<- factor(dat.merge$behav, levels = c(3,1,2))
levels(dat.merge$behav)<- c("Burrow","Foraging","Transit")


# All behavs and all IDs
ggplot() +
  geom_path(data = dat.merge, aes(x, y, color = behav, group = id), size = 0.75) +
  theme_bw() +
  scale_color_brewer("Behavior", palette = "Dark2", na.value = "grey75") +
  facet_wrap(~id, scales = "free")


# Focus on tm15
ggplot() +
  geom_path(data = dat.merge %>% filter(id == "tm15"), aes(x, y, color = behav, group = id),
            size = 0.75, alpha = 0.6) +
  geom_point(data = dat.merge %>% filter(id == "tm15", behav == "Burrow"),
             aes(x, y, color = behav), size = 2) +
  theme_bw() +
  scale_color_brewer("Behavior", palette = "Dark2", na.value = "grey75") +
  coord_equal()

ggplot() +
  geom_path(data = dat.merge %>% filter(id == "tm15"), aes(x, y, group = id), color = "grey75",
            size = 0.75) +
  geom_point(data = dat.merge %>% filter(id == "tm15", behav == "Burrow"),
             aes(x, y, color = behav), size = 2) +
  theme_bw() +
  scale_color_brewer("Behavior", palette = "Dark2", na.value = "grey75") +
  coord_equal()



## Extract time of day for each obs

dat.merge2<- dat.merge %>% 
  drop_na(date)
dat.merge2$date<- as.POSIXct(strptime(dat.merge2$date, format = "%Y-%m-%d %H:%M:%S"),
                             tz="Brazil/West")
dat.list<- bayesmove::df_to_list(dat.merge2, "id")

tod_res<- dat.list %>% 
  map(., ~make_track(tbl = ., x, y, date, id, crs = CRS("+init=epsg:32721"))) %>% 
  map(., time_of_day, include.crepuscule = TRUE) %>% 
  bind_rows()

dat.merge2$TOD<- tod_res$tod_


## Summarize behavior props per day

daily.behav<- dat.merge2 %>% 
  drop_na(behav) %>% 
  group_by(id, date = lubridate::ymd(strptime(date, format = "%Y-%m-%d")), behav, .drop = F) %>% 
  count() %>% 
  group_by(id, date) %>%
  mutate(prop = n/sum(n)) %>% 
  ungroup()

#check that all days have props that sum to 1
(foo<- daily.behav %>% 
  group_by(id, date) %>% 
  summarise(prop=sum(prop)))






fill_date_gaps = function(dat) {
  
  id1<- unique(dat$id)
  
  #create seq of dates (1 per day)
  date.seq<- seq(range(dat$date)[1], range(dat$date)[2], by = 1)
  date.seq.df<- data.frame(date = rep(date.seq, each = 3))
  
  #join data frames to include gaps
  merge.df<- right_join(dat, date.seq.df, by = "date")
  ind<- which(is.na(merge.df$id))
  merge.df$behav[ind]<- rep(levels(dat$behav), length(ind)/3)
  merge.df$id[ind]<- rep(id1, length(ind))
  
  #add columns to define upper and lower bounds for geom_ribbon plot
  merge.df$ymax <-merge.df$prop
  merge.df$ymin <- 0
  behav.ind <- rev(levels(merge.df$behav))
  for ( i in 2:length(behav.ind) ) {
    behav_i <- merge.df$behav==behav.ind[i]
    behav_i_1 <- merge.df$behav==behav.ind[i-1]
    merge.df$ymin[behav_i] <- merge.df$ymax[behav_i_1]
    merge.df$ymax[behav_i] <- merge.df$ymin[behav_i] + merge.df$ymax[behav_i]
  }
  
  merge.df
}


daily.list<- bayesmove::df_to_list(daily.behav, "id")
daily.behav2<- purrr::map(daily.list, fill_date_gaps) %>% 
  bind_rows()


ggplot(daily.behav2) +
  geom_ribbon(aes(x=date, ymin=ymin, ymax=ymax, fill = behav), color = "black", size = 0.25) +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")



ggplot(daily.behav2 %>% filter(id == "tm15")) +
  geom_ribbon(aes(x=date, ymin=ymin, ymax=ymax, fill = behav), color = "black", size = 0.25) +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")



#Assess relationship between behavior proportions and lunar illumination

daily.behav2$lunar<- lunar.illumination(daily.behav2$date, shift = 12)

ggplot(daily.behav2, aes(prop, lunar)) +
  geom_point(size = 2, alpha = 0.7, show.legend = F) +
  geom_smooth(method = "gam") +
  labs(x="Daily Proportion of Behavior", y="Lunar Illumination") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~ behav)


ggplot(daily.behav2 %>% filter(id == "tm15"), aes(prop, lunar)) +
  geom_point(size = 2, alpha = 0.7, show.legend = F) +
  geom_smooth(method = "gam") +
  labs(x="Daily Proportion of Behavior", y="Lunar Illumination") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~ behav)



### Import Environ Covars

setwd("~/Documents/Snail Kite Project/Data/armadillos/Environ Data")

rast<- dir(getwd(), "*.tif$")
for (i in rast) assign(i, raster(i))

#Need to project rasters same as tracks
EucDist_cerc_Copy.tif<- projectRaster(EucDist_cerc_Copy.tif, crs = "+init=epsg:32721")
ex.ras<- raster(ext = extent(EucDist_cerc_Copy.tif), crs = "+init=epsg:32721", res = 30)
dist2rdN_30m<- resample(EucDist_cerc_Copy.tif, ex.ras, method = "bilinear")
dist2rdN_30m.df<- as.data.frame(dist2rdN_30m, xy=T)

#resample DEMs to 30m from 18m; they will now share the same dimensions and extent
dem.N<- resample(dem_N.tif, dist2rdN_30m, method = "bilinear")
compareRaster(dist2rdN_30m, dem.N)
dem.N.df<- as.data.frame(dem.N, xy=T)


setwd("~/Documents/Snail Kite Project/Data/armadillos/NDVI")

#load files
ndvi.filenames<- list.files(getwd(), pattern = "*.grd$")
ndvi.N<- brick(ndvi.filenames[1])

#change extent and dimensions of RasterBricks using resample()
ndvi.N<- resample(ndvi.N, dist2rdN_30m, method = "bilinear")
compareRaster(dist2rdN_30m, ndvi.N)

ndvi.N.mean<- mean(ndvi.N, na.rm = T)
ndvi.N.mean.df<- as.data.frame(ndvi.N.mean, xy=T)



#Plot over dist2rd
ggplot() +
  geom_tile(data = dist2rdN_30m.df, aes(x, y, fill = EucDist_cerc_Copy)) +
  scale_fill_viridis_c("Distance to Road (m)", option = "inferno", na.value = "n") +
  geom_path(data = dat.merge %>% filter(id == "tm15"), aes(x, y, color = behav, group = id),
            size = 0.75, alpha = 0.6) +
  geom_point(data = dat.merge %>% filter(id == "tm15", behav == "Burrow"),
             aes(x, y, color = behav), size = 2) +
  scale_color_brewer("Behavior", palette = "Dark2", na.value = "grey75") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))



#Plot over elev
ggplot() +
  geom_tile(data = dem.N.df, aes(x, y, fill = dem_N)) +
  scale_fill_gradientn("Elevation (m)", colours = terrain.colors(10), na.value = "n") +
  geom_path(data = dat.merge %>% filter(id == "tm15"), aes(x, y, color = behav, group = id),
            size = 0.75, alpha = 0.6) +
  geom_point(data = dat.merge %>% filter(id == "tm15", behav == "Burrow"),
             aes(x, y, color = behav), size = 2) +
  scale_color_brewer("Behavior", palette = "Dark2", na.value = "grey75") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))




#Plot over NDVI
green_colors <- brewer.pal(9, "YlGn") %>%
  colorRampPalette()

ggplot() +
  geom_tile(data = ndvi.N.mean.df, aes(x, y, fill = layer), alpha = 0.8) +
  scale_fill_gradientn(name = "NDVI", limits = c(0,1), colours = green_colors(20),
                       na.value = "n") +
  geom_path(data = dat.merge %>% filter(id == "tm15"), aes(x, y, color = behav, group = id),
            size = 0.75, alpha = 0.6) +
  geom_point(data = dat.merge %>% filter(id == "tm15", behav == "Burrow"),
             aes(x, y, color = behav), size = 2) +
  scale_color_brewer("Behavior", palette = "Dark2", na.value = "grey75") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
