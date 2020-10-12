# exploratory analysis of NASH data


rm(list=ls()) 
library(nleqslv)
library(extraDistr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(rstan)
library(bayesplot)
library(loo)
library(VGAM)
library(evd)
# library(moments)

# datadir = file.path('..', '..','hbev', 'data', 'Data_GHCN_Daily', 'extracted_csv')
datadir = file.path('..', 'data', 'Data_GHCN_Daily', 'extracted_csv')
  dfpos = read.csv(file.path(datadir, '..',  'list_ushcn_stats.csv'))
  # dflakes <- subset(dfpos, dfpos$LATITUDE > 45 & dfpos$LONGITUDE > -94)
  dflakes <- subset(dfpos, dfpos$STATE == 'WI ')
  dften <- subset(dfpos, dfpos$STATE == 'TN ')



# library(tseries)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



outplot = file.path('..', 'nash_output_figures')
dir.create(outplot, showWarnings = FALSE)
outdata = file.path('..', 'nash_output', 'outdata')
dir.create(outdata, showWarnings = FALSE)


source("nashfun.R")    # main functions for data analysis
# source("hbev_module.R")    # main functions for data analysis
# source("hbev_functions.R") # other functions

iter = 2000
chains = 4
ndraws = iter*chains/2
thresh = 1

# load station data
# datadir = file.path('..','data', 'Data_GHCN_Daily', 'extracted_csv')
filenames = list.files(datadir)

# station_id = 'USW00094728' # NEW YORK CENTRAL PARK
# station_id = 'USC00044232' # TEXAS, HOUSTON
station_id = 'USC00045532' # stat numj 77
# station_id = 'USC00311677' # NORTH CAROLINA, CHAPEL HILL
# station_id = 'USC00478027' # WISCONSIN, SPOONER
# station_id = 'USC00475516' # WISCONSIN, MINOQUA
# station_id = 'USC00207812' # MICHIGAN, STAMBAUGH
# station_id = 'USC00409502' # TENNESSEE, WAYNESBORO
  
filepath = file.path(datadir, sprintf('%s.csv', station_id)) 


  
mystation = strsplit(as.character(dfpos$NAME[dfpos$ID == station_id]), ' ')[[1]][1]
filepath = file.path(datadir, sprintf('%s.csv', station_id)) # NYCP

df = read.csv(filepath)
df$YEAR = floor(df$DATE/10000)
df$MONTH = floor((df$DATE-df$YEAR*10000)/100)
df = subset(df, df$YEAR >= 1948 & df$YEAR <= 2007)
df = subset(df, df$MONTH %in% c(6, 7, 8)) # summer precipitation only
df$PRCP = df$PRCP/10 # to millimeter accumulations
# df$PRCP = pmax(df$PRCP-10, 0)
prcp_years = unique(df$YEAR)

# load nash data
ynash = read.table(file.path( 'nash_pos_data', 'ridge_lat_data.txt'), header=FALSE, col.names = c('YEAR', 'LAT'))
xnash = read.table(file.path( 'nash_pos_data', 'ridge_lon_data.txt'), header=FALSE, col.names = c('YEAR', 'LON'))
xnash$LON = -xnash$LON
plot(xnash$LON, ynash$LAT)

# dflakes = subset(df, df$)

# plot(ynash$LAT, )


xnash = subset(xnash, xnash$YEAR %in% prcp_years)
ynash = subset(ynash, ynash$YEAR %in% prcp_years)


xnn = (xnash$LON - mean(xnash$LON))/sd(xnash$LON)
ynn = (ynash$LAT - mean(ynash$LAT))/sd(ynash$LAT)

# xnn = xnash$LON
# ynn = ynash$LAT




# compute annual maxima statistics
Nt = 92
res = table_max(df, Nt=Nt, thresh = thresh)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue', 'red'))

#This adds a column of color values
# based on the y values
datCol <- rbPal(10)[as.numeric(cut(ynn,breaks = 10))]

# plot(dat$x,dat$y,pch = 20,col = dat$Col)
plot(xnn, res$max, col=datCol, pch=20)

# plot(ynn[xnn <0], res$max[xnn< 0], col=datCol, pch=20)
plot(ynn*xnn, res$max, col=datCol, pch=20)


# plot(ynn[xnn>1], res$max[xnn > 1], col=datCol, pch=20)


datCol2 <- rbPal(20)[as.numeric(cut(res$totals, breaks = 20))]
plot(xnn, ynn, col=datCol2, pch=10)

ndf = data.frame( list(xnn = xnn, ynn = ynn, max = res$max, 
                       N = res$N, totals = res$totals, 
                       means =res$totals/res$N, sdwets = res$sdwets, 
                       skews = res$skews))


# Nt = 92
N = ndf$N
M = length(N)
XN = ndf$xnn
YN = ndf$ynn
Mgen = 60 # number of years of data to generate
XNgen = rep(-1.0, Mgen) # west
YNgen = rep(1.0, Mgen) # north

# for comparison
XNgenS = rep(-1.0, Mgen) # west
YNgenS = rep(-1.0, Mgen) # south
XNgenN = rep(-1.0, Mgen) # west
YNgenN = rep(1.0, Mgen) # north
nwets_data = list(M = M, Mgen = Mgen, Nt = Nt, N=N, XN = XN, 
                  YN = YN, XNgen = XNgen, YNgen = YNgen)








gendata <- gener_2smc(wetf = 0.5, corr = 0.8, nyears=1, nobs=10000)
# trp <- trans_prob(gendata, thresh = 0)
# print(trp$corr)
# print(trp$wetf)



trans_prob_mc <- trans_prob(res$data, thresh = 0)

# test with synthetic data::
# gendata <- gener_2smc(wetf = 0.6, corr = 0.9, nyears=M, nobs=Nt)
# res$data = gendata


# thresh = 0
totnwets = length(res$data[res$data > 0])
nash_data = list(M = length(ndf$means),
                 Mgen = Mgen,
                  Nt = Nt,
                  y = res$data,
                  # WETF = trans_prob_mc$wetf,
                  # CORR = trans_prob_mc$corr,
                  XN = ndf$xnn, YN = ndf$ynn,
                  XNgen = XNgen, YNgen = YNgen,
                  totnwets = totnwets)




# sum(gendata)/length(gendata)

# length(res$data[res$data > 0])

fitmc_no <- stan(file = 'nash_nj_markov_depno.stan',  data = nash_data,
iter = iter, chains = 4)
fitmc_yy <- stan(file = 'nash_nj_markov_depyy.stan',  data = nash_data,
iter = iter, chains = chains)
fitmc <- stan(file = 'nash_nj_markov_depxy.stan',  data = nash_data,
iter = iter, chains = chains)

sim_mc <- rstan::extract(fitmc)

# print(mean(sim_mc$corr))
# print(mean(sim_mc$pwet))

plot(fitmc, pars=c('pwet', 'corr'))



pairs(fitmc_no, pars=c('pwet', 'corr'))
pairs(fitmc, pars=c('pw0', 'pwy', 'pwx', 'cr0', 'cry', 'crx'))
pairs(fitmc_yy, pars=c('pw0', 'pwy', 'cr0', 'cry'))


trp <- trans_prob(res$data, thresh = 1)
print(mean(trp$corr))
print(mean(trp$wetf))
