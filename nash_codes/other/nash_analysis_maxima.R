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



# library(tseries)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



outplot = file.path('..', 'nash_output', 'figures')
outdata = file.path('..', 'nash_output', 'outdata')

source("hbev_module.R")    # main functions for data analysis
source("hbev_functions.R") # other functions

iter = 2000
chains = 4
ndraws = iter*chains/2

# load station data
# datadir = file.path('..','data', 'Data_GHCN_Daily', 'extracted_csv')
filenames = list.files(datadir)
  
# filepath = file.path(datadir, 'USW00094728.csv') # NYCP
# filepath = file.path(datadir, 'USC00044232.csv') # HOUSTON
filepath = file.path(datadir, 'USC00311677.csv') # CHAPEL HILL
# filepath = file.path(datadir, 'USC00207812.csv') # MICHIGAN - STAMBAUGH
# filepath = file.path(datadir, 'USC00478027.csv') # WISCONSIN SPOONER AG RES STN 
# state = 'WISCONSIN'

station_id = 'USC00311677' # NORTH CAROLINA
# station_id = 'USC00478027' # WISCONSIN SPOONER
# station_id = 'USC00207812' # MICHIGAN
  
mystation = strsplit(as.character(dfpos$NAME[dfpos$ID == station_id]), ' ')[[1]][1]
filepath = file.path(datadir, sprintf('%s.csv', station_id)) # NYCP


# filenames = list.files(file.path('..','data', 'datasetsM'), pattern = "\\.csv$")
# datadir = file.path('..','data', 'datasetsM')
# filepath = file.path(datadir, 'Asheville.csv')
# filepath = file.path(datadir, 'Albany.csv')
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
res = table_max(df, Nt=Nt)

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

ndf = data.frame( list(xnn = xnn, ynn = ynn, max = res$max, N = res$N, totals = res$totals, 
                       means =res$totals/res$N, sdwets = res$sdwets, skews = res$skews))


ridgepos <- function(xlon, ylat){
  nelem = length(xlon)
  quadrant = rep(0, nelem)
  for (i in 1:nelem){
      x = xlon[i]
      y = ylat[i]
    if (x > 0 & y > 0){
      quadrant[i] = 'NE'
    } else if (x > 0 & y < 0 ) {
        quadrant[i] = 'SE'
    } else if (x < 0 & y > 0){
        quadrant[i] = 'NW'
    } else if (x < 0 & y < 0) {
        quadrant[i] = 'SW'
    }
  }
  return(quadrant)
}

ndf$resquad = ridgepos(xnn, ynn)

resquad = ridgepos(xnn, ynn)
nse = length(resquad[resquad == 'SE'])
nsw = length(resquad[resquad == 'SW'])
nne = length(resquad[resquad == 'NE'])
nnw = length(resquad[resquad == 'NW'])


# p0 <- ggplot(ndf)+
#   geom_point(aes(x=xnn, y=ynn, color = max, size = 1))

p1 <- ggplot(ndf, aes(x=resquad, y=totals, fill = resquad, alpha = 0.6)) +
  # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
  geom_boxplot() + 
  ggtitle('Total prcp.') + 
  geom_jitter(width=0.25, alpha=0.5) 

p2 <- ggplot(ndf, aes(x=resquad, y=N, fill = resquad, alpha = 0.6)) +
  # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
  geom_boxplot() + 
  ggtitle('Number of events') + 
  geom_jitter(width=0.25, alpha=0.5) 

p3 <- ggplot(ndf, aes(x=resquad, y=means, fill = resquad, alpha = 0.6)) +
  # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
  geom_boxplot() + 
  ggtitle('Mean of wet-day prcp') + 
  geom_jitter(width=0.25, alpha=0.5) 

p4 <- ggplot(ndf, aes(x=resquad, y=sdwets, fill = resquad, alpha = 0.6)) +
  # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
  geom_boxplot() + 
  ggtitle('Stdv of wet-day prcp') + 
  geom_jitter(width=0.25, alpha=0.5) 

p5 <- ggplot(ndf, aes(x=resquad, y=sdwets, fill = resquad, alpha = 0.6)) +
  # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
  geom_boxplot() + 
  ggtitle('Skewness of wet-day prcp') + 
  geom_jitter(width=0.25, alpha=0.5) 

p6 <- ggplot(ndf, aes(x=resquad, y=max, fill = resquad, alpha = 0.6)) +
  # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
  geom_boxplot() + 
  ggtitle('Maximum prcp') + 
  geom_jitter(width=0.25, alpha=0.5)  +
  theme_bw() + 
  labs( x= "Ridge position")

fig <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)
ggsave(file.path( outplot, sprintf('example_quandrants_%s.png', mystation)), plot=fig)


Nt = 92
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
nwets_data = list(M = M, Mgen = Mgen, Nt = Nt, N=N, XN = XN, YN = YN, XNgen = XNgen, YNgen = YNgen)


# nwets_data_origrep = list(M = M, Mgen = M, Nt = Nt, N=N, XN = XN, YN = YN, XNgen = XN, YNgen = YN)
# nwets_data_southrep = list(M = M, Mgen = Mgen, Nt = Nt, N=N, XN = XN, YN = YN, XNgen = XNgen2, YNgen = YNgen2)
# nwets_data_northrep = list(M = M, Mgen = Mgen, Nt = Nt, N=N, XN = XN, YN = YN, XNgen = XNgen, YNgen = YNgen)


nwets_depno_fit <- stan(file = 'nash_nj_bin_depno.stan', data = nwets_data, iter = iter)
nwets_depxy_fit <- stan(file = 'nash_nj_bin_depxy.stan', data = nwets_data, iter = iter)
nwets_depyy_fit  <- stan(file = 'nash_nj_bin_depyy.stan',  data = nwets_data, iter = iter)

nsim_depxy <- rstan::extract(nwets_depxy_fit)
nsim_depyy <- rstan::extract(nwets_depyy_fit)
nsim_depno <- rstan::extract(nwets_depno_fit)

pairs(nwets_depxy_fit, pars = c('n0', 'ny','nx',  'pn[1]'))

# nsim$pn
# pnm = colMeans(nsim$pn)
# plot(YN, pnm)
# plot(ndf$ynn, ndf$N, col = 'red')

infc_depxy = info_crit(nsim_depxy$log_lik, psis = FALSE, is_log_like = TRUE)
infc_depyy = info_crit(nsim_depyy$log_lik, psis = FALSE, is_log_like = TRUE)
infc_depno = info_crit(nsim_depno$log_lik, psis = FALSE, is_log_like = TRUE)

# compute information criterion
lpml_depxy = infc_depxy$lpml
lpml_depyy = infc_depyy$lpml
lpml_depno = infc_depno$lpml







# fit model for seasonal precipitation totals
ntots_data = list(M = length(ndf$totals), Mgen = Mgen,  
                  totals=ndf$totals, XN = ndf$xnn, YN = ndf$ynn,
                  XNgen = XNgen, YNgen = YNgen)
ntots_fit_depxy <- stan(file = 'nash_totals_wei_depxy.stan', data = ntots_data, iter = iter)
ntots_fit_depyy <- stan(file = 'nash_totals_wei_depyy.stan',  data = ntots_data, iter = iter)
ntots_fit_depno <- stan(file = 'nash_totals_wei_depno.stan', data = ntots_data, iter = iter)

tsim_depxy <- rstan::extract(ntots_fit_depxy)
tsim_depyy <- rstan::extract(ntots_fit_depyy)
tsim_depno <- rstan::extract(ntots_fit_depno)

pairs(ntots_fit_depxy, pars = c('p0', 'c0', 'cy','cx', 'w'))
pairs(ntots_fit_depyy, pars = c('p0', 'c0', 'cy', 'w'))
pairs(ntots_fit_depno, pars = c('p0','C', 'w'))



thresh = 0
totnwets = length(res$data[res$data > thresh])
nash_data = list(M = length(ndf$means),
                 Mgen = Mgen,
                  Nt = Nt,
                  y = res$data,
                  N = ndf$N,
                  XN = ndf$xnn, YN = ndf$ynn,
                  XNgen = XNgen, YNgen = YNgen,
                  totnwets = totnwets)

nash_fit_depxy <- stan(file = 'nash_xij_wei_depxy.stan', data = nash_data, iter = iter)
nash_fit_depyy <- stan(file = 'nash_xij_wei_depyy.stan', data = nash_data, iter = iter)
nash_fit_depno <- stan(file = 'nash_xij_wei_depno.stan', data = nash_data, iter = iter)

nashsim_depxy <- rstan::extract(nash_fit_depxy)
nashsim_depyy <- rstan::extract(nash_fit_depyy)
nashsim_depno <- rstan::extract(nash_fit_depno)

hist(nashsim_depxy$Nrep)
hist(nashsim_depno$Nrep)


nash_data_origrep = list(M = length(ndf$means),
                 Mgen = M,
                  Nt = Nt,
                  y = res$data,
                  N = ndf$N,
                  XN = ndf$xnn, YN = ndf$ynn,
                  XNgen = XN, YNgen = YN,
                  totnwets = totnwets)

nash_data_northrep = list(M = length(ndf$means),
                 Mgen = Mgen,
                  Nt = Nt,
                  y = res$data,
                  N = ndf$N,
                  XN = ndf$xnn, YN = ndf$ynn,
                  XNgen = XNgenN, YNgen = YNgenN,
                  totnwets = totnwets)

nash_data_southrep = list(M = length(ndf$means),
                 Mgen = Mgen,
                  Nt = Nt,
                  y = res$data,
                  N = ndf$N,
                  XN = ndf$xnn, YN = ndf$ynn,
                  XNgen = XNgenS, YNgen = YNgenS,
                  totnwets = totnwets)

nwets_depxy_fit_origrep <- stan(file = 'nash_xij_wei_depxy.stan',  data = nash_data_origrep, iter = iter)
nwets_depxy_fit_southrep <- stan(file = 'nash_xij_wei_depxy.stan', data = nash_data_southrep, iter = iter)
nwets_depxy_fit_northrep <- stan(file = 'nash_xij_wei_depxy.stan', data = nash_data_northrep, iter = iter)

nwets_sim_origrep = rstan::extract(nwets_depxy_fit_origrep)
nwets_sim_northrep = rstan::extract(nwets_depxy_fit_northrep)
nwets_sim_southrep = rstan::extract(nwets_depxy_fit_southrep)

plot(density(nwets_sim_origrep$maxrep), log = 'y')
lines(density(nwets_sim_southrep$maxrep), col = 'red')
lines(density(nwets_sim_northrep$maxrep), col = 'green')

dim(nwets_sim_origrep$maxrep)

# Weibull plotting position for annual maxima non exceedance frequency
FiM = 1:M/(1+M)
TrM = 1/(1-FiM)
FiMgen = 1:Mgen/(1+Mgen)
TrMgen = 1/(1-FiMgen)

XIor = apply(nwets_sim_origrep$maxrep, 1, sort)
XIno = apply(nwets_sim_northrep$maxrep, 1, sort)
XIso = apply(nwets_sim_southrep$maxrep, 1, sort)

XIor_mean = apply(XIor, 1, mean)
XIno_mean = apply(XIno, 1, mean)
XIso_mean = apply(XIso, 1, mean)
XIor_std = apply(XIor, 1, sd)
XIno_std = apply(XIno, 1, sd)
XIso_std = apply(XIso, 1, sd)

plot(TrM, res$Xi, log = 'x')
lines(TrM, XIor_mean)
polygon(c(TrM, rev(TrM)), c(XIor_mean + XIor_std, rev(XIor_mean - XIor_std)), border = NA, col = rgb(0, 0, 0, 0.3))
lines(TrMgen, XIno_mean, col ='blue')
lines(TrMgen, XIso_mean, col ='red')
polygon(c(TrMgen, rev(TrMgen)), c(XIno_mean + XIno_std, rev(XIno_mean - XIno_std)), border = NA, col = rgb(0.0, 0, 1, 0.3))
polygon(c(TrMgen, rev(TrMgen)), c(XIso_mean + XIso_std, rev(XIso_mean - XIso_std)), border = NA, col = rgb(1, 0, 0, 0.3))
abline(v = 20, col="green", lwd=2, lty=2)

# compute spread of the two cases over std for Tr = 50
# compute position of the return time closest to my target
target_tr = 20
postr = which.min(abs(TrMgen - target_tr))

sq50 = XIor_std[postr]
dq50 = XIno_mean[postr] - XIso_mean[postr] 
alpha50 = dq50/sq50 


color_scheme_set("orange")

ppc_intervals(
  # y = mtcars$mpg,
  y = ndf$N,
  # yrep = posterior_predict(nsim),
  yrep = nwets_sim_origrep$Nrep,
  x = ndf$ynn,
  prob = 0.5
) +
  labs(
    x = "(Normalized) ridge latitude",
    y = "Number of events N on JJA",
    title = "50% posterior predictive intervals for N ",
    subtitle = "\nvs ridge meridional position"
  ) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")



# plot(XIor[,1])


# generate posterior predictive distributions based on the parameters of the N-distribution



pairs(nash_fit_depxy, pars = c('w0', 'wy', 'wx', 'c0', 'cy', 'cx', 'n0','ny', 'nx'))

# 
# plot(nash_fit_depxy, pars = c('wy', 'cy','ny'))
# plot(nash_fit_depyy, pars = c('wy', 'cy','ny'))
# 
# plot(density(nashsim_depxy$Crep))
# lines(density(nashsim_depyy$Crep), col = 'green')
# lines(density(nashsim_depno$C), col = 'purple')
# lines(density(nashsim_depxy$C), col = 'red')
# lines(density(nashsim_depy$C), col = 'blue')
# 
# plot(density(nashsim_depxy$wrep))
# lines(density(nashsim_depy$wrep), col = 'purple')
# lines(density(nashsim_nodep$wrep), col = 'green')
# lines(density(nashsim_depxy$w), col = 'red')
# lines(density(nashsim_depy$w), col = 'blue')

# A = nashsim_depxy$xijrep
# A = A[A>0]
# B = nashsim_nodep$xijrep
# B = B[B>0]
# C = res$data[res$data > 0]
# 
# plot(density(B), col = 'red', log = 'xy')
# lines(density(C))
# # plot(density(A))
# lines(density(A), col = 'blue')

# mean(nashsim_depxy$Crep)
# sd(nashsim_depxy$Crep)
# mean(nashsim_depxy$C)
# mean(nashsim_depxy$wrep)
# sd(nashsim_depxy$wrep)
# mean(nashsim_depxy$w)
# 
# # 
# plot(density(log(nashsim_depxy$xijrep)))
# lines(density(log(nashsim_nodep$xijrep)), col = 'red')
# 
# dim(nashsim_depxy$xijrep)
# 
# plot(density(nashsim_depxy$maxrep), log = 'y')
# lines(density(nashsim_depy$maxrep), col = 'red')
# lines(density(nashsim_nodep$maxrep), col = 'green')

# Fi = 1:(Mgen*4000)/(1+Mgen*4000)
# xiy = sort(nashsim_depy$maxrep)
# xino = sort(nashsim_nodep$maxrep)
# xixy = sort(nashsim_depxy$maxrep)
# 
# plot(xino, 1-Fi, log = 'x', type = 'l', lty=1)
# lines(xiy, 1-Fi, col = 'red')
# lines(xixy, 1-Fi, col = 'blue')

infc_depxy = info_crit(nashsim_depxy$log_lik, psis = FALSE, is_log_like = TRUE)
infc_depyy = info_crit(nashsim_depyy$log_lik, psis = FALSE, is_log_like = TRUE)
infc_depno = info_crit(nashsim_depno$log_lik, psis = FALSE, is_log_like = TRUE)

# lppd_depxy = infc_depxy$lppd
lpml_depxy = infc_depxy$lpml
# effnumpar_depxy = (lpml_depxy - lppd_depxy)*totnwets # times sample size
# 
# lppd_depy = infc_depy$lppd
lpml_depyy = infc_depyy$lpml
# effnumpar_depy = (lpml_depy - lppd_depy)*totnwets # times sample size

# lppd_nodep = infc_nodep$lppd
lpml_depno = infc_depno$lpml
# effnumpar_nodep = (lpml_nodep - lppd_nodep)*totnwets # times sample size

# infcnj_depxy = info_crit(nashsim_depxy$log_lik_nj, psis = FALSE, is_log_like = TRUE)
# infcnj_depyy = info_crit(nashsim_depy$log_lik_nj, psis = FALSE, is_log_like = TRUE)
# infcnj_nodep = info_crit(nashsim_nodep$log_lik_nj, psis = FALSE, is_log_like = TRUE)
# lpml_depxy_nj = infcnj_depxy$lpml
# lpml_depy_nj = infcnj_depy$lpml
# lpml_nodep_nj = infcnj_nodep$lpml

  # ppc_intervals(
  #   # y = mtcars$mpg,
  #   y = ndf$N,
  #   # yrep = posterior_predict(nsim),
  #   yrep = nsim$Nrep,
  #   x = ndf$ynn,
  #   prob = 0.5
  # ) +
  #   labs(
  #     x = "(Normalized) ridge latitude",
  #     y = "Number of events N on JJA",
  #     title = "50% posterior predictive intervals for N ",
  #     subtitle = "\nvs ridge meridional position"
  #   ) +
  #   panel_bg(fill = "gray95", color = NA) +
  #   grid_lines(color = "white")


