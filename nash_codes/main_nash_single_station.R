# exploratory analysis of NASH data


rm(list=ls()) 
library(reshape2)
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
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

datadir = file.path('..', 'data', 'Data_GHCN_Daily', 'extracted_csv')
filenames = list.files(datadir)
dfpos = read.csv(file.path(datadir, '..',  'list_ushcn_stats.csv'))


# dflen = read.csv(file.path('..', 'nash_output', 'length_file.csv'),sep = '\t')
# mylen = subset(dflen, dflen$)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
dfpos$STATE = trim(dfpos$STATE)
dfpos$NAME = trim(dfpos$NAME)


# dflakes <- subset(dfpos, dfpos$LATITUDE > 45 & dfpos$LONGITUDE > -94)
dflakes <- subset(dfpos, dfpos$STATE == 'MN')
dfgeorgia <- subset(dfpos, dfpos$STATE == 'GA')
# dften <- subset(dfpos, dfpos$STATE == 'TN')


outplot1 = file.path('..', 'nash_output_figures')
outplot = file.path(outplot1, 'stats_figures')
dir.create(outplot, showWarnings = FALSE)
# outdata = file.path('..', 'nash_output', 'outdata')
# dir.create(outdata, showWarnings = FALSE)

source("nashfun.R")    # main functions for data analysis



#######################  parameters to set for MCMC  ###########################
Nt = 92 # (30 + 31 + 31 expected observations)
cores = 4 # let us try to set this for both PC and cluster
min_num_years = 20 # if less skip station (years with enough data & overlap with NASH data)
maxmiss_per_season = 4 # remove all years with more than 4 events /season missing
thresh = 1 # for ordinary events
iter = 2000
chains = 4
ndraws = iter*chains/2

  
  
# station_id = 'USW00094728' # NEW YORK CENTRAL PARK
# station_id = 'USC00044232' # CALIFORNIA
# station_id = 'USC00218419' # Minnesota, Two Harbors
# station_id = 'USC00311677' # NORTH CAROLINA, CHAPEL HILL
# station_id = 'USC00478027' # WISCONSIN, SPOONER
# station_id = 'USC00475516' # WISCONSIN, MINOQUA
# station_id = 'USC00207812' # MICHIGAN, STAMBAUGH
# station_id = 'USC00409502' # TENNESSEE, WAYNESBORO
# station_id = 'USC00093621' # Georgia, Gainesville
# station_id = "USC00229079" # stat with larger decrease in maxima going NW (MS)
# station_id = "USC00221962" # stat with second larger decrease in maxima going NW
# station_id = "USC00478110" # stat with larger INCREASE in maxima going NW (WI)
# station_id = "USC00041048" # stat with second larger INCREASE in maxima going NW ( BRAWLEY CA)
# station_id = "USC00475808" # stat with  -a- larger INCREASE in maxima going NW 
################################################################################  





#######################comment for running all the stations at once #########
STATID = c(
           'USW00094728',
           'USC00044232',
           # 'USC00218419',
           # 'USC00478027',
           'USC00475516')
           # 'USC00207812',
           # 'USC00409502',
           # 'USC00093621')
           # "USC00229079",
           # "USC00221962",
           # "USC00041048")
           # "USC00478110", # outlier in Wisconsin
           # "USC00311677", # CHAPEL HILL, NC
           # "USC00475808")
NSTATS = length(STATID)

for (iii in 1:NSTATS){
  station_id = STATID[iii]
################################################################################


# station_state = substr(dfpos$STATE[dfpos$ID==station_id], 1, 2)
  print( sprintf('reading the data for station %s', iii))
station_state = dfpos$STATE[dfpos$ID==station_id]
station_name = dfpos$NAME[dfpos$ID==station_id]

# ddx = subset(dfpos,dfpos$ID==station_id)

# read precipitation and nash data
filepath = file.path(datadir, sprintf('%s.csv', station_id)) 
df <- read_summer_data(filepath, nmaxmiss = maxmiss_per_season, Nt_JJA = Nt)
prcp_years = unique(df$YEAR)
dfnash <- read_nash_values(prcp_years = prcp_years)
xnn = dfnash$XNN
ynn = dfnash$YNN
dfnash$resquad = ridgepos(xnn=xnn, ynn=ynn) # compute quadrants



# compute annual maxima statistics

  print( sprintf('compute and plot annual maxima statistics for station %s', iii))

res = table_max(df, Nt=Nt, thresh = thresh)


source("nashfun.R")    # main functions for data analysis
quadplot <- plot_quadrant_stats(res, xnn, ynn, station_id, station_name, station_state)
  ggsave(file.path( outplot, sprintf('example_quadrants_%s_%s.png', station_name, station_state)),
         width = 6, height = 6.5, plot=quadplot)
  

# RIDGE position to simulate data from:
# Mgen = 60 # number of years of data to generate
Mgen = length(xnn)
XNgen = xnn
YNgen = ynn
# XNgenN = xnn
# YNgenN = ynn
# XNgenS = xnn
# YNgenS = ynn
XNgenS = rep(-1.0, Mgen) # west
YNgenS = rep(-1.0, Mgen) # south
XNgenN = rep(-1.0, Mgen) # west
YNgenN = rep(1.0, Mgen) # north


totnwets = length(res$data[res$data > 0]) #these are already excesses
nashdata = list(M = res$nyears, Mgen = Mgen, Nt = Nt, N=res$N, y = res$data,
                totals = res$totals,
                totnwets = totnwets,
                XN = dfnash$XNN, YN = dfnash$YNN, 
                XNgen = XNgen, YNgen = YNgen,
                XNgenS = XNgenS, YNgenS = YNgenS,
                XNgenN = XNgenN, YNgenN = YNgenN,
                NP3=2, NP5=4
                )

  print( sprintf('fit Stan model for station %s', iii))

fit  <- stan(file = 'nash_xij_weiNSO_xydep.stan', data = nashdata, iter = iter, chains = chains)
# fitmc  <- stan(file = 'nash_xij_weiNSOmarkov_xydep.stan', data = nashdata, iter = iter, chains = chains)

sim <- rstan::extract(fit)
# simmc <- rstan::extract(fitmc)

# pairs(fit, pars = c('n0', 'ny','nx'))
# pairs(fit, pars = c('c0', 'cy','cx'))
# pairs(fit, pars = c('pw0', 'pwy','pwx'))
# pairs(fit, pars = c('cr0', 'cry','crx'))

infc_depxy = info_crit(sim$log_lik, psis = FALSE, is_log_like = TRUE)
lpml_depxy = infc_depxy$lpml

# if maxima are simulated in 2 scenarios, compute the frequency of extremes:

      # ## To CHECK ##
      # XIso = apply(sim$maxrepS, 1, sort) + thresh
      # XIno = apply(sim$maxrepN, 1, sort) + thresh
      # XIso_mean = apply(XIso, 1, mean)
      # XIno_mean = apply(XIno, 1, mean)
      # XIso_std = apply(XIso, 1, sd)
      # XIno_std = apply(XIno, 1, sd)
      # postr20 = which.min(abs(TrMgen - target_tr_20))
      # sqTR20 = XIno_std[postr20]
      # dqTR20 = XIno_mean[postr20] - XIso_mean[postr20] 
      # 
      # postr50 = which.min(abs(TrMgen - target_tr_50))
      # sqTR50 = XIno_std[postr50]
      # dqTR50 = XIno_mean[postr50] - XIso_mean[postr50] 
      # 
      # qDq20 = dqTR20/sqTR20 
      # qDq50 = dqTR50/sqTR50 
      # 
      # qnorm20 = dqTR20/XIso_mean[postr20] 
      # qnorm50 = dqTR50/XIso_mean[postr50] 
      ###############

# dim(sim$maxrep)


  print( sprintf('compute return time quantiles for station %s', iii))

# Weibull plotting position for annual maxima non exceedance frequency
M = res$nyears
FiM = 1:M/(1+M)
TrM = 1/(1-FiM)
FiMgen = 1:Mgen/(1+Mgen)
TrMgen = 1/(1-FiMgen)
XIor = apply(sim$maxrepO, 1, sort) + thresh
XIno = apply(sim$maxrepN, 1, sort) + thresh
XIso = apply(sim$maxrepS, 1, sort) + thresh
XIor_mean = apply(XIor, 1, mean)
XIno_mean = apply(XIno, 1, mean)
XIso_mean = apply(XIso, 1, mean)
XIor_std = apply(XIor, 1, sd)
XIno_std = apply(XIno, 1, sd)
XIso_std = apply(XIso, 1, sd)

# compute spread of the two cases over std for Tr = 50
# compute position of the return time closest to my target
target_tr = 20
postr = which.min(abs(TrMgen - target_tr))
sq50 = XIno_std[postr]
dq50 = XIno_mean[postr] - XIso_mean[postr] 
qDq = dq50/sq50 
qnorm50 = (XIno_mean[postr] - XIso_mean[postr] )/XIso_mean[postr] 

# to plot pdf of daily data:
xij_obs = as.vector(res$data)
mydimz = dim(sim$xijrepO)
# xij_rep = matrix(0, nrow = mydimz[1] , ncol = mydimz[2]*mydimz[3] )

xij_rep <- matrix(sim$xijrepO, nrow = mydimz[1] , ncol = mydimz[2]*mydimz[3], byrow = 'False' )

# sum(xij_rep2[17, ])
# sum(sim$xijrepO[17, ,])

# xij_rep  = reshape(sim$xijrep0, nrow = mydimz[1] , ncol = mydimz[2]*mydimz[3] )
# library("pracma")
# xij_rep  = Reshape(sim$xijrep0, mydimz[1] , mydimz[2]*mydimz[3] )

# for (ii in 1:mydimz[1]){
#   for (jj in 1:mydimz[2]){
#     for (kk in 1:mydimz[3]){
#       # xij_rep[((jj-1)*mydimz[3] + kk), ii] <- sim$xijrepO[ii, jj, kk]  
#       xij_rep[ii, ((kk-1)*mydimz[2] + jj)] <- sim$xijrepO[ii, jj, kk]  
#     }
#   }
# }

# xij_rep2 <- array_reshape(sim$xijrepO, c(mydimz[1], mydimz[2]*mydimz[3]), order = "C")




# # comparing extremes: Markov Chain model vs Binomial model
# XIormc =apply(simmc$maxrepO, 1, sort) + thresh
# XIormc_mean =apply(XIormc, 1, mean)
# plot(TrM, XIormc_mean)
# lines(TrM, XIor_mean, col = 'red')

dfquantemp <- data.frame(TrM = TrM, Xi = res$Xi)
dfquant <- data.frame(TrMgen = TrMgen,
                      XIor_mean = XIor_mean,
                      XIor_std = XIor_std,
                      XIno_mean = XIno_mean,
                      XIno_std =   XIno_std,
                      XIso_mean = XIso_mean,
                      XIso_std =   XIso_std
                      )


  print( sprintf('plot results for station %s', iii))

ppcfontsize = 12
qfignam = file.path(outplot, sprintf('ppchecks_quantiles_%s_%s.png', station_name, station_state))
qplot2<- ggplot(dfquant)+
  # geom_ribbon(aes(x = TrMgen, ymin = XIor_mean - XIor_std, ymax = XIor_mean + XIor_std), fill = "grey70", alpha = 0.6) +
  # geom_line(aes(x = TrMgen, y = XIor_mean), color="black", size=1.2)+
    geom_ribbon(aes(x = TrMgen, ymin = XIno_mean - XIno_std, ymax = XIno_mean + XIno_std), fill = "darkred", alpha = 0.5) +
  geom_line(aes(x = TrMgen, y = XIno_mean, color='darkred'), size=1.2)+
    geom_ribbon(aes(x = TrMgen, ymin = XIso_mean - XIso_std, ymax = XIso_mean + XIso_std), fill = "darkblue", alpha = 0.5) +
  geom_line(aes(x = TrMgen, y = XIso_mean, color='darkblue'), size=1.2)+
  geom_point(data=dfquantemp, aes(x=TrM, y=Xi), size=2.8)+
  theme_bw()+
  scale_color_identity(name = "Quadrant",
                       breaks = c("darkred", "darkblue"),
                       labels = c("North-West", "South-West"),
                       guide = "legend")+
  ylab("Daily rainfall accumulation [mm]")+
  xlab("Return time [years]")+
 scale_x_continuous(trans = 'log10')+
 scale_y_continuous(trans = 'log10')+
  theme(text = element_text(size=ppcfontsize),
            axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      legend.position = c(0.7, 0.3))

qplot1 <- ggplot(dfquant) +
  geom_ribbon(aes(x = TrMgen, ymin = XIor_mean - XIor_std, ymax = XIor_mean + XIor_std), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = TrMgen, y = XIor_mean, color="black"), size=1.2)+
  geom_point(data=dfquantemp, aes(x=TrM, y=Xi, shape='Observed maxima'), size= 2.8) +
  #   geom_ribbon(aes(x = TrMgen, ymin = XIno_mean - XIno_std, ymax = XIno_mean + XIno_std), fill = "darkred", alpha = 0.6) +
  # geom_line(aes(x = TrMgen, y = XIno_mean), color='darkred', size=1.2)+
  #   geom_ribbon(aes(x = TrMgen, ymin = XIso_mean - XIso_std, ymax = XIso_mean + XIso_std), fill = "darkblue", alpha = 0.6) +
  # geom_line(aes(x = TrMgen, y = XIso_mean), color='darkblue', size=1.2)+
  theme_bw()+
    scale_color_identity(name = "",
                       breaks = c("black"),
                       labels = c("Replicates"),
                       guide = "legend")+
      # scale_shape_identity(name = "",
      #                  breaks = c("black"),
      #                  labels = c("Replicates"),
      #                  guide = "legend")+
  ylab("Daily rainfall accumulation [mm]")+
  xlab("Return time [years]")+
 scale_x_continuous(trans = 'log10')+
 scale_y_continuous(trans = 'log10')+
  theme(text = element_text(size=ppcfontsize),
            axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      legend.title=element_blank(),
      legend.position = c(0.7, 0.3))
# qplot = grid.arrange(qplot1, qplot2, ncol=2)

# ggsave(qfignam, width = 8, height = 4, dpi = 300, plot=qplot)
# plot(qplot)

# plot frequency of annual maxima in the two simulated scenarios
# png(file=qfignam, width = 4, height = 4)
# qfig <- plot(TrM, res$Xi, log = 'xy')
# lines(TrM, XIor_mean)
# polygon(c(TrM, rev(TrM)), c(XIor_mean + XIor_std, rev(XIor_mean - XIor_std)), 
#                       border = NA, col = rgb(0, 0, 0, 0.3))
# lines(TrMgen, XIno_mean, col ='blue', lwd = 2)
# lines(TrMgen, XIso_mean, col ='red', lwd = 2)
# polygon(c(TrMgen, rev(TrMgen)), c(XIno_mean + XIno_std, 
#         rev(XIno_mean - XIno_std)), border = NA, col = rgb(0.0, 0, 1, 0.3))
# polygon(c(TrMgen, rev(TrMgen)), c(XIso_mean + XIso_std, 
#         rev(XIso_mean - XIso_std)), border = NA, col = rgb(1, 0, 0, 0.3))
# dev.off()

# jpeg(file="saving_plot1.jpeg")
# png(file=, width = 4, height = 4)

# abline(v = 20, col="green", lwd=2, lty=2)

# if generatring 2D matrix (as in the case of the Markov model), flatten in out
# xijrepO <- t(matrix(simmc$xijrepO, prod(dim(simmc$xijrepO)[2:3]), dim(simmc$xijrepO)[1]))
# print(dim(xijrepO))


# ppcfontsize=1
# numrep = 20
# MAXO = melt(data.frame(t(sim$maxrepO[1:numrep,])))
# MAXN = melt(data.frame(t(sim$maxrepN[1:numrep,])))
# MAXS = melt(data.frame(t(sim$maxrepS[1:numrep,])))
# ppm1 <- ppc_dens_overlay(y = res$Xi,
#                  yrep = sim$maxrepO[1:numrep, ] + thresh
# ) +
# theme(text = element_text(size=ppcfontsize),
#             axis.text=element_text(size=12),
#       axis.title=element_text(size=14),
#       legend.position = "None")+
#   # scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))+
#   scale_colour_manual(values=c("black", "skyblue1"), labels = c("Replicates", "Sample"))+
#   ylab("probability density")+
#   xlab("Seasonal maxima")+
#   geom_density(data=MAXN, aes(x=value, group=variable), color='red', alpha=0.2)+
#   geom_density(data=MAXS, aes(x=value, group=variable), color='blue')

#   ppcfontsize=1
#   numrep = 20
# XIJO = melt(data.frame(t(sim$xijrepO[1:numrep,])))
# XIJN = melt(data.frame(t(sim$xijrepN[1:numrep,])))
# XIJS = melt(data.frame(t(sim$xijrepS[1:numrep,])))
# 
# ppm1 <- ppc_dens_overlay(y = as.vector(res$data),
#                  yrep = sim$xijrepO[1:numrep, ]
# ) +
# theme(text = element_text(size=ppcfontsize),
#             axis.text=element_text(size=12),
#       axis.title=element_text(size=14),
#       legend.position = "None")+
#   # scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))+
#   scale_colour_manual(values=c("black", "green"), labels = c("Replicates", "Sample"))+
#   ylab("probability density")+
#   xlab("xij")+
#   scale_x_continuous(trans = 'log10',
#                     breaks = c(0.1, 1, 10, 100),
#                    limits = c(0.01, max(res$data)+100)
#                     ) +
#   geom_density(data=XIJN, aes(x=value, group=variable), color='red', alpha=0.2)+
#   geom_density(data=XIJS, aes(x=value, group=variable), color='blue', alpha = 0.2)


  ppcfontsize=1
  numrep = 100
ppc3 <- ppc_dens_overlay(y = res$Xi,
                 yrep = sim$maxrepO[1:numrep, ] + thresh
) +
theme(text = element_text(size=ppcfontsize),
            axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      legend.position = "None")+
  # scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))+
  scale_colour_manual(values=c("black", "skyblue1"), labels = c("Replicates", "Sample"))+
  # xlim( c(pmax(min(res$Xi)-30, 0), max(res$Xi))+50)+
  ylab("Probability density")+
  xlab("Seasonal maxima [mm]")

# ppc1 <- ppc_dens_overlay(y = res$N,
#                  yrep = sim$NrepO[1:numrep, ]
# ) +
# ppc1 <- ppc_bars(y = res$N,
#                  yrep = sim$NrepO[1:numrep, ],
#                  size=0.5,
#                  prob = 0.90,
#                  fatten=1
# mybreaks = 11:41 + 0.5
simvalues = sim$NrepO[1:numrep,]
maxval = max(max(simvalues), max(res$N))
minval = min(min(simvalues), min(res$N))
mybreaks = minval:maxval
ncounts = length(mybreaks)-1
ndraws = dim(sim$NrepO)[1]
mycounts = matrix(0, nrow = ncounts, ncol = numrep)
myhist_obs = hist(res$N, breaks = mybreaks, freq = TRUE)
for (ix in 1:numrep){
  print(ix)
  myhist_ix = hist(sim$NrepO[ix, ], breaks = mybreaks, freq = TRUE)
  mycounts[,ix ] = myhist_ix$counts
}

mean_counts = rowMeans(mycounts)
upper_quants = apply(mycounts, 1, quantile, probs = 0.75)
lower_quants = apply(mycounts, 1, quantile, probs = 0.25)
# compute median and 90% intervals

# myhist0 = hist(res$N, breaks = mybreaks, freq = TRUE)
# points(myhist_ix$mids, mean_counts)
# points(myhist_ix$mids, upper_quants)
# points(myhist_ix$mids, lower_quants)
# hist(sim$NrepO[1:numrep, 1], breaks = 19)
#  ppc1 <- ppc_intervals(y = res$N,
#                  yrep = sim$NrepO[1:numrep, ],
#                  size=0.5,
#                  prob = 0.50
#                  # prob = 0.0,
# )+
    ppc1 <- bayesplot::ppc_intervals(y = myhist_obs$counts,
                 yrep = t(mycounts),
                 size=1,
                 prob_outer = 0.90,
                 prob = 0.50
                 # prob = 0.0,
)+
    # bayesplot:::scale_color_ppc_dist(values=c("black", "skyblue1"),labels = c("Sample", "Replicates")) +
  # bayesplot:::scale_fill_ppc_dist(values=c("black", "skyblue1"),labels = c("Sample", "Replicates"))+
      # bayesplot:::scale_fill_ppc_dist(values=c("skyblue1"),labels = c("Sample")) +
      # bayesplot:::scale_color_ppc_dist(values=c("grey20"),labels = c("Replicates")) +
        bayesplot:::scale_color_ppc_dist(values = c("grey20", "skyblue1"), labels = c("Sample", "Replicates")) +
  bayesplot:::scale_fill_ppc_dist(values = c("grey20", "skyblue1"), labels = c("Sample", "Replicates"))+
  # bayesplot:::scale_fill_ppc_dist(values=c("", "skyblue1"),labels = c("Sample", "Replicates"))+
  xlab("Number of events")+
  ylab("Count")+
  # scale_colour_manual(values=c("black", "skyblue1"), labels = c("Replicates", "Sample"))+
  # scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))+
theme( 
              axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      # legend.position = "None")
       legend.position = c(0.85, 0.78))
    
      # theme(
      #       axis.text=element_text(size=12),
      # axis.title=element_text(size=14))+
  # bayesplot:::scale_color_ppc_dist(labels = c("Sample", "Replicates")) +
  # bayesplot:::scale_fill_ppc_dist(labels = c("Sample", "Replicates"))
    # ppc1
      # legend.position = c(0.2, 0.8))
# ppc2 <- ppc_dens_overlay(y = as.vector(res$data),
#                  yrep = sim$xijrepO[1:numrep, ]
                 ppc2 <- ppc_dens_overlay(y = xij_obs,
                 yrep = xij_rep
) + 
# ppc3 <- ppc_dens_overlay(y = as.vector(res$data),     # for MARKOV
#                  yrep = xijrepO[1:numrep, ]
# ) + 
  xlab("Rainfall accumulations [mm]")+
  ylab("Probability density")+
theme(text = element_text(size=ppcfontsize),
      axis.text=element_text(size=12),
      axis.title=element_text(size=14),
       legend.position = c(0.3, 0.78))+
       # legend.position = "None")+
scale_x_continuous(trans = 'log10',
                    breaks = c(0.1, 1, 10, 100),
                   limits = c(0.01, max(res$data)+100)
                    ) +
  scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))
# ppc123 <- grid.arrange(ppc1, ppc2, ppc3, ncol=3)
# ggsave(file.path(outplot, sprintf('ppchecks_pdfs_%s.png', station_state)),
#         width = 9, height = 3, dpi = 300, plot=ppc123)


ppc6 <- ppc_ecdf_overlay(y = res$Xi,
                 yrep = sim$maxrepO[1:numrep, ] + thresh
) +
theme(text = element_text(size=ppcfontsize),
            axis.text=element_text(size=12),
      axis.title=element_text(size=14),
                  # axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
      legend.position = "None")+
  scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))+
  ylab("Cumulative probability")+
  xlab("Seasonal maxima [mm]")
ppc4 <- ppc_ecdf_overlay(y = res$N,
                 yrep = sim$NrepO[1:numrep, ]
)+
  xlab("Number of events")+
  ylab("Cumulative probability")+
  scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))+
theme(text = element_text(size=ppcfontsize), 
            axis.text=element_text(size=12),
      axis.title=element_text(size=14),
                  # axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
      legend.position = "None")
# ppc6 <- ppc_ecdf_overlay(y = as.vector(res$data), # for MARKOV
#                  yrep = xijrepO[1:numrep, ]
# ) + 
  # ppc5 <- ppc_ecdf_overlay(y = as.vector(res$data),
  #                yrep = sim$xijrepO[1:numrep, ]
                   ppc5 <- ppc_ecdf_overlay(y = xij_obs,
                 yrep = xij_rep
) + 
  xlab("Rainfall accumulations [mm]")+
  ylab("Cumulative probability")+
theme(text = element_text(size=ppcfontsize),
      axis.text=element_text(size=12),
      axis.title=element_text(size=14),
            # axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
       # legend.position = c(0.3, 0.8))+
       legend.position = "None")+
       # legend.position = "None")+
scale_x_continuous(trans = 'log10',
                    breaks = c(0.1, 1, 10, 100),
                   limits = c(0.01, max(res$data)+100)
                    ) +
  scale_colour_manual(values=c("black", "skyblue1"), labels = c("Sample", "Replicates"))
# ppc456 <- grid.arrange(ppc4, ppc5, ppc6, ncol=3)
# ggsave(file.path(outplot, sprintf('ppchecks_cdfs_%s.png', station_id)),
#         width = 9, height = 3, dpi = 300, plot=ppc456)

# ppc123456 <- grid.arrange(ppc1, ppc2, ppc3, ppc4, ppc5, ppc6, ncol=3)
# ggsave(file.path(outplot, sprintf('ppchecks_pdfscdfs_%s_%s.png', station_name, station_state)),
#         width = 9, height = 6, dpi = 300, plot=ppc123456)


# plot(density(res$Xi))
# lines(density(XIno_mean))

# plot posterior predictive distribution for number of events (Binomial)


# BINOMIAL
ppn1 <- ppc_intervals(
  # y = mtcars$mpg,
  y = res$N,
  # yrep = posterior_predict(nsim),
  # yrep = nwets_sim_origrep$Nrep,
  # yrep = simmc$NrepO,
  yrep = sim$NrepO[1:numrep, ],
  x = dfnash$LAT,
  prob = 0.5,
  prob_outer = 0.90
) +
  theme_bw()+
  legend_move(c(0.8, 0.8))+
  labs(
    x = "NASH Ridge latitude [North]",
    y = "Number of events in JJA"
    # title = "50% posterior predictive intervals for nj ",
    # subtitle = "\nvs ridge meridional position"
  ) +
  # panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")+

  theme(
            axis.text=element_text(size=12),
      axis.title=element_text(size=14))+
  bayesplot:::scale_color_ppc_dist(labels = c("Sample", "Replicates")) +
  bayesplot:::scale_fill_ppc_dist(labels = c("Sample", "Replicates"))

# ggsave(file.path(outplot, sprintf('ppchecks_nj_binomial_%s_%s.png', station_name, station_state)),
#         width = 5, height = 5, dpi = 300, plot=ppn1)

ppc11 <- ppc1 + annotate(geom = 'text', label = '(a)', x = -Inf, y = Inf,    hjust = -1, vjust = 1, size=8)
ppc22 <- ppc2 + annotate(geom = 'text', label = '(b)', x =0.03, y = Inf,      hjust = 0, vjust = 1, size=8)
ppc33 <- ppc3 + annotate(geom = 'text', label = '(c)', x = -Inf, y = Inf,    hjust = -1, vjust = 1, size=8)
ppc44 <- ppc4 + annotate(geom = 'text', label = '(d)', x = -Inf, y = Inf,    hjust = -1, vjust = 2, size=8)
ppc55 <- ppc5 + annotate(geom = 'text', label = '(e)', x = 0.02, y = Inf,     hjust = 0, vjust = 2, size=8)
ppc66 <- ppc6 + annotate(geom = 'text', label = '(f)', x = -Inf, y = Inf,    hjust = -1, vjust = 2, size=8)
ppn11 <- ppn1 + annotate(geom = 'text', label = '(g)', x = 26.5, y = Inf,    hjust = 0, vjust = 1.4, size=8)
qplot11 <- qplot1 + annotate(geom = 'text', label = '(h)', x = 1.2, y = Inf, hjust = 0, vjust = 1.4, size=8)
qplot22 <- qplot2 + annotate(geom = 'text', label = '(i)', x = 1.2, y = Inf, hjust = 0, vjust = 1.4, size=8)
# ppc11 <- ppc1+ annotate(geom = 'text', label = '(a)', x = Inf, y = Inf, hjust = 1, vjust = 1, size=8)
# plot(ppc11)
comp_stat_plot = grid.arrange(ppc11, ppc22, ppc33, ppc44, ppc55, ppc66, ppn11, qplot11, qplot22, 
                              ncol=3, top = sprintf("Station %s %s (%s)", station_id, station_name, station_state))

ggsave(file.path(outplot, sprintf('complete_station_%s_%s.png', station_name, station_state)),
        width = 12, height = 12, dpi = 300, plot=comp_stat_plot)


nquantsplot <-  grid.arrange(ppn1, qplot1, qplot2, 
                              ncol=3, top = sprintf("Station %s %s (%s)", station_id, station_name, station_state))

ggsave(file.path(outplot, sprintf('n_quants_%s_%s.png', station_name, station_state)),
        width = 12, height = 4, dpi = 300, plot=nquantsplot)



# remove fit and sim from the memory
rm(fit, sim, nquantsplot, qplot1, qplot2, comp_stat_plot,
   ppc11, ppc22, ppc33, ppc44, ppc55, ppc66, ppn11, qplot11, qplot22 )


#######################comment for running all the stations at once ###########
} # end loop on stations
################################################################################