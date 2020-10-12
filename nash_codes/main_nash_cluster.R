# exploratory analysis of NASH data


rm(list=ls()) 
library(nleqslv)
library(extraDistr)
# library(ggplot2)
# library(tidyr)
# library(dplyr)
library(gridExtra)
library(rstan)
library(bayesplot)
library(loo)
# library(VGAM)
# library(evd)
# library(tseries)
library(data.table)


if (Sys.getenv('SLURM_CPUS_PER_TASK') == ""){
  use_cluster = FALSE
} else {
  use_cluster = TRUE
}

args=(commandArgs(TRUE))
if(length(args)==0){
  # print("No arguments supplied. Default is cross validation")
  ##supply default values: no cross validation
  test = TRUE
  dataset = 'G' # default dataset
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(sprintf("using cluster = %s", use_cluster))
print(sprintf('is this a test = %s', test))
print(sprintf('dataset = %s', dataset))


if (use_cluster){
  # setwd( file.path('~','hbev','codes'))
  # rstan_options(auto_write = TRUE)
  mc.cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
  numj = Sys.getenv("SLURM_ARRAY_TASK_ID")
  numj = as.integer(numj)
  datadir = file.path('..', '..', 'hmev2020si', 'data', 'Data_GHCN_Daily')
} else {
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  numj = 14  # index of station record to read
  # numj = 1111  # index of station record to read (WISCONSIN - SPOONER)
  # datadir = file.path('..', '..', 'hbev', 'data', 'Data_GHCN_Daily')
  datadir = file.path('..', 'data', 'Data_GHCN_Daily')
}




source("nashfun.R")    # main functions for data analysis
# source("hbev_module.R")    # main functions for data analysis
# source("hbev_functions.R") # other functions

dfs = read.csv(file.path(datadir, 'list_ushcn_stats.csv'))
datadir_stats = file.path(datadir, 'extracted_csv')
filenames = list.files(datadir_stats)

check_station_length = FALSE # do it only once

Nt = 92 # (30 + 31 + 31 expected observations)
#######################  parameters to set for MCMC  ###########################
cores = 4 # let us try to set this for both PC and cluster
min_num_years = 20 # if less skip station (years with enough data & overlap with NASH data)
maxmiss_per_season = 4 # remove all years with more than 4 events /season missing
thresh = 1 # for ordinary events
iter = 2000
chains = 4
ndraws = iter*chains/2
# thresh = 0
# thresh_hbev = 1
# adapt_delta = 0.8
# Mgen = 100
target_tr_20 = 20 # Tr for which alpha is computed 
target_tr_50 = 50 # Tr for which alpha is computed 
NEA = c(2, 4) # number of excesses / year to compute thresholds

# FiM = 1:M/(1+M)
# TrM = 1/(1-FiM)

################################################################################


# setwd( file.path('~','Projects','hbev','nash_codes'))
# outplot = file.path('..', 'nash_output', 'figures')
outdata1 = file.path('..', 'nash_output') # create it if does not exosts
dir.create(outdata1, showWarnings = FALSE)
outdata = file.path('..', 'nash_output', sprintf('stats_out_th_%s', thresh)) # create it if does not exosts
dir.create(outdata, showWarnings = FALSE)
outfilename = file.path(outdata, sprintf('res_stat_%s.csv', numj))

# if (file.exists(outfilename)){
#   print('Already done!')
# } else {



# load list of stations (and compute their length if not done)
# change it with "if exists file = ..."
# if (check_station_length){
# 
#   # load station data
#   dfs = read.csv(file.path('..', 'data', 'Data_GHCN_Daily', 'list_ushcn_stats.csv'))
#   nstats = dim(dfs)[1]
#   nyears_complete = rep(0, nstats)
#   for (is in 1:nstats){
#     print(is)
#     df = read_summer_data(file.path(datadir, sprintf('%s.csv', dfs$ID[is])), nmaxmiss = 4, Nt_JJA = Nt)
#     years = unique(df$YEAR)
#     nyears_complete[is] = length(years)
#   }
#   dfs$COMPLETE_YEARS = nyears_complete
#   write.csv(dfs, file.path(outdata, 'list_stats_complete_years.csv'))
# } else {
#   dfs = read.csv(file.path(outdata, 'list_stats_complete_years.csv'))
# }


  
  nstats = 1
  myID = as.character(dfs$ID[numj])
  myLAT = as.character(dfs$LATITUDE[numj])
  myLON = as.character(dfs$LONGITUDE[numj])

  
if (use_cluster){
  models = c(
             'nash_totals_wei_xydep',
             'nash_totals_wei_yydep' ,
             'nash_totals_wei_nodep' ,
             'nash_totals_wei_xxdep' ,
             # 'nash_xij_weiNS_xydep' ,
             'nash_xij_weiNSO_xydep' ,
             'nash_xij_wei_xydep' ,
             'nash_xij_wei_xxdep' ,
             'nash_xij_wei_yydep' ,
             'nash_xij_wei_nodep',
             # 'nash_xij_wei_xydepw' ,
             # 'nash_xij_wei_yydepw' ,
             # 'nash_xij_weiNS_xydepw' ,
             'nash_nj_bin_xydep' ,
             'nash_nj_bin_yydep' ,
             'nash_nj_bin_xxdep' ,
             'nash_nj_bin_nodep',
             'nash_njmc_markov_nodep',
             'nash_njmc_markov_yydep',
             'nash_njmc_markov_xxdep',
             'nash_njmc_markov_xydep',
             'nash_ne3j_potbin3_xydep' ,
             'nash_ne3j_potbin3_yydep' ,
             'nash_ne3j_potbin3_xxdep' ,
             'nash_ne3j_potbin3_nodep',
             'nash_ne5j_potbin5_xydep' ,
             'nash_ne5j_potbin5_yydep' ,
             'nash_ne5j_potbin5_xxdep' ,
             'nash_ne5j_potbin5_nodep'
             )
} else {
  # to test only
  models = c(
             'nash_totals_wei_xydep',
             'nash_totals_wei_yydep' ,
             'nash_totals_wei_nodep' ,
             'nash_totals_wei_xxdep' ,
             'nash_xij_weiNS_xydep' ,
             'nash_xij_weiNSO_xydep',
             'nash_xij_wei_xydep' ,
             'nash_xij_wei_xxdep' ,
             'nash_xij_wei_yydep' ,
             'nash_xij_wei_nodep'
             # 'nash_xij_wei_xydepw' ,
             # 'nash_xij_wei_yydepw' ,
             # 'nash_xij_weiNS_xydepw' ,
             # 'nash_nj_bin_xydep' ,
             # 'nash_nj_bin_yydep' ,
             # 'nash_nj_bin_xxdep' ,
             # 'nash_nj_bin_nodep',
             # 'nash_njmc_markov_nodep',
             # 'nash_njmc_markov_yydep',
             # 'nash_njmc_markov_xxdep',
             # 'nash_njmc_markov_xydep',
             # 'nash_ne3j_potbin3_xydep' ,
             # 'nash_ne3j_potbin3_yydep' ,
             # 'nash_ne3j_potbin3_xxdep' ,
             # 'nash_ne3j_potbin3_nodep',
             # 'nash_ne5j_potbin5_xydep' ,
             # 'nash_ne5j_potbin5_yydep' ,
             # 'nash_ne5j_potbin5_xxdep' 
             # 'nash_ne5j_potbin5_nodep'
             )

}

nmodels = length(models)



params = c("lppd", "lpml", "qDqTR1", "qDqTR2","qnormTR1", "qnormTR2",
           'c0', 'cy', 'cx', 'w0', 'wy', 'wx', 'n0', 'ny', 'nx',
           'n30', 'n3y', 'n3x', 'n50', 'n5y', 'n5x',
           'pw0', 'pwy', 'pwx', 'cr0', 'cry', 'crx'
           )

nparams = length(params)


mycube <- array(NaN,
          dim = c(nstats, nstats, nstats, nmodels, nparams),
          dimnames = list(lat = myLAT, lon = myLON, id = myID, model = models, param = params))


dfr = as.data.table(mycube, value.name = 'score', na.rm = FALSE)


  df = read_summer_data(file.path(datadir_stats, sprintf('%s.csv', dfs$ID[numj])), nmaxmiss = maxmiss_per_season)
  prcp_years = unique(df$YEAR)
  
  # write on csv the length of the shared dataset
    lengthfile = file.path(outdata1, 'length_file.csv')
  xxxx = data.frame(list('numj'=numj, 'ID'=myID, 'NASHYEARS'=length(prcp_years)))
      if (!(file.exists(lengthfile))) {
      file.create(lengthfile)
      write.table(xxxx, file=lengthfile, append=T, row.names = F, col.names = T)
      } else {
      write.table(xxxx, file=lengthfile, append=T, row.names = F, col.names = F)
      }
  
 # write down the file number if not enough years, if not continue
  if (length(prcp_years) < min_num_years){
    outwarningfile = file.path(outdata1, 'skipped_short_stats.csv')
    if (!(file.exists(outwarningfile))) {
      file.create(outwarningfile)
    }
      write(numj, file=outwarningfile, append=T)
  } else {
  # skip all the following analysis
  
  # load nash data and keep only years also in the PRCP time series
    
dfnash <- read_nash_values(prcp_years = prcp_years)
xnn = dfnash$XNN
ynn = dfnash$YNN

  # ynash = read.table(file.path('nash_pos_data', 'ridge_lat_data.txt'), header=FALSE, col.names = c('YEAR', 'LAT'))
  # xnash = read.table(file.path('nash_pos_data', 'ridge_lon_data.txt'), header=FALSE, col.names = c('YEAR', 'LON'))
  # xnash$LON = -xnash$LON
  # compute mean and sd of NASH before subsetting the years
  # xnash = subset(xnash, xnash$YEAR %in% prcp_years)
  # ynash = subset(ynash, ynash$YEAR %in% prcp_years)
  # xnn = (xnash$LON - mean(xnash$LON))/sd(xnash$LON)
  # ynn = (ynash$LAT - mean(ynash$LAT))/sd(ynash$LAT)
  
  # compute annual maxima statistics
  res = table_max(df, Nt=Nt, thresh=thresh)
  ndf = data.frame( list(xnn = xnn, ynn = ynn, max = res$max, 
                         N = res$N, totals = res$totals, 
                         means =res$totals/res$N, 
                         sdwets = res$sdwets, skews = res$skews))
    
  # compute frequency of excesses over some thresholds

  NP = list()
  npots = length(NEA)
  xiall = sort(df$PRCP, decreasing=TRUE)
  thresh_pot = rep(0, npots)
  for (i in 1:npots){
  thresh_pot[i] = xiall[res$nyears*NEA[i]+1]
  
  # compute number of excesses / year 
  Np = array(0, res$nyears)
  for (k in 1:res$nyears){
    all_events = df$PRCP[df$YEAR == res$years[k]]
    wets = all_events[all_events > thresh_pot[i]]
    Np[k] = length(wets)
  } 
  NP[[i]] = Np
  }
  
  for (im in 1:nmodels){
    
    
                      # if (models[im] == 'nash_xij_weiNS_xydep' ){
      if (models[im] %in% c('nash_xij_weiNS_xydep', 'nash_xij_weiNSO_xydep') ){
                        Mgen = 100
                      } else {
                        Mgen= 2
                      }
                      XNgen =  rep(-0.0, Mgen) # center
                      YNgen =  rep(0.0,  Mgen) # center
                      XNgenS = rep(-1.0, Mgen) # west
                      YNgenS = rep(-1.0, Mgen) # south
                      XNgenN = rep(-1.0, Mgen) # west
                      YNgenN = rep(1.0,  Mgen) # north
                      FiMgen = 1:Mgen/(1+Mgen)
                      TrMgen = 1/(1-FiMgen)
                  
    totnwets = length(res$data[res$data > 0]) #these are already excesses
    nash_data = list(
              M = length(ndf$means),
              Mgen = Mgen, # generate same length 
              Nt = Nt,
              y = res$data,
              totals = res$totals,
              N = ndf$N,
              XN = ndf$xnn, YN = ndf$ynn,
              # XNgen = xnash60n, YNgen = ynash60n, # generate from same values
              # XNgen = , YNgen = ynash60n, 
              XNgen = XNgen, YNgen = YNgen, # generate from same values 
              XNgenS = XNgenS, YNgenS = YNgenS, # generate from same values 
              XNgenN = XNgenN, YNgenN = YNgenN, # generate from same values 
              # XNgen = xnash60n, YNgen = ynash60n, # generate from same values 
              totnwets = totnwets,
              NP3 = NP[[1]],
              NP5 = NP[[2]]
              )
    
    # check for divergences
    n_err_max = 1
    while (n_err_max > 0){
      print(sprintf("station %s - fitting model %s", numj, models[im]))

      
      fit <- stan(file = sprintf('%s.stan', models[im]), data = nash_data, 
                    iter =iter, chains = chains, refresh = 0, cores = cores)
    
      # save diagnostic quantities:
      n_diverg = get_num_divergent(fit)
      n_max_tree = get_num_max_treedepth(fit)
      n_low_bfmi = length(get_low_bfmi_chains(fit))
      # check no divergences occurred and repeat if necessary
      n_err_max = max(n_diverg, n_max_tree, n_low_bfmi)
    } # end of while - section of code to repeat if divergence occur
    
    sim <- rstan::extract(fit)
    infc = info_crit(sim$log_lik, psis = FALSE, is_log_like = TRUE)
    
    # plot(fit)
              
    dfr[model == models[im] & param == 'lppd']$score = infc$lppd
    dfr[model == models[im] & param == 'lpml']$score = infc$lpml
    
    parnames = names(sim)
    for (ip in 1:nparams){
      if (params[ip] %in% parnames){
        dfr[model == models[im] & param == params[ip]]$score =  mean(sim[[params[ip]]])
      }
    }
    
      # if (models[im] %in% quantile_models){
      if (models[im] %in% c('nash_xij_weiNS_xydep', 'nash_xij_weiNSO_xydep') ){
      
      XIso = apply(sim$maxrepS, 1, sort) + thresh
      XIno = apply(sim$maxrepN, 1, sort) + thresh
        
      # XIso_mean = apply(sim$maxrepS, 2, mean)
      # XIno_mean = apply(sim$maxrepN, 2, mean)
      # XIso_std = apply(sim$maxrepS, 2, sd)
      # XIno_std = apply(sim$maxrepN, 2, sd)
      XIso_mean = apply(XIso, 1, mean)
      XIno_mean = apply(XIno, 1, mean)
      XIso_std = apply(XIso, 1, sd)
      XIno_std = apply(XIno, 1, sd)
      
      # mean(sim$xijrepN)
      # mean(sim$xijrepS)

      # XIno_mean = apply(XIno, 1, mean)
      # XIso_mean = apply(XIso, 1, mean)
      # XIno_std = apply(XIno, 1, sd)
      # XIso_std = apply(XIso, 1, sd)
      
      # plot(TrMgen, XIso_mean)
      # lines(TrMgen, XIno_mean, col='red')
      
      postr20 = which.min(abs(TrMgen - target_tr_20))
      sqTR20 = XIno_std[postr20]
      dqTR20 = XIno_mean[postr20] - XIso_mean[postr20] 
      
      postr50 = which.min(abs(TrMgen - target_tr_50))
      sqTR50 = XIno_std[postr50]
      dqTR50 = XIno_mean[postr50] - XIso_mean[postr50] 
      
      qDq20 = dqTR20/sqTR20 
      qDq50 = dqTR50/sqTR50 
      
      qnorm20 = dqTR20/XIso_mean[postr20] 
      qnorm50 = dqTR50/XIso_mean[postr50] 
      
      dfr[model == models[im] & param == 'qDqTR1']$score = qDq20
      dfr[model == models[im] & param == 'qDqTR2']$score = qDq50
      dfr[model == models[im] & param == 'qnormTR1']$score = qnorm20
      dfr[model == models[im] & param == 'qnormTR2']$score = qnorm50

      } # end if on model depxy
    } # end loop on models
  write.csv(dfr, outfilename) # save results
  } # if not enough years of data
# } # if output already exists
  
# pairs(fit, pars=c('w','c0', 'cx', 'cy'))
# pairs(fit, pars=c('n50', 'n5x'))
# pairs(fit, pars=c('w','C', 'p0'))

  # plot(fit, pars=('w'))

  
  # sum(is.na(sim$w))
