

read_summer_data <- function(filepath, nmaxmiss = 10, Nt_JJA = 92){
  df = read.csv(filepath)
  df$YEAR = floor(df$DATE/10000)
  df$MONTH = floor((df$DATE-df$YEAR*10000)/100)
  df = subset(df, df$YEAR >= 1948 & df$YEAR <= 2019)
  df = subset(df, df$MONTH %in% c(6, 7, 8)) # summer precipitation only
  df$PRCP = df$PRCP/10 # to millimeter accumulations
  
  # remove years with less than nmaxmiss missing data:
  df = subset(df, df$PRCP > -1 & !is.na(df$PRCP))
  # Nt_JJA = 92 # (30 + 31 + 31 expected observations)
  min_number_obs = Nt_JJA -nmaxmiss
  years = unique(df$YEAR)
  nyears = length(years)
  for (iy in 1:nyears){
    df_iy <- subset(df, df$YEAR == years[iy])
    # print(dim(df_iy)[1])
    if (dim(df_iy)[1] < min_number_obs){
     df = subset(df, df$YEAR != years[iy]) 
    }
  }
return(df)
}

skewness <- function(x){
  m <- mean(x)
  s = sd(x)
  skew = mean((x-m)^3)/(s^3)
}


read_nash_values <- function(prcp_years=1948:2019){
  # """-------------------------------------------------------------------------
  # Read Nash values and return those corresponding to the years in prcp_years
  # (by default all years)
  # Return data frame with years, ridge LAT and LON and normalized position
  # (XNN=normalized long, YNN=normalized lat)
  # using ALL YEARS 1948-2007 for the normalization
  # ------------------------------------------------------------------------"-""
  ynash = read.table(file.path( 'nash_pos_data', 'ridge_lat_data.txt'), 
                     header=FALSE, col.names = c('YEAR', 'LAT'))
  xnash = read.table(file.path( 'nash_pos_data', 'ridge_lon_data.txt'), 
                     header=FALSE, col.names = c('YEAR', 'LON'))
  xnash$LON = -xnash$LON # use West negative (x increases toward the East)
  # plot(xnash$LON, ynash$LAT)
  # compute mean and sd of NASH before subsetting the years
  mean_xnash = mean(xnash$LON)
  mean_ynash = mean(ynash$LAT)
  stdv_xnash = sd(xnash$LON)
  stdv_ynash = sd(ynash$LAT)
  xnash = subset(xnash, xnash$YEAR %in% prcp_years)
  ynash = subset(ynash, ynash$YEAR %in% prcp_years)
  # compute the normalized ridge position - use all years though
  # xnn = (xnash$LON - mean(xnash$LON))/sd(xnash$LON)
  # ynn = (ynash$LAT - mean(ynash$LAT))/sd(ynash$LAT)
  xnn = (xnash$LON - mean_xnash)/stdv_xnash
  ynn = (ynash$LAT - mean_ynash)/stdv_ynash
  dfn = data.frame(YEAR=prcp_years, LON=xnash$LON, LAT=ynash$LAT, XNN=xnn, YNN=ynn)
return(dfn)
}


ridgepos <- function(xnn=0, ynn=0){
  # compute quadrants based on ridge position
  # use normalized ridge positions in input
  nelem = length(xnn)
  quadrant = rep(0, nelem)
  for (i in 1:nelem){
      x = xnn[i]
      y = ynn[i]
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



table_max <- function(df, Nt = 92, thresh=0){
  "---------------------------------------------------
  df = dataframe with field $PRCP with the variable
  of interest, and with a field $YEAR integer
  in format YYYY used to divide the dataset in blocks.
  ----------------------------------------------------
  Nt (default Nt = 92) Number of events / season
  in the output matrix

  ----------------------------------------------------
  from a dataframe df compute matrix of yearly data
  of dimension nyears*Nt,
  the vector N with the number of events > thresh / year
  the array max of annual maxima,
  Xi = sorted maxima,
  Fi = their non exceedance frequencies,
  Tr = their empirical return times
  ----------------------------------------------------
  NOTE: The 366th day year is thrown away if Nt = 365
  ----------------------------------------------------"
  years  = unique(df$YEAR)
  nyears = length(years)
  maxima = array(0, nyears)
  minima = array(0, nyears)
  totals = array(0, nyears)
  sdwets = array(0, nyears)
  skews = array(0, nyears)
  N = array(0, nyears)
  datamat = matrix(0, nyears, Nt)
  count = 0
  all_data = rep(0, Nt)
  for (i in 1:nyears){
    all_events = df$PRCP[df$YEAR == years[i]]
    if (length(all_events)>Nt){
      print('table_max ERROR: some years of record have more than Nt obs.')
    }
    wets = all_events[all_events > thresh]-thresh # excesses
    N[i] = length(wets)
    if (N[i] > 0){
      maxima[i] = max(wets) + thresh
      minima[i] = min(wets) + thresh
      totals[i] = sum(wets) + N[i]*thresh # totals only counting events exceeding
      sdwets[i] = sd(wets)
      skews[i] = skewness(wets)
      nevents = length(all_events)
      datamat[i,1:nevents] = pmax(all_events-thresh, 0) # write excesses or zeros
    }
    else{
      count = count + 1
      sprintf('table_max WARNING: some %s years have no non-zero observations!', count)
    }
  }
  Fi = (1:nyears)/(1+nyears)
  Xi = sort(maxima)
  Tr = 1/(1-Fi)
  res = list(data = datamat, max = maxima, Fi = Fi, Xi = Xi, Tr = Tr, 
             N = N, years = years, nyears = nyears, min = minima,
             totals = totals, sdwets = sdwets, skews = skews)
  return(res)
}





info_crit <- function(pdfmat, psis = TRUE, is_log_like = FALSE){
  "-----------------------------------------------------------------------------
  given a matrix of p(yi | \theta^s) compute various information criteria
  pdfmat -> matrix with shape S*M, where
  S -> number of draw from the posterior = number of rows
  M -> sample size used for evaluation = number of columns
  if psis = TRUE: compute also pareto smoothed loo
  
  NB: In current version, lppd-lpml-elpd are computed with a -1/M factor,
  where M is the sample size, to make them positive (smallest = better)
  and more narrowly distributed
  
  NB: if turn is_log_like = TRUE, provide already log likelihood instead of pdf
  -----------------------------------------------------------------------------"
  S = dim(pdfmat)[1]
  M = dim(pdfmat)[2]
  if (is_log_like){
    loglik = pdfmat
    pdfmat = exp(loglik)
  } else {
  loglik = log(pdfmat)
  }
  
  cpoi = rep(0, M)
  mean_pdfi = rep(0, M)
  mean_logl = rep(0, M)
  var_logl = rep(0, M)
  # for (i in 1:M){
  #     cpoi[i] = 1/mean(1/pdfmat[,i])
  #     mean_pdfi[i] = mean(pdfmat[,i])
  #     mean_logl[i] = mean(loglik[,i])
  #     var_logl[i] = var(loglik[,i])
  # }
    for (i in 1:M){
      cpoi[i] = 1/mean(1/pdfmat[,i])
      mean_pdfi[i] = mean(pdfmat[,i])
      mean_logl[i] = mean(loglik[,i])
      var_logl[i] = var(loglik[,i])
  }
  lpml_sum = sum(log(cpoi))      
  lppd_sum = sum(log(mean_pdfi)) 
  lpml = -mean(log(cpoi))       # should be +sums
  lppd = -mean(log(mean_pdfi))  # should be +sums
  
  p_waic1 = 2*(lppd_sum - sum(mean_logl)) # effective number of parameters, waic 1
  p_waic2 = sum(var_logl) # effective number of parameters, waic 2
               
  # elpd_waic1 = lppd_sum - p_waic1
  # elpd_waic2 = lppd_sum - p_waic2
  elpd_waic1 = -1/M*(lppd_sum - p_waic1) # should be +sums
  elpd_waic2 = -1/M*(lppd_sum - p_waic2) # should be +sums
  
  if (psis == TRUE){
    psis_loo = loo(loglik)
    elpd_loo = psis_loo$estimates[1] # elpd, psis-loo
    p_loo = psis_loo$estimates[2] # eff number of parameters, psis-loo
  } else {
    elpd_loo = 0
    p_loo = 0
  }
  
 infocrits = list(lpml = lpml, lppd = lppd, elpd_waic1 = elpd_waic1, 
                  elpd_waic2 = elpd_waic2,
                  p_waic1 = p_waic1, p_waic2 = p_waic2, 
                  elpd_loo = elpd_loo, p_loo = p_loo)
 return(infocrits)
 }



# compute transition frequencies for each year of the dataset
# data = res$data
# nyears=res$nyears

gener_2smc <- function(wetf=0.5, corr=0, nyears=1, nobs=10){
  # generate nyears of data, each of length nobs
  # return matrix nyears*nobs
  # according to 2-states markov chain 
  # with probability of ones = wetf, serial correlation = 'corr'
  p11 = corr + wetf*(1-corr)
  p01 = p11 - corr
  data = matrix(0, nyears, nobs)
  for (i in 1:nyears){
    sample = rep(0, nobs)
    # generate markov chain data
    # random starting point
    sample[1] = rbern(1, wetf)
    for (j in 2:nobs){
      if (sample[j-1] == 0){ # zero
      sample[j] = rbern(1, p01)
      } else { # one
      sample[j] = rbern(1, p11) 
      }
    }
    data[i, ] = sample
  }
  return(data)
}

trans_prob <- function(data, thresh=0){
# compute transition prob for each year of data
# above / below threshold
# data = nyears*nobs matrix
nyears = dim(data)[1]
nobs = dim(data)[2]
n01 = rep(0, nyears)
n11 = rep(0, nyears)
n10 = rep(0, nyears)
n00 = rep(0, nyears)
# loop on years:
for (i in 1:nyears){
  datai = data[i, ]
  for (j in 2:Nt){
    if (datai[j] > thresh){
      # one
      if (datai[j-1] > thresh){ # one one
        n11[i] = n11[i] + 1
      } else { # zero one
        n01[i] = n01[i] + 1
      }
    } else {
      # zero
      if (datai[j-1] > thresh){ # one zero
        n10[i] = n10[i] + 1
      } else { # zero zero
        n00[i] = n00[i] + 1
      }
    }
  }
}
p01 = n01/(n01+n00)
p00 = n00/(n01+n00)
p10 = n10/(n10+n11)
p11 = n11/(n10+n11)
corr = p11-p01
wetf = p01/(1-corr)
transp = list('p01'=p01, 'p10'=p10, 'p00'=p00, 'p11'=p11, 'wetf'=wetf, 'corr'=corr)
return(transp)
}


plot_quadrant_stats <- function(res, xnn, ynn, 
                                station_id, station_name, station_state){
  
  resquad = ridgepos(xnn=xnn, ynn=ynn)
  nse = length(resquad[resquad == 'SE'])
  nsw = length(resquad[resquad == 'SW'])
  nne = length(resquad[resquad == 'NE'])
  nnw = length(resquad[resquad == 'NW'])
  
  ndf = data.frame( list(xnn = xnn, ynn = ynn, max = res$max, 
                       N = res$N, totals = res$totals, 
                       means =res$totals/res$N, 
                       sdwets = res$sdwets, skews = res$skews))


  # p0 <- ggplot(ndf)+
  #   geom_point(aes(x=xnn, y=ynn, color = max, size = 1))
  
  p1 <- ggplot(ndf, aes(x=resquad, y=totals, fill = resquad, alpha = 0.6)) +
    # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
    geom_boxplot() + 
    theme_bw() + 
    ggtitle('Seasonal totals') + 
    xlab('NASH Quadrant')+
    ylab(expression(S[j] *  " [mm]"))+
    guides(fill=FALSE, alpha=FALSE)+
    geom_jitter(width=0.25, alpha=0.5) 
  
  p2 <- ggplot(ndf, aes(x=resquad, y=N, fill = resquad, alpha = 0.6)) +
    # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
    geom_boxplot() + 
    theme_bw() + 
    ggtitle('Seasonal number of events') + 
        xlab('NASH Quadrant')+
    # ylab(paste(expression(n[j]), ' [mm]'))+
    ylab(expression(n[j]))+
    guides(fill=FALSE, alpha=FALSE)+
    geom_jitter(width=0.25, alpha=0.5) 
  
  p3 <- ggplot(ndf, aes(x=resquad, y=means, fill = resquad, alpha = 0.6)) +
    # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
    geom_boxplot() + 
    theme_bw() + 
    ggtitle('Average event magnitude') + 
            xlab('NASH Quandrant')+
    ylab(expression("<" * x[ij] * ">" * " [mm/day]"))+
    guides(fill=FALSE, alpha=FALSE)+
    geom_jitter(width=0.25, alpha=0.5) 
  
  p4 <- ggplot(ndf, aes(x=resquad, y=sdwets, fill = resquad, alpha = 0.6)) +
    # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
    geom_boxplot() + 
    theme_bw() + 
    guides(fill=FALSE, alpha=FALSE)+
    ggtitle('Events standard deviation') + 
    xlab('NASH Quadrant')+
    ylab(expression(sigma[x[ij]] * " [mm/day]"))+
    geom_jitter(width=0.25, alpha=0.5) 
  
  p5 <- ggplot(ndf, aes(x=resquad, y=sdwets, fill = resquad, alpha = 0.6)) +
    # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
    geom_boxplot() + 
    theme_bw() + 
    guides(fill=FALSE, alpha=FALSE)+
    ggtitle('Events skewness') + 
                xlab('NASH Quadrant')+
    ylab(expression("<" * x[ij]^3 * ">" / sigma[x[ij]]^3))+
    geom_jitter(width=0.25, alpha=0.5) 
  
  p6 <- ggplot(ndf, aes(x=resquad, y=max, fill = resquad, alpha = 0.6)) +
    # geom_boxplot(aes(x=resquad, y=N, color = resquad)) + 
    geom_boxplot() + 
    ggtitle('Seasonal maxima') + 
                    xlab('NASH Quadrant')+
    # ylab(expression(" max x"[ij]  *  "[mm/day]"))+
    ylab(expression(x[ij]^"(m)"  *  " [mm/day]"))+
    geom_jitter(width=0.25, alpha=0.5)  +
    theme_bw() + 
    guides(fill=FALSE, alpha=FALSE)+
    xlab('NASH Quadrant')
    # ylab(expression( max[j] * " {" * x[ij] * "}"))
    # labs( x= "Ridge position")
  
  fig <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2,
                  top = sprintf("Station %s %s (%s)", 
                  station_id, station_name, station_state))
  return(fig)
}