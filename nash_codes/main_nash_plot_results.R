
rm(list=ls()) 
library(dplyr)
library(mapproj)
library(grid)
library(scales)
library(ggplot2)
library(data.table)
library(gridExtra)
library(dplyr)
library(viridis)
# Get the world polygon and extract the USA
library(maps)
library(stringr)
library(mapdata)
library(latex2exp)
  library("sf")
library("rnaturalearth")
library("rnaturalearthdata")


source("nashfun.R")    # main functions for data analysis

# ALP = -0.1
# XXX <- seq(from = -5, to = 5, by = 0.2)
# YYY1 <- 1/(1+exp(+XXX))
# YYY <- 1/(1+exp(-ALP*XXX))
# plot(XXX,YYY1, lty=1)
# lines(XXX,YYY, lty=1, col = 'red')

# mymodel = 'nash_nj_bin'
# var = strsplit(mymodel, split = '_')[[1]][2]


  states <- map_data("state")
  usa <- map_data("usa")
  

theme_set(theme_bw())

dfnash <- read_nash_values()

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# setwd( file.path('~','Projects','hbev','nash_codes'))


outplot = file.path('..', 'nash_output_figures')
dir.create(outplot, showWarnings = FALSE)
outdata = file.path('..', 'nash_output', 'outdata')
# outdata_cluster = file.path('..', 'nash_output', 'stats_out')
outdata_cluster = file.path('..', 'nash_output', 'stats_out_th_1')



datalist <- list()
files = list.files(outdata_cluster)
nfiles = length(files)
# read in sequence and concat data
for (i in 1:nfiles) {
  dat <-  read.table( file.path(outdata_cluster, files[i]), sep = ',', 
                    header = TRUE, row.names = 1)
  datalist[[i]] <- dat # add it to your list
}
dfr <- data.table::rbindlist(datalist)


# check for missintg data:
NUMJ = rep(NaN, 1218)
# MISSING = 1:1218
MISSING = rep(NaN, 1218)
for (i in 1:1218){
  myfileii <- file.path(outdata_cluster, sprintf('res_stat_%s.csv',i))
  if (!(file.exists(myfileii))){
    print(i)
    MISSING[i] <- i  
  }
  # NUMJ[i] <- as.numeric(strsplit(strsplit(as.character(files[i]), split ='_')[[1]][3], '.csv')[[1]])
}
MISSING <- MISSING[na.omit(MISSING)]
SHORT <- read.table(file.path(outdata_cluster,'..', 'skipped_short_stats_20.csv'))$V1
# REALMISS = set(MISSING) - set(SHORT)
REALMISS = setdiff(MISSING, SHORT)
print(REALMISS)

names(dfr)[names(dfr) == 'id'] <- 'ID'
names(dfr)[names(dfr) == 'lat'] <- 'LATITUDE'
names(dfr)[names(dfr) == 'lon'] <- 'LONGITUDE'





# dfr$model[dfr$model == 'nash_nej_potbin3_depxy'] <- 'nash_ne3j_potbin3_depxy'
# dfr$model[dfr$model == 'nash_nej_potbin3_depyy'] <- 'nash_ne3j_potbin3_depyy'
# dfr$model[dfr$model == 'nash_nej_potbin3_depno'] <- 'nash_ne3j_potbin3_depno'
# dfr$model[dfr$model == 'nash_nej_potbin5_depxy'] <- 'nash_ne5j_potbin5_depxy'
# dfr$model[dfr$model == 'nash_nej_potbin5_depyy'] <- 'nash_ne5j_potbin5_depyy'
# dfr$model[dfr$model == 'nash_nej_potbin5_depno'] <- 'nash_ne5j_potbin5_depno'
# 
# dfr$model[dfr$model == 'nash_nj_markov_depxy'] <- 'nash_njmc_markov_depxy'
# dfr$model[dfr$model == 'nash_nj_markov_depyy'] <- 'nash_njmc_markov_depyy'
# dfr$model[dfr$model == 'nash_nj_markov_depno'] <- 'nash_njmc_markov_depno'


# dfr = read.csv(file.path(outdata, 'res_allstats_188.csv'))
# read position (LAT, LON) of all stations/ runs
# names(dfr)[names(dfr) == 'stat'] <- 'ID'
# # load stations positions and subset them
# posfile = file.path('..', 'data', 'Data_GHCN_Daily', 'list_ushcn_stats.csv')
# dfpos = read.table(posfile, sep=',', header = TRUE)
# dfpos <- dfpos[c("LATITUDE", "LONGITUDE", "ID")]
# dfr <- merge(dfpos, dfr, by = c("ID"))

# df3 <- reshape2::dcast(dfr, ID + LATITUDE + LONGITUDE +  numj 
#                 + gof + ssize + test ~ model, value.var = 'score')

varfun <- function(x) strsplit(as.character(x), split = '_')[[1]][2]
modfun <- function(x) tail(strsplit(as.character(x), split = '_')[[1]],1)

# dfrr = na.omit(dfr) # remove later should not change
dfr$var <- mapply(varfun, dfr$model)
dfr$mod <- mapply(modfun, dfr$model)

print(unique(dfr$param))

dfr$mod[dfr$mod == 'nodep'] <- 'NOD'
dfr$mod[dfr$mod == 'xydep'] <- 'LWLD'
dfr$mod[dfr$mod == 'yydep'] <- 'LATD'
dfr$mod[dfr$mod == 'xxdep'] <- 'LOND'

# dfr$mod <- factor(dfr$mod, levels = c("NOD", "LOND", "LATD", "LWLD"))
dfr$mod <- factor(dfr$mod, levels = c("LWLD", "LATD", "LOND", "NOD"))
# scale_fill_manual(values=c("darkred", "darkgreen", "lightblue", "grey40")) +
# scale_fill_manual(values=c("grey40", "darkgreen", "darkred", "lightblue")) +
# dfr$var <- as.factor(dfr$var)
# dfr$param <- as.factor(dfr$param)

sapply(dfr, typeof)

dfrr = dfr

# drops <- c("model")
# dfrr = dfrr[ , !(names(dfrr) %in% drops)]

dfrr = subset(dfrr, dfrr$param %in% c("lppd", "lpml") )

dfrr = subset(dfrr, !(dfrr$model %in% c("nash_xij_weiNSO_xydep") ))

dfrr$model <- NULL
df3 <- reshape2::dcast(dfrr, ID + LATITUDE + LONGITUDE + var
                + param ~ mod, value.var = 'score')
# df3 = subset(df3, df3$param %in% c("lppd", "lpml") )
# dfval = df3[,c("depxy", "depyy", "depno")]
# dfval = df3[,c("xydep", "yydep", "nodep", "xxdep")]
dfval = df3[,c("LWLD", "LATD", "NOD", "LOND")]
df3$best = colnames(dfval)[max.col(-dfval, ties.method="first")] #minimum = best




# unique(df3$var)
df3$varname = as.factor(df3$var)
levels(df3$varname) <- c(
                    "atop(n[j]~average~2, excesses/year)",
                    "atop(n[j]~average~4, excesses/year)",
                    "atop(n[j], Binomial)",
                    "atop(n[j], Markov~Chain)",
                    "atop(S[j]~seasonal, totals)",
                    "atop(x[ij]~daily, magnitudes)"
                    )


ggplot(df3) +
      # ggtitle( sprintf("%s sample size=%s, %s ", MYTEST[j], MYSSIZE[i], MYGOF[k])) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
    # geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill="grey", alpha=0.6) +
    # coord_quickmap() +  
    # geom_point(aes(x=LONGITUDE, y=LATITUDE,  color=best, shape=best)) +
    geom_point(aes(x=LONGITUDE, y=LATITUDE,  fill=best), alpha = 0.9, size = 0.8, color = 'black', shape = 21) +
  # ylab("Sample size")+
    labs(fill = "Best model")+
  guides(fill = guide_legend(override.aes = list(size = 4)))+
   
      # scale_color_manual(values=c("grey40", "darkgreen", "darkblue", "darkred")) +
# scale_color_manual(values=c("grey40", "darkgreen", "darkred", "lightblue")) +
scale_fill_manual(values=c("NOD"="grey40", "LATD"="darkgreen", "LOND"="darkred", "LWLD" = "lightblue")) +
  
    # annotate('text', x=-115, y=26, label=sprintf('bm = %s', round(freq[modelfrac], 2))) +
    # scale_fill_discrete(name = "best model")   + 
  #   su <- signup('enrico89', 'enrico.zorzetto@duke.edu', save = TRUE)
    # scale_size_continuous(range=c(1,12)) +
    # scale_color_viridis(trans="log") +
    theme_bw()+
                   theme(
               axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(size=15),
               strip.background =element_rect(fill="white"),
               legend.position = 'bottom') +
  facet_grid(varname~ param, labeller = label_parsed)+
    ggsave(file.path(outplot, sprintf('bestmodel_maps.png')),
          width = 6.5, height = 9, dpi = 300)

# report the best lplm model:
dfrefbest = subset(df3[c("ID", 'var', "best")], df3$param == "lpml")
dfr <- merge(dfrefbest, dfr, by = c("ID", 'var'))

extract_param <- function(dfk, myvar = 'totals', myparamz = c('cx', 'cy') ){
  # note must provide paramz in order x y
  dfk2 = subset(dfk, dfk$var == myvar & dfk$param %in% myparamz)
  dfk2 <- subset(dfk2, dfk2$best == dfk2$mod)
  # keep only the rows with the parameter of interest (nan, cx or cy)
  # dfk3 <- subset(dfk2, (dfk2$best == 'nodep' & dfk2$param == myparamz[1]) |
  #               (dfk2$best == 'xxdep' & dfk2$param == myparamz[1]) |
  #               (dfk2$best %in% c('yydep', 'xydep') & dfk2$param == myparamz[2])
  #                )
    dfk3 <- subset(dfk2, (dfk2$best == 'NOD' & dfk2$param == myparamz[1]) |
                (dfk2$best == 'LOND' & dfk2$param == myparamz[1]) |
                (dfk2$best %in% c('LATD', 'LWLD') & dfk2$param == myparamz[2])
                 )
  dfk3 <- transform(dfk3, effect = ifelse(score < 0, 'Negative', 'Positive'))
  dfk3$effect[dfk3$mod == 'NOD'] = 'No dependence'
  return(dfk3)
}

# seasonal totals
dfpartot <- extract_param(dfr, myvar = 'totals', myparamz = c('cx', 'cy'))


dfpartot_xx <- subset(dfpartot, dfpartot$best == 'LOND')
dfpartot_yy <- subset(dfpartot, dfpartot$best == 'LATD')
dfpartot_xy <- subset(dfpartot, dfpartot$best == 'LWLD')
dfpartot_yyxy <- subset(dfpartot, (dfpartot$best == 'LWLD' | dfpartot$best == 'LATD'))
dfpartot_yyxy$bestsign = paste(dfpartot_yyxy$best, dfpartot_yyxy$effect, sep=' ') 

par1a <- ggplot(dfpartot) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfpartot_xx, aes(x=LONGITUDE, y=LATITUDE,
                  fill = effect, shape=effect), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_manual(values = c('darkred', 'blue'))+
# scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            # 'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    # guides(# remove legend for triangles direction
    #  fill=guide_legend(title.position="top", title = 'Best model',
    #  override.aes = list(shape = 21, size = 4)))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
        geom_segment(
    x = -76, y = 34,
    xend = -65, yend = 34,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkgreen" # Also accepts "red", "blue' etc
  ) + 
  annotate(geom="text", x=-72, y= 36, label="RIDGE", color="darkgreen", size = 5)+
  # guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  annotate(geom="text", x=-73, y=30, label="a) LOND model", color="black", size=4) +
  annotate(geom="text", x=-72, y=28, label="Seasonal totals", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="Weibull~scale~delta[x]^(S)", color="black", size=4, parse=TRUE) +
    theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

par1b <- ggplot(dfpartot) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfpartot_yyxy, aes(x=LONGITUDE, y=LATITUDE,
                  fill = bestsign, shape=bestsign), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24, 25, 24))+
scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    # guides(# remove legend for triangles direction
    #  fill=guide_legend(title.position="top", title = 'Best model',
    #  override.aes = list(shape = 21, size = 4)))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
          geom_segment(
    x = -66, y = 34,
    xend = -66, yend = 40,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkgreen" # Also accepts "red", "blue' etc
  ) + 
  annotate(geom="text", x=-70, y= 36, label="RIDGE", color="darkgreen", size = 5)+
  annotate(geom="text", x=-72, y=30, label="b) LATD and LWLD", color="black", size=4) +
  annotate(geom="text", x=-72, y=28, label="Seasonal totals", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="Weibull~scale~delta[y]^(S)", color="black", size=4, parse=TRUE) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.position = c(0.20, 0.13), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())



    # ggsave(file.path(outplot, sprintf('res_maps_%s_%s_binary.png', myparam, mymodel)),
    #       width = 6.5, height = 4.5)

# myparam = 'ny'
# mymodel = 'nash_nj_bin_xydep'
# dfrrr = subset(dfr, dfr$model == mymodel & dfr$param == myparam)
# 
# # dfrr = subset(dfr, dfr$model == 'nash_nj_bin_depxy' & dfr$param == 'mcy')
# # dfrrr = subset(dfr, dfr$model == 'nash_xij_wei_depxy' & dfr$param == 'mcy')
# dfrrr <- transform(dfrrr, effect = ifelse(score < 0, 'Negative', 'Positive'))
# 
# # add info on the best model
# dfrefbest = subset(df3[c("ID", "best")], df3$param == "lpml" & df3$var == "nj")
# dfrrr <- merge(dfrefbest, dfrrr, by = c("ID"))


dfparnj <- extract_param(dfr, myvar = 'nj', myparamz = c('nx', 'ny'))

dfparnj_xx <- subset(dfparnj, dfparnj$best == 'LOND')
dfparnj_yy <- subset(dfparnj, dfparnj$best == 'LATD')
dfparnj_xy <- subset(dfparnj, dfparnj$best == 'LWLD')
dfparnj_yyxy <- subset(dfparnj, (dfparnj$best == 'LWLD' | dfparnj$best == 'LATD'))
dfparnj_yyxy$bestsign = paste(dfparnj_yyxy$best, dfparnj_yyxy$effect, sep=' ') 

par2a <- ggplot(dfparnj) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparnj_xx, aes(x=LONGITUDE, y=LATITUDE,
                  fill = effect, shape=effect), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_manual(values = c('darkred', 'blue'))+
# scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            # 'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    # guides(# remove legend for triangles direction
    #  fill=guide_legend(title.position="top", title = 'Best model',
    #  override.aes = list(shape = 21, size = 4)))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  # annotate(geom="text", x=-72, y=28, label="c) Seasonal totals", color="black", size=4) +
  # annotate(geom="text", x=-72, y=26, label="Weibull~scale~delta[y]^(S)", color="black", size=4, parse=TRUE) +
  annotate(geom="text", x=-73, y=30, label="c) LOND model", color="black", size=4) +
      annotate(geom="text", x=-72, y=28, label="Number of events", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="Binomial~rate~eta[x]", color="black", size=4, parse=TRUE)+
    theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

par2b <- ggplot(dfparnj) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparnj_yyxy, aes(x=LONGITUDE, y=LATITUDE,
                  fill = bestsign, shape=bestsign), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24, 25, 24))+
scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    # guides(# remove legend for triangles direction
    #  fill=guide_legend(title.position="top", title = 'Best model',
    #  override.aes = list(shape = 21, size = 4)))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  annotate(geom="text", x=-72, y=30, label="d) LATD and LWLD", color="black", size=4) +
    annotate(geom="text", x=-72, y=28, label="Number of events", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="Binomial~rate~eta[y]", color="black", size=4, parse=TRUE)+
  # annotate(geom="text", x=-72, y=28, label="a) Seasonal totals", color="black", size=4) +
  # annotate(geom="text", x=-72, y=26, label="Weibull~scale~delta[y]^(S)", color="black", size=4, parse=TRUE) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.position = c(0.20, 0.13), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

# par2 <- ggplot(dfparnj) +
#     geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
#     # geom_polygon(data = usa, aes(x=long, y = lat), fill="grey", alpha=0.9) +
#     # geom_point(aes(x=LONGITUDE, y=LATITUDE,
#     #                fill=best, shape = effect), size = 2.0, color = 'black') +
#       geom_point(aes(x=LONGITUDE, y=LATITUDE,
#                    fill=effect, shape = best), size = 2.0, color = 'black') +
#   # scale_shape_manual(values=c(25, 21, 24))+
#   scale_shape_manual(values=c(24, 22, 23, 21))+
# scale_fill_manual(values=c( "lightblue", "grey60", "red")) +
#     xlab('Longitude')+
#     ylab('Latitude')+
# # scale_fill_manual(values=c("grey40", "darkgreen", "lightblue", "darkred")) +
# # scale_fill_manual(values=c("darkred", "darkgreen", "lightblue", "grey40")) +
# # scale_fill_manual(values=c("NOD"="grey40", "LATD"="darkgreen", "LOND"="darkred", "LWLD" = "lightblue")) +
#   guides(shape=FALSE,  # remove legend for triangles direction
#      fill=guide_legend(title.position="top", title = 'Best model',
#      override.aes = list(shape = 21, size = 4)))+
#     theme_bw()+
#     theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
#           legend.background = element_rect(fill=alpha(0.4)))+
#   annotate(geom="text", x=-72, y=28, label="b) Number of events", color="black", size=4) +
#   annotate(geom="text", x=-72, y=26, label="Binomial~rate~eta[y]", color="black", size=4, parse=TRUE) 
#     # ggsave(file.path(outplot, sprintf('res_maps_%s_%s_binary.png', myparam, mymodel)),
#     #       width = 6.5, height = 4.5)

# daily rainfall intensities
dfr_noNSO = subset(dfr, dfr$model != 'nash_xij_weiNSO_xydep')
dfparcj <- extract_param(dfr_noNSO, myvar = 'xij', myparamz = c('cx', 'cy'))
dfparcj_xx <- subset(dfparcj, dfparcj$best == 'LOND')
dfparcj_yy <- subset(dfparcj, dfparcj$best == 'LATD')
dfparcj_xy <- subset(dfparcj, dfparcj$best == 'LWLD')
dfparcj_yyxy <- subset(dfparcj, (dfparcj$best == 'LWLD' | dfparcj$best == 'LATD'))
dfparcj_yyxy$bestsign = paste(dfparcj_yyxy$best, dfparcj_yyxy$effect, sep=' ') 

par3a <- ggplot(dfparcj) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparcj_xx, aes(x=LONGITUDE, y=LATITUDE,
                  fill = effect, shape=effect), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_manual(values = c('darkred', 'blue'))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+

  
  annotate(geom="text", x=-73, y=30, label="e) LOND model", color="black", size=4) +
      annotate(geom="text", x=-72, y=28, label="Daily magnitudes", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label=paste("Weibull~scale~delta[x]"), color="black", size=4, parse=TRUE)+ 
    theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

par3b <- ggplot(dfparcj) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparcj_yyxy, aes(x=LONGITUDE, y=LATITUDE,
                  fill = bestsign, shape=bestsign), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24, 25, 24))+
scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
  annotate(geom="text", x=-72, y=30, label="f) LATD and LWLD", color="black", size=4) +
    annotate(geom="text", x=-72, y=28, label="Daily magnitudes", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label=paste("Weibull~scale~delta[y]"), color="black", size=4, parse=TRUE)+ 
    theme(legend.position = c(0.20, 0.13), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())


  
  figpar <- grid.arrange(par1a, par1b, par2a, par2b, par3a, par3b, ncol=2)
    ggsave(file.path(outplot, sprintf('nash_param_ncs0.png')),
          width = 13, height = 12, dpi = 300, plot=figpar)
    
    
    
        
    ### MARKOV-CHAIN
# myparam = 'pwy'
mymodel = 'nash_njmc_markov_xydep'
# dfrrr = subset(dfr, dfr$model == mymodel & dfr$param == myparam)
# dfrefbest = subset(df3[c("ID", "best")], df3$param == "lpml" & df3$var == "njmc")
# dfrrr <- merge(dfrefbest, dfrrr, by = c("ID"))
# dfrrr <- transform(dfrrr, effect = ifelse(score < 0, 'Negative', 'Positive'))



dfparnjmc1 <- extract_param(dfr, myvar = 'njmc', myparamz = c('pwx', 'pwy'))

dfparnjmc1_xx <- subset(dfparnjmc1, dfparnjmc1$best == 'LOND')
dfparnjmc1_yy <- subset(dfparnjmc1, dfparnjmc1$best == 'LATD')
dfparnjmc1_xy <- subset(dfparnjmc1, dfparnjmc1$best == 'LWLD')
dfparnjmc1_yyxy <- subset(dfparnjmc1, (dfparnjmc1$best == 'LWLD' | dfparnjmc1$best == 'LATD'))
dfparnjmc1_yyxy$bestsign = paste(dfparnjmc1_yyxy$best, dfparnjmc1_yyxy$effect, sep=' ') 


dfparnjmc2 <- extract_param(dfr, myvar = 'njmc', myparamz = c('crx', 'cry'))

dfparnjmc2_xx <- subset(dfparnjmc2, dfparnjmc2$best == 'LOND')
dfparnjmc2_yy <- subset(dfparnjmc2, dfparnjmc2$best == 'LATD')
dfparnjmc2_xy <- subset(dfparnjmc2, dfparnjmc2$best == 'LWLD')
dfparnjmc2_yyxy <- subset(dfparnjmc2, (dfparnjmc2$best == 'LWLD' | dfparnjmc2$best == 'LATD'))
dfparnjmc2_yyxy$bestsign = paste(dfparnjmc2_yyxy$best, dfparnjmc2_yyxy$effect, sep=' ') 

njmc1a <- ggplot(dfparnjmc1) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparnjmc1_xx, aes(x=LONGITUDE, y=LATITUDE,
                  fill = effect, shape=effect), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_manual(values = c('darkred', 'blue'))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
          geom_segment(
    x = -76, y = 34,
    xend = -65, yend = 34,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkgreen" # Also accepts "red", "blue' etc
  ) + 
  annotate(geom="text", x=-72, y= 36, label="RIDGE", color="darkgreen", size = 5)+
  annotate(geom="text", x=-72, y=30, label="a) LOND model", color="black", size=4) +
      annotate(geom="text", x=-72, y=28, label="Wet~fraction", color="black", size=4, parse=TRUE) +
  annotate(geom="text", x=-74, y=26, label="beta[x]^(MC)", color="black", size=4, parse=TRUE) +
    theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

njmc1b <- ggplot(dfparnjmc1) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparnjmc1_yyxy, aes(x=LONGITUDE, y=LATITUDE,
                  fill = bestsign, shape=bestsign), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24, 25, 24))+
scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
            geom_segment(
    x = -66, y = 34,
    xend = -66, yend = 40,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkgreen" # Also accepts "red", "blue' etc
  ) + 
  annotate(geom="text", x=-70, y= 36, label="RIDGE", color="darkgreen", size = 5)+
  annotate(geom="text", x=-72, y=30, label="b) LATD and LWLD", color="black", size=4) +
        annotate(geom="text", x=-72, y=28, label="Wet~fraction", color="black", size=4, parse=TRUE) +
  annotate(geom="text", x=-74, y=26, label="beta[y]^(MC)", color="black", size=4, parse=TRUE) +
  #   annotate(geom="text", x=-72, y=26, label="(b)~Correlation", color="black", size=4, parse=TRUE) +
  # annotate(geom="text", x=-74, y=24, label="rho[y]^(MC)", color="black", size=4, parse=TRUE) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.position = c(0.20, 0.13), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

njmc2a <- ggplot(dfparnjmc2) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparnjmc2_xx, aes(x=LONGITUDE, y=LATITUDE,
                  fill = effect, shape=effect), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_manual(values = c('darkred', 'blue'))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  annotate(geom="text", x=-72, y=30, label="c) LOND model", color="black", size=4) +
      annotate(geom="text", x=-72, y=28, label="Correlation", color="black", size=4, parse=TRUE) +
  annotate(geom="text", x=-74, y=26, label="rho[x]^(MC)", color="black", size=4, parse=TRUE) +
    theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

njmc2b <- ggplot(dfparnjmc2) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparnjmc2_yyxy, aes(x=LONGITUDE, y=LATITUDE,
                  fill = bestsign, shape=bestsign), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24, 25, 24))+
scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  annotate(geom="text", x=-72, y=30, label="d) LATD and LWLD", color="black", size=4) +
    annotate(geom="text", x=-72, y=28, label="Correlation", color="black", size=4, parse=TRUE) +
  annotate(geom="text", x=-74, y=26, label="rho[y]^(MC)", color="black", size=4, parse=TRUE) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.position = c(0.20, 0.13), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

  figmcnj <- grid.arrange(njmc1a, njmc1b, njmc2a, njmc2b, ncol=2)
    ggsave(file.path(outplot, sprintf('res_maps_%s.png', mymodel)),
          width = 13, height = 9, dpi = 300, plot=figmcnj)




# pot1 <- ggplot(dfparnjmc1) +
#     geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
#     geom_point(aes(x=LONGITUDE, y=LATITUDE,
#                    fill = best, shape=effect), size = 2.0, color = 'black') +
#   # scale_shape_manual(values=c(6, 1, 2))+
#   scale_shape_manual(values=c(25, 21, 24))+
#       # guides(shape = FALSE,color=guide_legend(title.position="top", title = 'Best model'))+ 
#   guides(shape=FALSE,  # remove legend for triangles direction
#      fill=guide_legend(title.position="top", title = 'Best model',
#      override.aes = list(shape = 21, size = 4)))+
#     xlab('Longitude')+
#     ylab('Latitude')+
#     theme_bw()+
#     # scale_fill_manual(values = c("grey40", "darkgreen", "lightblue", "darkred"))+
# scale_fill_manual(values=c("NOD"="grey40", "LATD"="darkgreen", "LOND"="darkred", "LWLD" = "lightblue")) +
#     # scale_colour_manual(values = c("grey40", "darkgreen", "darkblue", "darkred"))+
#   annotate(geom="text", x=-72, y=26, label="(a)~Wet~fraction", color="black", size=4, parse=TRUE) +
#   annotate(geom="text", x=-74, y=24, label="beta[y]^(MC)", color="black", size=4, parse=TRUE) +
#     theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
#           legend.background = element_rect(fill=alpha(0.4)))

# myparam = 'cry'
# dfrrr = subset(dfr, dfr$model == mymodel & dfr$param == myparam)
# dfrefbest = subset(df3[c("ID", "best")], df3$param == "lpml" & df3$var == "njmc")
# dfrrr <- merge(dfrefbest, dfrrr, by = c("ID"))
# dfrrr <- transform(dfrrr, effect = ifelse(score < 0, 'Negative', 'Positive'))


# pot2 <- ggplot(dfparnjmc2) +
#     geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
#     geom_point(aes(x=LONGITUDE, y=LATITUDE,
#                    fill = best, shape=effect), size=2.0, color = 'black') +
#   # scale_shape_manual(values=c(6, 1, 2))+
#   scale_shape_manual(values=c(25, 21, 24))+
#       # guides(shape = FALSE,color=guide_legend(title.position="top", title = 'Best model'))+ 
#   guides(shape=FALSE,  # remove legend for triangles direction
#      fill=guide_legend(title.position="top", title = 'Best model',
#      override.aes = list(shape = 21, size = 4)))+
#     xlab('Longitude')+
#     ylab('Latitude')+
#     theme_bw()+
#   annotate(geom="text", x=-72, y=26, label="(b)~Correlation", color="black", size=4, parse=TRUE) +
#   annotate(geom="text", x=-74, y=24, label="rho[y]^(MC)", color="black", size=4, parse=TRUE) +
#     # scale_fill_manual(values = c("grey40", "darkgreen", "lightblue", "darkred"))+
# scale_fill_manual(values=c("NOD"="grey40", "LATD"="darkgreen", "LOND"="darkred", "LWLD" = "lightblue")) +
#     theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
#           legend.background = element_rect(fill=alpha(0.4)))








################## Frequency of Peaks over threshold ###########################
    
myparam3 = 'n3y'
mymodel3 = 'nash_ne3j_potbin3_xydep'
dfparne3j <- extract_param(dfr, myvar = 'ne3j', myparamz = c('n3x', 'n3y'))

dfparne3j_xx <- subset(dfparne3j, dfparne3j$best == 'LOND')
dfparne3j_yy <- subset(dfparne3j, dfparne3j$best == 'LATD')
dfparne3j_xy <- subset(dfparne3j, dfparne3j$best == 'LWLD')
dfparne3j_yyxy <- subset(dfparne3j, (dfparne3j$best == 'LWLD' | dfparne3j$best == 'LATD'))
dfparne3j_yyxy$bestsign = paste(dfparne3j_yyxy$best, dfparne3j_yyxy$effect, sep=' ') 

myparam5 = 'n5y'
mymodel5 = 'nash_ne5j_potbin5_xydep'
dfparne5j <- extract_param(dfr, myvar = 'ne5j', myparamz = c('n5x', 'n5y'))

dfparne5j_xx <- subset(dfparne5j, dfparne5j$best == 'LOND')
dfparne5j_yy <- subset(dfparne5j, dfparne5j$best == 'LATD')
dfparne5j_xy <- subset(dfparne5j, dfparne5j$best == 'LWLD')
dfparne5j_yyxy <- subset(dfparne5j, (dfparne5j$best == 'LWLD' | dfparne5j$best == 'LATD'))
dfparne5j_yyxy$bestsign = paste(dfparne5j_yyxy$best, dfparne5j_yyxy$effect, sep=' ') 



pot1a <- ggplot(dfparne3j) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparne3j_xx, aes(x=LONGITUDE, y=LATITUDE,
                  fill = effect, shape=effect), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_manual(values = c('darkred', 'blue'))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
          geom_segment(
    x = -76, y = 34,
    xend = -65, yend = 34,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkgreen" # Also accepts "red", "blue' etc
  ) + 
  annotate(geom="text", x=-72, y= 36, label="RIDGE", color="darkgreen", size = 5)+
  annotate(geom="text", x=-72, y=28, label="a) LOND model", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="eta[x]~aey==2", color="black", size=5, parse=TRUE) +
      # annotate(geom="text", x=-72, y=28, label="Wet~fraction", color="black", size=4, parse=TRUE) +
  # annotate(geom="text", x=-74, y=26, label="beta[y]^(MC)", color="black", size=4, parse=TRUE) +
    theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

pot1b <- ggplot(dfparne3j) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparne3j_yyxy, aes(x=LONGITUDE, y=LATITUDE,
                  fill = bestsign, shape=bestsign), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24, 25, 24))+
scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
            geom_segment(
    x = -66, y = 34,
    xend = -66, yend = 40,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkgreen" # Also accepts "red", "blue' etc
  ) + 
  annotate(geom="text", x=-70, y= 36, label="RIDGE", color="darkgreen", size = 5)+
  annotate(geom="text", x=-72, y=28, label="b) LATD and LWLD", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="eta[y]~aey==2", color="black", size=5, parse=TRUE) +
        # annotate(geom="text", x=-72, y=28, label="Wet~fraction", color="black", size=4, parse=TRUE) +
  # annotate(geom="text", x=-74, y=26, label="beta[y]^(MC)", color="black", size=4, parse=TRUE) +
  #   annotate(geom="text", x=-72, y=26, label="(b)~Correlation", color="black", size=4, parse=TRUE) +
  # annotate(geom="text", x=-74, y=24, label="rho[y]^(MC)", color="black", size=4, parse=TRUE) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.position = c(0.20, 0.13), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

pot2a <- ggplot(dfparne5j) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparne5j_xx, aes(x=LONGITUDE, y=LATITUDE,
                  fill = effect, shape=effect), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_manual(values = c('darkred', 'blue'))+
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  annotate(geom="text", x=-72, y=28, label="c) LOND model", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="eta[x]~aey==4", color="black", size=5, parse=TRUE) +
      # annotate(geom="text", x=-72, y=28, label="Correlation", color="black", size=4, parse=TRUE) +
  # annotate(geom="text", x=-74, y=26, label="rho[y]^(MC)", color="black", size=4, parse=TRUE) +
    theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

pot2b <- ggplot(dfparne5j) +
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
      geom_point(aes(x=LONGITUDE, y=LATITUDE,
                  ), size= 2.0, color = 'black', fill = 'grey', shape=21) +
        geom_point(data = dfparne5j_yyxy, aes(x=LONGITUDE, y=LATITUDE,
                  fill = bestsign, shape=bestsign), size= 2.0, color = 'black') +
  scale_shape_manual(values=c(25, 24, 25, 24))+
scale_fill_manual(values=c( 'LWLD Negative'="blue", "LWLD Positive"="darkred", 
                            'LATD Negative'="lightblue", "LATD Positive"="orange")) +
    xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  annotate(geom="text", x=-72, y=28, label="d) LATD and LWLD", color="black", size=4) +
  annotate(geom="text", x=-72, y=26, label="eta[y]~aey==4", color="black", size=5, parse=TRUE) +
    # annotate(geom="text", x=-72, y=28, label="Correlation", color="black", size=4, parse=TRUE) +
  # annotate(geom="text", x=-74, y=26, label="rho[y]^(MC)", color="black", size=4, parse=TRUE) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.position = c(0.20, 0.13), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)), legend.title=element_blank())

  figpot <- grid.arrange(pot1a, pot1b, pot2a, pot2b, ncol=2)
    ggsave(file.path(outplot, sprintf('res_maps_%s.png', mymodel5)),
          width = 13, height = 9, dpi = 300, plot=figpot)





# pot1 <- ggplot(dfparne3j) +
#     geom_polygon(data = states, aes(x=long, y = lat, group = group), 
#                  color = 'black', fill = 'grey', alpha=0.7) +
#     geom_point(aes(x=LONGITUDE, y=LATITUDE,
#                   fill = best, shape=effect), size = 2.0, color = 'black') +
#   # scale_shape_manual(values=c(6, 1, 2))+
#   scale_shape_manual(values=c(25, 21, 24))+
# # scale_fill_manual(values=c("grey40", "darkgreen", "lightblue", "darkred")) +
# scale_fill_manual(values=c("NOD"="grey40", "LATD"="darkgreen", "LOND"="darkred", "LWLD" = "lightblue")) +
#   # guides(shape=FALSE,  # remove legend for triangles direction
#   #    color=guide_legend(title.position="top", title = 'Best model'))+ 
#   guides(shape=FALSE,  # remove legend for triangles direction
#      fill=guide_legend(title.position="top", title = 'Best model',
#      override.aes = list(shape = 21, size = 4)))+
#     xlab('Longitude')+
#     ylab('Latitude')+
#     theme_bw()+
#   annotate(geom="text", x=-70, y=26, label="eta[y]~aey==2", color="black", size=5, parse=TRUE) +
#     theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
#           legend.background = element_rect(fill=alpha(0.4)))
#     # ggsave(file.path(outplot, sprintf('res_maps_%s_%s_binary.png', 
#     #       myparam, mymodel)), width = 6.5, height = 4)
# 
# 
# 
# 
# pot2 <- ggplot(dfparne5j) +
#     geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
#     geom_point(aes(x=LONGITUDE, y=LATITUDE,
#                    fill = best, shape=effect), size = 2.0, color = 'black') +
#   # scale_shape_manual(values=c(6, 1, 2))+
#   scale_shape_manual(values=c(25, 21, 24))+
#     # scale_fill_manual(values=c("grey50", "darkgreen", "lightblue", "darkred")) +
# scale_fill_manual(values=c("NOD"="grey40", "LATD"="darkgreen", "LOND"="darkred", "LWLD" = "lightblue")) +
#     # guides(shape=FALSE, color=FALSE)+ # remove legend for triangles direction
#     #  # color=guide_legend(title.position="top", title = 'Best model'))+ 
#   guides(shape=FALSE,  # remove legend for triangles direction
#      fill=guide_legend(title.position="top", title = 'Best model',
#      override.aes = list(shape = 21, size = 4)))+
#     xlab('Longitude')+
#     ylab('Latitude')+
#     theme_bw()+
#   # annotate(geom="text", x=-70, y=26, label="AEY = 4", color="black", size=5) +
#   annotate(geom="text", x=-70, y=26, label="eta[y]~aey==4", color="black", size=5, parse=TRUE) +
#     theme(legend.position = c(0.28, 0.1), legend.direction = "horizontal",
#           legend.background = element_rect(fill=alpha(0.4)))
#     # ggsave(file.path(outplot, sprintf('res_maps_%s_%s_binary.png', myparam, mymodel)),
#     #       width = 6.5, height = 4)
# 
#   figpot <- grid.arrange(pot1, pot2, ncol=2)
#     ggsave(file.path(outplot, sprintf('res_maps_pot_%s.png', mymodel)),
#           width = 13, height = 5, dpi = 300, plot=figpot)


    # ANNUAL MAXIMA - RESULTS FOR XYDEP ONLY
    
    

trvalue = 'qDqTR1'
maxmodel = 'nash_xij_weiNSO_xydep' 
dfmax <- subset(dfr, dfr$model == maxmodel & dfr$param == trvalue)

dfmax$Size = abs(dfmax$score)
dfmax$Sign = as.factor(sign(dfmax$score))
thr = 1
dfmax$effect = '1'
dfmax$effect[abs(dfmax$score) < thr] = 'Small'
dfmax$effect[dfmax$score > thr] = 'Positive'
dfmax$effect[dfmax$score < -thr] = 'Negative'

max1 <- ggplot(dfmax) +
  # geom_polygon(aes(x = long, y = lat, group = group), fill = "palegreen", color = "black") + 
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
    # geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill="grey", color = 'black', alpha=0.9) +
  # coord_quickmap()+
    geom_point(aes(x=LONGITUDE, y=LATITUDE,
                   color = effect, size =  Size, shape = Sign), alpha = 2) +
  # scale_shape_manual(values=c(6, 2))+
  # scale_color_gradient(low="blue", high="red")+
    scale_colour_manual(values = c("blue", "red", "black"))+
    guides(shape=FALSE, size=FALSE, # remove legend for triangles direction
    color = guide_legend(title.position="top", title = ''))+
    # guides(shape=FALSE, color=FALSE, size=FALSE)+ # remove legend for triangles direction
     # color=guide_legend(title.position="top", title = 'Best model'))+ 
    scale_size(range = c(0, 2))+
    xlab('Longitude')+
    ylab('Latitude')+
  scale_shape_manual(values=c(6, 2))+
    theme_bw()+
  annotate(geom="text", x=-120, y=27, label="Delta~q[Tr]/sigma[q[Tr]]", color="black", size=5, parse=TRUE) +
  annotate(geom="text", x=-72, y=26, label="a) Tr=20", color="black", size=5) +
    theme(legend.position = c(0.25, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)))

trvalue = 'qDqTR2'
maxmodel = 'nash_xij_weiNSO_xydep' 
dfmax <- subset(dfr, dfr$model == maxmodel & dfr$param == trvalue)

dfmax$Size = abs(dfmax$score)
dfmax$Sign = as.factor(sign(dfmax$score))
# thr = 1
dfmax$effect = '1'
dfmax$effect[abs(dfmax$score) < thr] = 'Small'
dfmax$effect[dfmax$score > thr] = 'Positive'
dfmax$effect[dfmax$score < -thr] = 'Negative'

max2 <- ggplot(dfmax) +
  # geom_polygon(aes(x = long, y = lat, group = group), fill = "palegreen", color = "black") + 
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey') +
    # geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill="grey", color = 'black', alpha=0.9) +
  # coord_quickmap()+
    geom_point(aes(x=LONGITUDE, y=LATITUDE,
                   color = effect, size =  Size, shape = Sign), alpha = 1) +
  # scale_shape_manual(values=c(6, 2))+
  # scale_color_gradient(low="blue", high="red")+
    scale_colour_manual(values = c("blue", "red", "black"))+
    # guides(shape=FALSE, color=FALSE, size=FALSE)+ # remove legend for triangles direction
      guides(shape=FALSE, size=FALSE, # remove legend for triangles direction
    color = guide_legend(title.position="top", title = ''))+
    scale_size(range = c(0, 2))+
    xlab('Longitude')+
    ylab('Latitude')+
  scale_shape_manual(values=c(6, 2))+
    theme_bw()+
  annotate(geom="text", x=-120, y=27, label="Delta~q[Tr]/sigma[q[Tr]]", color="black", size=5, parse=TRUE) +
  annotate(geom="text", x=-72, y=26, label="b) Tr=50", color="black", size=5) +
    theme(legend.position = c(0.25, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)))
#   
# figmax <- grid.arrange(max1, max2, ncol=1)
#     ggsave(file.path(outplot, sprintf('max_norm_over_width_%s.png', maxmodel)),
#           width = 6.5, height = 8.5, dpi = 300, plot=figmax)
#     
####################### to imrove ggplot color scales ##########################
    # log_both <- function(x){ifelse(x == 0, 0, log(abs(x)) * sign(x))}
    # exp_both <- function(x){exp(abs(x)) * sign(x)} # this is the inverse of log_both
    # 
    # log_both_trans <- 
    #   function(){
    #     trans_new(name = 'log_both', 
    #               transform = log_both,
    #               inverse = exp_both)}
    
trvalue = 'qnormTR1'
maxmodel = 'nash_xij_weiNSO_xydep' 
dfmax <- subset(dfr, dfr$model == maxmodel & dfr$param == trvalue)


# dfmaxord <- dfmax[order(dfmax$score),]

dfmax$Size = abs(dfmax$score)
dfmax$Sign = as.factor(sign(dfmax$score))
# thr = 1
# dfmax$effect = '1'
# dfmax$effect[abs(dfmax$score) < thr] = 'Small'
# dfmax$effect[dfmax$score > thr] = 'Positive'
# dfmax$effect[dfmax$score < -thr] = 'Negative'

# check max:
# dfmaxord = dplyr::arrange(dfmax, -score)
# dfmaxpos= dfmax[dfmax$score > 0]
# dfmaxneg= dfmax[dfmax$score < 0]

max3 <- ggplot(dfmax) +
  # geom_polygon(aes(x = long, y = lat, group = group), fill = "palegreen", color = "black") + 
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey') +
    # geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill="grey", color = 'black', alpha=0.9) +
  # coord_quickmap()+
    geom_point(aes(x=LONGITUDE, y=LATITUDE,
                    fill = score), size = 2, colour="black",pch=21 ) +
      # guides(fill=guide_legend(title.position="top", title = eval(parse(text="delta"))))+ 
      # guides(fill=guide_legend(title.position="top", title = ""))+ 


  # scale_fill_gradient2(low='red', high='blue', midpoint=0, mid='white')+ 
  # scale_fill_gradient2(low='red', high='blue', midpoint=0, mid='white')+ 
  scale_fill_gradient2(low='red', high='blue', name="")+ 
  # scale_shape_manual(values=c(6, 2))+
  # scale_color_gradient(low="blue", high="red")+
    # scale_colour_manual(values = c("blue", "red", "black"))+
    # scale_size(range = c(0, 2))+
    xlab('Longitude')+
    ylab('Latitude')+
  # scale_shape_manual(values=c(6, 2))+
  # scale_shape_manual(values=c(24, 25))+
  # scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  # scale_shape_manual(values=c("\u25BC","\u25B2"))+
    theme_bw()+
  annotate(geom="text", x=-72, y=26, label="c) Tr=20", color="black", size=5) +
  # annotate(geom="text", x=-120, y=27, label="Delta~q[Tr]/q[Tr]", color="black", size=5, parse=TRUE) +
  # annotate(geom="text", x=-123, y=26.5, label="Delta~q[Tr]/q[Tr]", color="black", size=5, parse=TRUE) +
  annotate(geom="text", x=-118, y=29, label="Delta~q[Tr]/q[Tr]", color="black", size=5, parse=TRUE) +
    theme(legend.position = c(0.18, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)))

#   grid.force()  # To make the grobs visible to grid editing tools
# grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 3))
# max1

trvalue = 'qnormTR2'
maxmodel = 'nash_xij_weiNSO_xydep' 
dfmax <- subset(dfr, dfr$model == maxmodel & dfr$param == trvalue)

dfmax$Size = abs(dfmax$score)
dfmax$Sign = as.factor(sign(dfmax$score))


dfmaxord = dplyr::arrange(dfmax, -score)
maxoutfilename <- file.path(outdata, 'dfmaxqnormTR2.csv')
write.csv(dfmaxord, maxoutfilename) # save results
# thr = 1
# dfmax$effect = '1'
# dfmax$effect[abs(dfmax$score) < thr] = 'Small'
# dfmax$effect[dfmax$score > thr] = 'Positive'
# dfmax$effect[dfmax$score < -thr] = 'Negative'

max4 <- ggplot(dfmax) +
  # geom_polygon(aes(x = long, y = lat, group = group), fill = "palegreen", color = "black") + 
    geom_polygon(data = states, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
    # geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill="grey", color = 'black', alpha=0.9) +
  # coord_quickmap()+
    geom_point(aes(x=LONGITUDE, y=LATITUDE,
                   fill = score), size = 2, colour="black",pch=21 ) +
  # scale_colour_gradient2(midpoint = 0)+
  # scale_fill_gradient2(low='red', high='blue', midpoint=0, mid='white')+ 
  scale_fill_gradient2(low='red', high='blue', guide="colourbar", name="")+ 
  # guide_colourbar(title = 'Hello', title.position = 'top')+
  
      # guides(fill=guide_legend(title.position="top", title = ''))+ 
  
  # scale_color_gradientn(colors = colorRampPalette(colors = c("blue", "white"))(nrow(dfmax)), 
  #                       values = scales::rescale(log(sort(dfmax$score))))+

    # scale_size(range = c(0, 2))+
    xlab('Longitude')+
    ylab('Latitude')+
  # scale_shape_manual(values=c(6, 2))+
  # scale_shape_manual(values=c("\u25BC","\u25B2"))+
    theme_bw()+
  annotate(geom="text", x=-72, y=26, label="d) Tr=50", color="black", size=5) +
  annotate(geom="text", x=-118, y=29, label="Delta~q[Tr]/q[Tr]", color="black", size=5, parse=TRUE) +
    theme(legend.position = c(0.18, 0.1), legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)))
  
figmax <- grid.arrange(max1, max3, max2, max4, ncol=2)
    ggsave(file.path(outplot, sprintf('max_norm_diff_%s.png', maxmodel)),
          width = 13, height = 9, dpi = 300, plot=figmax)



# load nash data

# ynash = read.table(file.path('nash_pos_data', 'ridge_lat_data.txt'), header=FALSE, col.names = c('YEAR', 'LAT'))
# xnash = read.table(file.path('nash_pos_data', 'ridge_lon_data.txt'), header=FALSE, col.names = c('YEAR', 'LON'))
# xnash$LON = -xnash$LON
# # plot(xnash$LON, ynash$LAT)
# 
# dfnash = xnash
# dfnash$LAT = ynash$LAT

    # dfnash = 
# min(xnash)

# southeast<- subset(states, region %in% c("florida", "georgia", "louisiana", 
#                                          "mississippi", "alabama", "south carolina",
#                                          "north carolina", "tennessee", "arkansas"
#                                          ))



mapplot1 <- ggplot(data = world) +
  geom_sf()+

# world = map_data("world")
# ggplot(dfnash) +
  # geom_polygon(aes(x = long, y = lat, group = group), fill = "palegreen", color = "black") + 
    # geom_polygon(data = world, aes(x=long, y = lat, group = group), color = 'black', fill = 'grey', alpha=0.7) +
  # geom_sf() +
  coord_sf(xlim = c(-100, -65), ylim = c(20, 38), expand = FALSE)+
  # coord_map(xlim = c(-98, -65),ylim = c(24, 40))+
  # geom_point(data=dfnash, aes(x=LON, y=LAT, color=YEAR), size=2)+
  geom_hline(yintercept=mean(dfnash$LAT), linetype = 'dashed')+
  geom_vline(xintercept=mean(dfnash$LON), linetype = 'dashed')+
      xlab('Longitude')+
    ylab('Latitude')+
  # xlim(-100, -60)+
  # ylim(25, 38)+
    theme_bw()+
    geom_point(data=dfnash, aes(x=LON, y=LAT, color = YEAR), size = 3)+
    theme(legend.position = c(0.85, 0.75),
          legend.background = element_rect(fill=alpha(0.4)))+
    ggsave(file.path(outplot, 'nashpos.png'),
          width = 6.5, height = 4.5)

mny = mean(dfnash$LAT)
mnx = mean(dfnash$LON)
sny = sd(dfnash$LAT)
snx = sd(dfnash$LON)
arrowsmap <- ggplot(data = world) +
  geom_sf()+
  coord_sf(xlim = c(-100, -65), ylim = c(20, 38), expand = FALSE)+
  geom_hline(yintercept=mnx, linetype = 'dashed')+
  geom_vline(xintercept=mny, linetype = 'dashed')+
      xlab('Longitude')+
    ylab('Latitude')+
    theme_bw()+
  
  geom_segment(
    x = mnx, y = mny,
    xend = mnx, yend = mny + 3*sny,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 2, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "lightblue" # Also accepts "red", "blue' etc
  ) + 
    geom_segment(
    x = mnx + 0.3, y = mny,
    xend = mnx + 0.3, yend = mny + 3*sny,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 2, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkgreen" # Also accepts "red", "blue' etc
  ) + 
    geom_segment(
    x = mnx, y = mny,
    xend = mnx + 1*snx, yend = mny,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 2, 
    arrow = arrow(length = unit(0.3, "inches")),
    colour = "darkred" # Also accepts "red", "blue' etc
  ) + 
  annotate(geom="text", x=mnx + snx, y=mny + sny, label="LOND", color="darkred", size = 7)+
  annotate(geom="text", x=mnx + snx/4, y=mny + 3.5*sny, label="LWLD", color="lightblue", size = 7)+
  annotate(geom="text", x=mnx + snx/4, y=mny + 4.2*sny, label="LATD", color="darkgreen", size = 7)+
  
    geom_point(data=dfnash, aes(x=LON, y=LAT, color = YEAR), size = 3, alpha = 0.8)+
    theme(legend.position = c(0.85, 0.75),
          legend.background = element_rect(fill=alpha(0.4)))
    # ggsave(file.path(outplot, 'nashpos.png'),
          # width = 6.5, height = 4.5)

  # figpar <- grid.arrange(par1, par2, par3, arrowsmap, ncol=2)
  #   ggsave(file.path(outplot, sprintf('nash_param_ncs.png')),
  #         width = 13, height = 9, dpi = 300, plot=figpar)
  # 
  # 



#