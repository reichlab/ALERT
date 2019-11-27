
source("inst/SIM-paper-code/ALERT-simulations-fxns.R")

library("dplyr")
library("xtable")
require(lubridate); year <- lubridate::year


dates <- read.csv("inst/SIM-paper-data/chco-trigger-dates.csv", stringsAsFactors = F)

dates$startDate <- ymd(dates$startDate)

dates$endDate <- ymd(dates$endDate)

##alternatively, use the following.
#real data
##load the data
fluA <- read.csv("inst/SIM-paper-data/CHCO-fluA.csv", 
                 stringsAsFactors=F)
fluA$Date <- ymd(fluA$Date)
fluA <- arrange(fluA, Date)

dates <- ymd(fluA$Date) ##need this later for the simulations

fluB <- read.csv("inst/SIM-paper-data/CHCO-fluB.csv", 
                 stringsAsFactors=F)
fluB$Date <- ymd(fluB$Date)
fluB <- arrange(fluB, Date)

RSV <- read.csv("inst/SIM-paper-data/CHCO-RSV.csv", 
                stringsAsFactors=F)
RSV$Date <- ymd(RSV$Date)
RSV <- arrange(RSV, Date)

Cases <- fluA$Cases + fluB$Cases + RSV$Cases
Date <- fluA$Date

data <- data.frame(Date, Cases)
all <- data



#bugnames <- colnames(all[2:4])

#all <- all[1:4]

#all$date <- mdy(all$Date..Month.Year.)

all <- arrange(all, by=Date)

all$state <- 0

#specify the model for endemic and autoregressive components
f.end <- surveillance::addSeason2formula(f = ~ 1+t, S=1, period=52)
model1 <- list(ar = list(f = ~ 1), end = list(f =f.end),
               family = "NegBinM", subset=2:469)

DisProg.dat <- surveillance::create.disProg(week = 1:nrow(all), 
                             observed = as.matrix(all[2]),
                             state = all$state, start = c(2001, 35))# convert to sts class
chco.dat <- surveillance::disProg2sts(DisProg.dat)# convert to sts class

#  # run model
chco.dat_res <- surveillance::hhh4(chco.dat, model1)
#  j <- as.integer(i-1)
holder <- chco.dat_res$coefficients


all <- data.frame(data.table::rbindlist(lapply(holder,as.list)))

#coefnames <- names(holder)
rownames(all) <- c("lambda", "alpha", "beta",
                   "gamma", "delta", "psi")

all[,2] <- data.frame(data.table::rbindlist(lapply(chco.dat_res$se,as.list)))

colnames(all) <- c("Estimated coefficients", "SE")


######################
## TABLE 2 in paper ##
######################

print(xtable(all, type="html", digits=4))

### use param values within 1 sd for upper and lower bounds
upper <- all[1]+(2*chco.dat_res$se)
lower <- all[1]-(2*chco.dat_res$se)

#upper <- rep(1.75, nrow(all))
#lower <- rep(-3.0, nrow(all))

limits <- data.frame(lower, upper)

#seq(from=limits[1,1], to=limits[1,2], length.out=50)

holder <- list()

for (i in 1:nrow(limits)){
  holder[[i]] <- seq(from=limits[i,1], to=limits[i,2], length.out=50
                     )
}

ranges <- t(data.frame(data.table::rbindlist(lapply(holder,as.list))))

names(holder) <- rownames(all) 

ranges <- data.frame(ranges)

#holder <- list()

#names(holder) <- coefnames

coefs <- chco.dat_res$coefficients
colnames(ranges) <- rownames(all) 

#######################################
########   GET PARAMS READY!!  ########
#######################################


###just change one param at a time.

## first end.1
alpha <- ranges
alpha$ar.1 <- coefs[1]  ## not going to mess with the AR param lambda
alpha$end.t <- coefs[3]  ###or the slope beta
alpha[4] <- coefs[4] ##or this seasonal term gamma
alpha[5] <- coefs[5] ##delta 
alpha[6] <- coefs[6] ##psi
alpha$param <- "alpha"

## overdisp
psi <- ranges
psi$ar.1 <- coefs[1]  ## not going to mess with the AR param lambda
psi$end.1 <- coefs[2]  ## alpha
psi$end.t <- coefs[3]  ###or the slope beta
psi[4] <- coefs[4] ##or this seasonal term gamma
psi[5] <- coefs[5] ##delta 
psi$param <- "psi"

## the delta cos term
delta <- ranges
delta$ar.1 <- coefs[1]  ## not going to mess with the AR param lambda
delta$end.1 <- coefs[2]  ## alpha
delta$end.t <- coefs[3]  ###or the slope beta
delta[4] <- coefs[4] ##or this seasonal term gamma
delta[6] <- coefs[6] ##psi
delta$param <- "delta"

##lambda
lambda <- ranges
lambda$end.1 <- coefs[2]  ## alpha
lambda$end.t <- coefs[3]  ###or the slope beta
lambda[4] <- coefs[4] ##or this seasonal term gamma
lambda[5] <- coefs[5] ##delta 
lambda[6] <- coefs[6] ##psi
lambda$param <- "lambda"

##beta
beta <- ranges
beta$ar.1 <- coefs[1]  ## not going to mess with the AR param lambda
beta$end.1 <- coefs[2]  ## alpha
beta[4] <- coefs[4] ##or this seasonal term gamma
beta[5] <- coefs[5] ##delta 
beta[6] <- coefs[6] ##psi
beta$param <- "beta"

##gamma
gamma <- ranges
gamma$ar.1 <- coefs[1]  ## not going to mess with the AR param lambda
gamma$end.1 <- coefs[2]  ## alpha
gamma$end.t <- coefs[3]  ###or the slope beta
gamma[5] <- coefs[5] ##delta 
gamma[6] <- coefs[6] ##psi
gamma$param <- "gamma"

ranges <- bind_rows(lambda, alpha, beta, delta, gamma, psi)

##########################################
######## PREP FOR SIMULATION #############
##########################################

#dates <- chco$Date

library("ALERT")

#this is for when there is a firstMonth error. It will be filtered in the following steps.
filler <- (createALERT(data, firstMonth=9))

#choose the number of simulations to perform

snum <- 20

simulations <- list()

#create the empty data.frame to allocate memory.

num <- (nrow(ranges)*snum)

alertstats2 <- data.frame(threshold=rep(NA, num),
                          median.dur=rep(NA, num),
                          median.pct.cases.captured=rep(NA, num),
                          min.pct.cases.captured=rep(NA, num),
                          max.pct.cases.captured=rep(NA, num),
                          pct.peaks.captured=rep(NA, num),
                          pct.ext.peaks.captured=rep(NA, num),
                          mean.low.weeks.incl=rep(NA, num),
                          parameter=rep(NA, num),
                          value=rep(NA, num),
                          success=rep(NA, num),
                          firstMonth=rep(NA, num),
                          sim_number=rep(NA, num))



##Make a column with the changed values

ranges$value <- NA

ranges$value[ranges$param=="lambda"] <- ranges$ar.1[ranges$param=="lambda"]
ranges$value[ranges$param=="alpha"] <- ranges$end.1[ranges$param=="alpha"]
ranges$value[ranges$param=="beta"] <- ranges$end.t[ranges$param=="beta"]
ranges$value[ranges$param=="delta"] <- ranges$`end.cos(2 * pi * t/52)`[ranges$param=="delta"]
ranges$value[ranges$param=="gamma"] <- ranges$`end.sin(2 * pi * t/52)`[ranges$param=="gamma"]
ranges$value[ranges$param=="psi"] <- ranges$`-log(overdisp)`[ranges$param=="psi"]

########################################
##########  SIMULATE DATA!  ############
########################################

for (i in 1:nrow(ranges)){
  print(i)
  params <- c(ranges$param[i], ranges$value[i])
  toycoef <- ranges[i,1:6]
    res2 <- surveillance::hhh4(chco.dat, model1)
    res2$coefficients <- unlist(c(toycoef))
     simulation <- as.vector(simulate(res2, nsim=snum, seed=123,
                                     y.start=NULL, simplify=TRUE))
     minisim <- split(simulation, as.factor(rep(seq(1,snum,1), each=length(dates))))
     for(k in 1:snum){
      simset <- data.frame(minisim[[k]], dates)
      colnames(simset) <- c("Cases", "Date")
      simset$Date <- as.Date(simset$Date)
      write.csv(simset, file=paste("inst/SIM-paper-simdata/", 
                                   ranges$param[i], ranges$value[i],
                                   "sim#", k, ".csv", sep="_"))
      simset_train <- simset[1:260,]
      simset_test <- simset[261:nrow(simset),]
      ##to get the FirstMonth
      gamma <- as.numeric(res2$coefficients[4])
      delta <- as.numeric(res2$coefficients[5])
      end.fxn <- function(y) gamma*sin(2*pi*y) + delta*cos(2*pi*y)
      xmax <- optimize(f=end.fxn, lower=0, upper=1, maximum=TRUE)$maximum
      xmin <- optimize(f=end.fxn, lower=0, upper=1, maximum=FALSE)$minimum
      firstMonth_value <- month(simset_train$Date[round(uniroot(end.fxn, c(0,1), 
                                                                 extendInt="upX", trace=TRUE)$root*52)])
      minWeeks_value <- 8
      alertstats2[((i-1)*snum)+k,] <- c(get_stats(simset_train, firstMonth=firstMonth_value,
                                                  minWeeks=minWeeks_value,
                                          params=params), firstMonth_value, k)
     # print(alertstats2[((i-1)*snum)+k,])
     }
}

alertstats3 <- filter(alertstats2, !is.na(threshold))

success_num <- nrow(alertstats3)

alerttest <- data.frame(test.threshold=rep(NA, success_num),
                        test.median.dur=rep(NA, success_num),
                        test.median.pct.cases.captured=rep(NA, success_num),
                        test.min.pct.cases.captured=rep(NA, success_num),
                        test.max.pct.cases.captured=rep(NA, success_num),
                        test.pct.peaks.captured=rep(NA, success_num),
                        test.pct.ext.peaks.captured=rep(NA, success_num), 
                        test.mean.low.weeks.incl=rep(NA, success_num))

for (i in 1:nrow(alertstats3)){
  alerttest[i,] <- thresholdtestALERT(simset_test,
  whichThreshold=as.numeric(alertstats3[i,1]), firstMonth=firstMonth_value,
  minWeeks = minWeeks_value)$out
  print(alerttest[i,])
}

alertstats <- data.frame(alertstats3, alerttest)

if(nrow(filter(alertstats, threshold!=test.threshold))>0){
  stop("alertstats dataframe construction ERROR")
}

med.stats <- alertstats %>% mutate(threshold=as.numeric(threshold), 
                        train.median.dur=as.numeric(median.dur),
                        median.pct.cases.captured = as.numeric(median.pct.cases.captured),
    test.median.pct.cases.captured=as.numeric(test.median.pct.cases.captured), 
    test.median.dur=as.numeric(test.median.dur),
    test.mean.low.weeks.incl=as.numeric(test.mean.low.weeks.incl),
    mean.low.weeks.incl=as.numeric(mean.low.weeks.incl),
    firstMonth=as.numeric(firstMonth)) %>% 
  filter(success==TRUE) %>% group_by(parameter, value) %>%
    summarise(test.median.pct.cases.captured=median(test.median.pct.cases.captured), 
            threshold=median(threshold),
            test.median.dur=median(test.median.dur), 
            train.median.dur=median(train.median.dur),
            median.pct.cases.captured = median(median.pct.cases.captured),
            median.test.low.weeks=median(test.mean.low.weeks.incl),
            median.train.low.weeks=median(mean.low.weeks.incl),
            median.firstMonth=median(firstMonth),
            test.peaks=median(test.pct.peaks.captured, na.rm=T),
            train.peaks=median(pct.peaks.captured, na.rm=T)) %>%
  mutate(duration.diff=train.median.dur-test.median.dur,
         median.pct.diff=median.pct.cases.captured-test.median.pct.cases.captured,
         median.low.weeks.diff=median.train.low.weeks-median.test.low.weeks,
         peaks.pct.diff=train.peaks-test.peaks) %>%
  data.frame()

##alertstats2 holds the uncollapsed ALERT results
## med.stats has the collapsed results
write.csv(med.stats, "inst/SIM-paper-simdata/med_stats.csv")
write.csv(alertstats, "inst/SIM-paper-simdata/alertstats.csv")
write.csv(alertstats2, "inst/SIM-paper-simdata/alertstats2.csv")

med.stats <- read.csv("inst/SIM-paper-simdata/med_stats.csv", stringsAsFactors = F)
alertstats <- read.csv("inst/SIM-paper-simdata/alertstats.csv", stringsAsFactors = F)
alertstats2 <- read.csv("inst/SIM-paper-simdata/alertstats2.csv", stringsAsFactors = F)

#it looks like ar.1, end.1, and end.t are going to be the most interesting.
#maybe overdispersion?

#plot this bit with each parameter adjusted while holding the others constant.

#I'm filtering the last row (overdispersion) because the threshold is crazy high
#and medians are either NA or 0.

#med.stats <- med.stats[-60,]

library("ggplot2")
#library("gridExtra")

alpha.plot <- ggplot(filter(alertstats2, parameter=="alpha")) +
  geom_point(aes(x=as.numeric(value), y=median.pct.cases.captured, color=success)) +
##  scale_color_gradient2(guide=guide_colorbar(direction="vertical"),
  #                     limits=c(3, (8)), low = "white", mid = "black", high = "black", midpoint = 6, na.value='grey') +
  scale_size(range = c(1, 6)) +
  #facet_wrap(~threshold, ncol=3)+
  theme_classic() +
  geom_hline(aes(yintercept=80), linetype="longdash", show.legend=FALSE) +
  xlab("parameter value")+
  ylab("median percent cases captured")

delta.plot <- ggplot(filter(alertstats2, parameter=="delta")) +
  geom_point(aes(x=as.numeric(value), y=median.pct.cases.captured, color=success)) +
  ##  scale_color_gradient2(guide=guide_colorbar(direction="vertical"),
  #                     limits=c(3, (8)), low = "white", mid = "black", high = "black", midpoint = 6, na.value='grey') +
  scale_size(range = c(1, 6)) +
  #facet_wrap(~threshold, ncol=3)+
  theme_classic() +
  geom_hline(aes(yintercept=80), linetype="longdash", show.legend=FALSE) +
  xlab("parameter value")+
  ylab("median percent cases captured")

lambda.plot <- ggplot(filter(alertstats2, parameter=="lambda")) +
  geom_point(aes(x=as.numeric(value), y=median.pct.cases.captured, color=success)) +
  ##  scale_color_gradient2(guide=guide_colorbar(direction="vertical"),
  #                     limits=c(3, (8)), low = "white", mid = "black", high = "black", midpoint = 6, na.value='grey') +
  scale_size(range = c(1, 6)) +
  #facet_wrap(~threshold, ncol=3)+
  theme_classic() +
  geom_hline(aes(yintercept=80), linetype="longdash", show.legend=FALSE) +
  xlab("parameter value")+
  ylab("median percent cases captured")

beta.plot <- ggplot(filter(alertstats2, parameter=="beta")) +
  geom_point(aes(x=as.numeric(value), y=median.pct.cases.captured, color=success)) +
  ##  scale_color_gradient2(guide=guide_colorbar(direction="vertical"),
  #                     limits=c(3, (8)), low = "white", mid = "black", high = "black", midpoint = 6, na.value='grey') +
  scale_size(range = c(1, 6)) +
  #facet_wrap(~threshold, ncol=3)+
  theme_classic() +
  geom_hline(aes(yintercept=80), linetype="longdash", show.legend=FALSE) +
  xlab("parameter value")+
  ylab("median percent cases captured")

gamma.plot <- ggplot(filter(alertstats2, parameter=="gamma")) +
  geom_point(aes(x=as.numeric(value), y=median.pct.cases.captured, color=success)) +
  ##  scale_color_gradient2(guide=guide_colorbar(direction="vertical"),
  #                     limits=c(3, (8)), low = "white", mid = "black", high = "black", midpoint = 6, na.value='grey') +
  scale_size(range = c(1, 6)) +
  #facet_wrap(~threshold, ncol=3)+
  theme_classic() +
  geom_hline(aes(yintercept=80), linetype="longdash", show.legend=FALSE) +
  xlab("parameter value")+
  ylab("median percent cases captured")

psi.plot <- ggplot(filter(alertstats2, parameter=="psi")) +
  geom_point(aes(x=as.numeric(value), y=median.pct.cases.captured, color=success)) +
  ##  scale_color_gradient2(guide=guide_colorbar(direction="vertical"),
  #                     limits=c(3, (8)), low = "white", mid = "black", high = "black", midpoint = 6, na.value='grey') +
  scale_size(range = c(1, 6)) +
  #facet_wrap(~threshold, ncol=3)+
  theme_classic() +
  geom_hline(aes(yintercept=80), linetype="longdash", show.legend=FALSE) +
  xlab("parameter value")+
  ylab("median percent cases captured")

#The sin term manages the width of the seasonal epidemic/peak height/baseline noise

#cos is the width of the sinus, Larger cos term should make the sinus more narrow. 
#smaller makes the sinus wider (less baseline noise)
med.stats$parameter2 <- factor(med.stats$parameter, 
                              labels = c(unique(med.stats$parameter)))

med.stats$value <- as.numeric(med.stats$value)

med.perc.diff.performance <- ggplot(med.stats, group=parameter2) +
 # geom_point(aes(x=value, y=median.low.weeks.diff, color=threshold)) +
  #scale_colour_gradient2(guide=guide_colorbar(direction="vertical"),
   #                      limits=c(min(med.stats$threshold), (max(med.stats$threshold))), 
    #                     low = "green", mid = "purple", high = "purple", 
     #                    midpoint = max(med.stats$threshold)-1.5, na.value='grey') +
  geom_smooth(aes(x=value, y=test.median.pct.cases.captured), colour="black") +
  facet_wrap(~parameter2, ncol=2, scales = "free_x", labeller = label_parsed)+
  theme_classic() +
  geom_hline(aes(yintercept=0), linetype="longdash", show.legend=FALSE, alpha=0.35) +
  xlab("parameter value")+
  ylab("median percent cases captured difference")

med.perc.diff.performance

med.lowweeks.diff.performance <- ggplot(med.stats, group=parameter2) +
  # geom_point(aes(x=value, y=median.low.weeks.diff, color=threshold)) +
  #scale_colour_gradient2(guide=guide_colorbar(direction="vertical"),
  #                      limits=c(min(med.stats$threshold), (max(med.stats$threshold))), 
  #                     low = "green", mid = "purple", high = "purple", 
  #                    midpoint = max(med.stats$threshold)-1.5, na.value='grey') +
  geom_smooth(aes(x=value, y=median.low.weeks.diff), colour="black") +
  facet_wrap(~parameter2, ncol=2, scales = "free_x", labeller = label_parsed)+
  theme_classic() +
  geom_hline(aes(yintercept=0), linetype="longdash", show.legend=FALSE, alpha=0.35) +
  xlab("parameter value")+
  ylab("median low weeks captured difference")

med.lowweeks.diff.performance

peaks.pct.diff.performance <- ggplot(alertstats, group=parameter) +
  # geom_point(aes(x=value, y=median.low.weeks.diff, color=threshold)) +
  #scale_colour_gradient2(guide=guide_colorbar(direction="vertical"),
  #                      limits=c(min(med.stats$threshold), (max(med.stats$threshold))), 
  #                     low = "green", mid = "purple", high = "purple", 
  #                    midpoint = max(med.stats$threshold)-1.5, na.value='grey') +
  geom_smooth(aes(x=value, y=test.pct.peaks.captured), colour="black") +
  facet_wrap(~parameter, ncol=2, scales = "free_x", labeller = label_parsed)+
  theme_classic() +
#  geom_hline(aes(yintercept=0), linetype="longdash", show.legend=FALSE, alpha=0.35) +
  xlab("parameter value")+
  ylab("median percent peaks captured")

peaks.pct.diff.performance

######################################
###### PLOT SOME FAKE DATA!!!! #######
######################################

## plot a set of randomly selected data

directory <- dir("inst/SIM-paper-simdata/")

sim_plot <- sample(directory,20,replace=FALSE)

start_and_end_real <- list()
selected_sims_real <- list()
sims_metadata_real <- list()

for (i in 1:length(sim_plot)){
  print(sim_plot[i])
  holder <- read.csv(paste("inst/SIM-paper-simdata/", sim_plot[i], sep=''))
  sim_components <- unlist(strsplit(sim_plot[i], "_"))
##extract the parameter info and simulation number  
  parameter_realized <- sim_components[2]
  value_realized <- sim_components[3]
  sim_number_realized <- as.integer(sim_components[5])
## set firstMonth and threshold to what was determined to be optimal
  ## to prep for thresholdtestALERT
  first_Month <- filter(alertstats2, parameter==parameter_realized &
                          value==value_realized & 
                          sim_number==sim_number_realized)$firstMonth
  threshold <- filter(alertstats2, parameter==parameter_realized &
                              value==value_realized & 
                              sim_number==sim_number_realized)$threshold
##make sure Dates are Dates
  holder$Date <- as.Date(holder$Date)
  ALERT_dates <- select(data.frame(
    thresholdtestALERT(holder, whichThreshold = threshold, 
                       firstMonth = first_Month)$details), start, end)
##convert dates to something human readable for plotting
  ALERT_dates$start <- as.Date(ALERT_dates$start, origin="1970-01-01")
  ALERT_dates$end <- as.Date(ALERT_dates$end, origin="1970-01-01")
  colnames(ALERT_dates) <- c("xstart", "xstop")
##replace case counts>100 with 100 to keep them from being dropped or needing to rescale
  holder$Cases[holder$Cases>100] <- 100
  start_and_end_real[[i]] <- ALERT_dates
  selected_sims_real[[i]] <- holder  
  sims_metadata_real[[i]] <- c(parameter_realized, value_realized, sim_number_realized, first_Month, threshold)
}

##for example

sim_compare_figs(start_and_end_real, selected_sims_real, 6, sims_metadata = sims_metadata_real)

library(gridExtra)

grid.arrange(sim_compare_figs(start_and_end_real, selected_sims_real, 3, sims_metadata_real), 
             sim_compare_figs(start_and_end_real, selected_sims_real, 5, sims_metadata_real), 
             sim_compare_figs(start_and_end_real, selected_sims_real, 7, sims_metadata_real),
             sim_compare_figs(start_and_end_real, selected_sims_real, 1, sims_metadata_real), 
             sim_compare_figs(start_and_end_real, selected_sims_real, 10, sims_metadata_real),
             sim_compare_figs(start_and_end_real, selected_sims_real, 2, sims_metadata_real),
             ncol=1, left="simulated cases", 
             bottom="Date")


###################################################
######### TIME TO ANALYZE THE RESULTS  ############
###################################################

##compare training to test statistics
require(xtable)

summary.stats <- alertstats %>% group_by(parameter) %>% summarize(
  min(median.dur, na.rm=T),
  max(median.dur, na.rm=T),
  median(median.dur), 
  min(test.median.dur, na.rm=T),
  max(test.median.dur, na.rm=T),
  median(test.median.dur), 
  min(as.numeric(median.pct.cases.captured), na.rm=T),
  max(as.numeric(median.pct.cases.captured), na.rm=T),
  median(as.numeric(median.pct.cases.captured)), 
  min(test.median.pct.cases.captured, na.rm=T),
  max(test.median.pct.cases.captured, na.rm=T),
  median(test.median.pct.cases.captured)) %>% data.frame

rownames(summary.stats) <- summary.stats$parameter
summary.stats$parameter <- NULL

summary.stats.overall <- alertstats %>% summarize(
  min(median.dur, na.rm=T),
  max(median.dur, na.rm=T),
  median(median.dur), 
  min(test.median.dur, na.rm=T),
  max(test.median.dur, na.rm=T),
  median(test.median.dur), 
  min(as.numeric(median.pct.cases.captured), na.rm=T),
  max(as.numeric(median.pct.cases.captured), na.rm=T),
  median(as.numeric(median.pct.cases.captured)), 
  min(test.median.pct.cases.captured, na.rm=T),
  max(test.median.pct.cases.captured, na.rm=T),
  median(test.median.pct.cases.captured)) %>% data.frame


summary_holder1 <- list()
summary_holder2 <- list()
summary_holder3 <- list()
summary_holder4 <- list()
for (i in 1:nrow(summary.stats)){
  summary_holder1[i] <- print(table_pretty(summary.stats[i,1], summary.stats[i, 2], summary.stats[i, 3]))
}
for (i in 1:nrow(summary.stats)){
  summary_holder2[i] <- print(table_pretty(summary.stats[i,4], summary.stats[i, 5], summary.stats[i, 6]))
}
for (i in 1:nrow(summary.stats)){
  summary_holder3[i] <- print(table_pretty(summary.stats[i,7], summary.stats[i, 8], summary.stats[i, 9]))
}
for (i in 1:nrow(summary.stats)){
  summary_holder4[i] <- print(table_pretty(summary.stats[i,10], summary.stats[i, 11], summary.stats[i, 12]))
}
tab <- cbind(summary_holder1, summary_holder2, summary_holder3, summary_holder4)

colnames(tab) <- c("training; duration (weeks)", "testing; duration (weeks)", 
                "training; cases (%)", "testing; cases (%)")

rownames(tab) <- rownames(summary.stats)

print(xtable(tab), include.rownames=T)

##compare low weeks captured between training and testing datasets

lowweek.perform.stats <- alertstats %>% group_by(parameter) %>% summarize(
  min(as.numeric(mean.low.weeks.incl), na.rm=T),
  max(as.numeric(mean.low.weeks.incl), na.rm=T),
  median(as.numeric(mean.low.weeks.incl), na.rm=T),
  min(test.mean.low.weeks.incl, na.rm=T),
  max(test.mean.low.weeks.incl, na.rm=T),
  median(test.mean.low.weeks.incl, na.rm=T),
  min(pct.peaks.captured, na.rm=T), 
  max(pct.peaks.captured, na.rm=T), 
  median(pct.peaks.captured, na.rm=T), 
  min(test.pct.peaks.captured, na.rm=T),
  max(test.pct.peaks.captured, na.rm=T),
  median(test.pct.peaks.captured, na.rm=T)) %>% data.frame

lowweek.perform.stats.overall <- alertstats %>% summarize(
  min(as.numeric(mean.low.weeks.incl), na.rm=T),
  max(as.numeric(mean.low.weeks.incl), na.rm=T),
  median(as.numeric(mean.low.weeks.incl), na.rm=T),
  min(test.mean.low.weeks.incl, na.rm=T),
  max(test.mean.low.weeks.incl, na.rm=T),
  median(test.mean.low.weeks.incl, na.rm=T),
  min(pct.peaks.captured, na.rm=T), 
  max(pct.peaks.captured, na.rm=T), 
  median(pct.peaks.captured, na.rm=T), 
  min(test.pct.peaks.captured, na.rm=T),
  max(test.pct.peaks.captured, na.rm=T),
  median(test.pct.peaks.captured, na.rm=T)) %>% data.frame

rownames(lowweek.perform.stats) <- lowweek.perform.stats$parameter
lowweek.perform.stats$parameter <- NULL

summary_holder1 <- list()
summary_holder2 <- list()
summary_holder3 <- list()
summary_holder4 <- list()
for (i in 1:nrow(lowweek.perform.stats)){
  summary_holder1[i] <- print(table_pretty(lowweek.perform.stats[i,1], lowweek.perform.stats[i, 2], lowweek.perform.stats[i, 3]))
}
for (i in 1:nrow(lowweek.perform.stats)){
  summary_holder2[i] <- print(table_pretty(lowweek.perform.stats[i,4], lowweek.perform.stats[i, 5], lowweek.perform.stats[i, 6]))
}
for (i in 1:nrow(lowweek.perform.stats)){
  summary_holder3[i] <- print(table_pretty(lowweek.perform.stats[i,7], lowweek.perform.stats[i, 8], lowweek.perform.stats[i, 9]))
}
for (i in 1:nrow(lowweek.perform.stats)){
  summary_holder4[i] <- print(table_pretty(lowweek.perform.stats[i,10], lowweek.perform.stats[i, 11], lowweek.perform.stats[i, 12]))
}
tab <- cbind(summary_holder1, summary_holder2, summary_holder3, summary_holder4)

colnames(tab) <- c("training; low weeks (weeks)", "testing; low weeks (weeks)", 
                   "training; peaks captured (%)", "testing; peaks captured (%)")

rownames(tab) <- rownames(lowweek.perform.stats)

print(xtable(tab), include.rownames=T)




