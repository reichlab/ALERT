##this is for the applied ALERT paper.

require(ALERT)
require(surveillance)
require(lubridate); year <- lubridate::year
require(ggplot2)
require(dplyr)

#real data
##load the data
#chco <- read.csv("~/Desktop/applied-ALERT-data/CHCO-fluA.csv", 
#                 stringsAsFactors=F)

##FLU A ONLY

chco <- fulldata
chco$Date <- ymd(chco$Date)
chco <- arrange(chco, Date)
#dates <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/chco-trigger-dates.csv", stringsAsFactors = F)
dates$startDate <- ymd(dates$startDate)
dates$endDate <- ymd(dates$endDate)
#truncate to get years we have cutoffs for and cutoffs we have years for.
chco <- filter(chco, Date >= head(dates$startDate, 1)-weeks(2))
dates <- filter(dates, endDate <= tail(chco$Date, 1)+weeks(6))


chco_train <- chco[1:210,]
chco_test <- chco[211:nrow(chco),]

train_dates <- dates[1:4,]




median(difftime(train_dates$startDate, train_dates$endDate))  # median of 18 weeks long in the real data
128.5/7


alertholder <- createALERT(chco_train, firstMonth=7, lag=0, allThresholds = T) #createALERT maybe not calculating median duration corectly?
#I slacked Steve about it.

threstest <- thresholdtestALERT(chco_train, firstMonth = 7, whichThreshold = 10)

ALERT_dates <- data.frame(threstest$details)
##convert dates to something human readable for plotting
ALERT_dates$start <- as.Date(ALERT_dates$start, origin="1970-01-01")
ALERT_dates$end <- as.Date(ALERT_dates$end, origin="1970-01-01")

median(difftime(ALERT_dates$start, ALERT_dates$end))  # median of 13.6 weeks long
mean(difftime(ALERT_dates$start, ALERT_dates$end))  # mean of 13 weeks long
91.5/7
173.25/7

##choose a threshold of 3 to get 24.75 weeks.



thres2 <- thresholdtestALERT(chco, firstMonth = 7, whichThreshold = 10)[[2]]

ALERT_dates <- data.frame(thres2)
##convert dates to something human readable for plotting
ALERT_dates$start <- as.Date(ALERT_dates$start, origin="1970-01-01")
ALERT_dates$end <- as.Date(ALERT_dates$end, origin="1970-01-01")

median(difftime(ALERT_dates$start, ALERT_dates$end))# median of 17 weeks long
94.7/7


ALERT_dates <- ALERT_dates[,8:9]
colnames(ALERT_dates) <- c("xstart", "xstop")

library(gridExtra)

grid.arrange(real_compare_figs(ALERT_dates, chco, "ALERT, threshold=2"),
             real_compare_figs(dates, chco, "actual dates"))


createALERT(chco_test)




median(difftime(dates$start, dates$end))
128.5/7


sum(difftime(dates$start, dates$end))
1088/7

threstest <- thresholdtestALERT(chco, firstMonth = 7, whichThreshold = 10)

ALERT_dates <- data.frame(threstest$details)
##convert dates to something human readable for plotting
ALERT_dates$start <- as.Date(ALERT_dates$start, origin="1970-01-01")
ALERT_dates$end <- as.Date(ALERT_dates$end, origin="1970-01-01")

sum(difftime(ALERT_dates$start, ALERT_dates$end))
1015/7







dates$yearIdx <- 1:nrow(dates)

vecs <- list()

for (i in 1:nrow(dates)){
  print(i)
  vecs[[i]] <- ifelse(dates$startDate[i] < chco$Date & dates$endDate[i] > chco$Date, TRUE, FALSE)
}

holder <- t(data.frame(Reduce(rbind, vecs)))

alertperiod <- list()

for (i in 1:nrow(holder)){
  alertperiod[[i]] <- any(holder[i,]==TRUE)
}

chco$alertperiod <- as.vector(Reduce(rbind, alertperiod))

#this doesn't work because flu doesn't care what year it is
#o$year <- factor(substr(as.character(o$Date), 1,4))

A <- rep(1, 42)

for(i in 2:7){
  A <- append(rep(i, 51), A)
}

A <- append(rep(8, nrow(chco)-length(A)), A)

chco$year <- as.factor(rev(A))


tot <- chco %>% group_by(year) %>% summarise(total_cases=sum(Cases)) #total cases per year

alertcases <- chco %>% filter(alertperiod==TRUE) %>% group_by(year) %>% 
  summarise(alert_cases=sum(Cases)) #"ALERT" cases per year

calc <- full_join(tot, alertcases)

calc[is.na(calc)] <- 0

calc$perc_alert <- calc$alert_cases/calc$total_cases *100 #percent alert cases captured

median(calc$perc_alert) #median cases captured

alert_dura <- chco %>% filter(alertperiod==TRUE) %>% group_by(year) %>% #duration
  summarise(duration=as.numeric(length(Date)))

calc <- unique(full_join(calc, alert_dura))

holder <- chco %>% group_by(year) %>% summarise(Cases=max(Cases)) %>%
  data.frame()

holder <- left_join(holder, chco)

holder <- holder %>% group_by(year) %>% summarise(Date=first(Date))

holder <- left_join(holder, chco) %>% data.frame()

calc$peak_captured <- holder$alertperiod

peaky <- table(calc$peak_captured) %>% data.frame() 

if(nrow(peaky)==2){
  perc_peaks_captured <- peaky$Freq[2] / peaky$Freq[1]*100
} else {
  if (peaky$Var1==TRUE){
    perc_peaks_captured <- 100
  } else {
    perc_peaks_captured <- 0
  }
}

holder <- chco %>% group_by(year) %>% summarise(Cases=min(Cases)) %>%
  data.frame()

holder <- left_join(holder, chco)

#number of 0 weeks included
holder <- filter(holder, alertperiod==TRUE) %>% group_by(year) %>% 
  summarise(low_weeks_incl=length(alertperiod))

calc <- left_join(calc, holder) %>% data.frame()

calc$year <- NULL

calc[is.na(calc)] <- 0

cutoffs <- dates %>% arrange(-yearIdx)

calc <- cbind(cutoffs, calc)

calc$yearIdx <- NULL

#get the colnames for eval statistics

alertholder <- data.frame(alertholder$out)

snames <- colnames(alertholder)[-7]

stats_real <- t(data.frame(c(NA, 
                             round(median(calc$duration), 1),
                             round(median(calc$perc_alert), 1), 
                             round(min(calc$perc_alert), 1), 
                             round(max(calc$perc_alert), 1), 
                             round(perc_peaks_captured, 3), 
                             round(mean(calc$low_weeks_incl), 1))))

colnames(stats_real) <- snames

comparison <- full_join(data.frame(stats_real), alertholder)

require(xtable)

print(xtable(comparison[1:7]), include.rownames=F)

str(stats_real)

















#plotting
#choose a threshold

threshold <- 3

chco$smooththres <- ifelse (lowess(chco$Cases, f=0.001)$y>=threshold, TRUE, FALSE)


ALERT_dates <- chco %>% group_by(year, smooththres) %>% arrange(Date) %>%
  summarise(xstop=last(Date)-1, xstart=first(Date)-1) %>% 
  data.frame() %>% 
  filter(smooththres==TRUE) %>% select(-year, -smooththres)

ymin <- -2

datetrig <- ggplot() + #this is the real dates
  theme_classic() +
  labs(y = "Influenza A cases", x = "Date-based intervention periods") +
  geom_bar(aes(y=chco$Cases, x=chco$Date), stat="identity") +
  geom_rect(aes(xmin = dates$startDate[1], 
                xmax = dates$endDate[1], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[2], 
                xmax = dates$endDate[2], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[3], 
                xmax = dates$endDate[3], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[4], 
                xmax = dates$endDate[4], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[5], 
                xmax = dates$endDate[5], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[6], 
                xmax = dates$endDate[6], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[7], 
                xmax = dates$endDate[7], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[8], 
                xmax = dates$endDate[8], 
                ymin = -40, ymax = -15), alpha = 0.2) 

threstrig <- ggplot() + #this is the ALERT dates
  theme_classic() +
  labs(y = "Influenza A cases", x = "Threshold-based intervention periods") +
  geom_bar(aes(y=chco$Cases, x=chco$Date), stat="identity") +
  geom_hline(aes(yintercept=threshold), linetype="dashed", 
             alpha=0.5, show.legend = FALSE) +
  geom_rect(aes(xmin = ALERT_dates$xstart[1], 
                xmax = ALERT_dates$xstop[1], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = ALERT_dates$xstart[2], 
                xmax = ALERT_dates$xstop[2], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = ALERT_dates$xstart[3], 
                xmax = ALERT_dates$xstop[3], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = ALERT_dates$xstart[4], 
                xmax = ALERT_dates$xstop[4], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = ALERT_dates$xstart[5], 
                xmax = ALERT_dates$xstop[5], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = ALERT_dates$xstart[6], 
                xmax = ALERT_dates$xstop[6], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = ALERT_dates$xstart[7], 
                xmax = ALERT_dates$xstop[7], 
                ymin = -40, ymax = -15), alpha = 0.2) +
  geom_rect(aes(xmin = ALERT_dates$xstart[8], 
                xmax = ALERT_dates$xstop[8], 
                ymin = -40, ymax = -15), alpha = 0.2) 

require(gridExtra)

grid.arrange(datetrig,threstrig)






















######################################
##### MODEL FORMULATION ##############
######################################

#model formulation

#real data
##load the data
fluA <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/CHCO-fluA.csv", 
                 stringsAsFactors=F)
fluA$Date <- ymd(fluA$Date)
fluA <- arrange(fluA, Date)

fluB <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/CHCO-fluB.csv", 
                 stringsAsFactors=F)
fluB$Date <- ymd(fluB$Date)
fluB <- arrange(fluB, Date)

RSV <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/CHCO-RSV.csv", 
                stringsAsFactors=F)
RSV$Date <- ymd(RSV$Date)
RSV <- arrange(RSV, Date)

Cases <- fluA$Cases + fluB$Cases + RSV$Cases
Date <- fluA$Date

data <- data.frame(Date, Cases)
fulldata <- data


##load the data
#chco <- read.csv("~/Desktop/applied-ALERT-data/CHCO-fluA.csv", 
#                 stringsAsFactors=F)

#chco$Date <- ymd(chco$Date)

chco <- arrange(fulldata, Date)

#train <- filter(chco, Date<ymd("2010-07-01"))

#test <- filter(chco, Date>ymd("2010-07-01"))

chco$state <- 0

fluDisProg <- create.disProg(week = 1:nrow(chco), 
                                    observed = chco$Cases,
                                    state = chco$state, start = c(2001, 35))
# convert to sts class
flu <- disProg2sts(fluDisProg)

#specify the model for endemic and autoregressive components
f.end <- addSeason2formula(f = ~ 1+t, S=1, period=52)
model1 <- list(ar = list(f = ~ 1), end = list(f =f.end),
               family = "NegBinM", subset=2:469)

# run model
res <- hhh4(flu, model1)

summaryres <- summary(res, idx2Exp=1, amplitudeShift=TRUE)

res$coefficients

res$se

#print(xtable(res$coefficients, type="html"))

#check the model

onesteppred <- oneStepAhead(res, nrow(chco)-1, type="rolling",
             which.start="current", verbose=FALSE)

#modelpredict <- data.frame(test, rev(onesteppred$pred), rev(onesteppred$observed))

colnames(modelpredict)[3:4] <- c("predicted", "observed")

grid.arrange(ggplot(modelpredict, aes(y=Cases, x=Date)) + 
  geom_bar(stat="identity", alpha=.5) + theme_classic() +
  geom_line(aes(x=Date, y=predicted)), 
ggplot(modelpredict, aes(x=Cases, y=predicted)) + 
  geom_point() + theme_classic() )

  

#try simulation
sim1 <- simulate(res, nsim=1, seed=NULL,
                y.start=NULL, simplify=TRUE)

sim2 <- simulate(res, nsim=1, seed=NULL,
                y.start=NULL, simplify=TRUE)

sim3 <- simulate(res, nsim=1, seed=NULL,
                y.start=NULL, simplify=TRUE)

fakeflu1 <- data.frame(chco$Date, as.vector(sim1))

fakeflu2 <- data.frame(chco$Date, as.vector(sim2))

fakeflu3 <- data.frame(chco$Date, as.vector(sim3))

grid.arrange(
  ggplot(fakeflu1, aes(x=chco.Date, y=as.vector.sim1.)) +
  geom_bar(stat="identity") +
  theme_classic(), 
  ggplot(fakeflu2, aes(x=chco.Date, y=as.vector.sim2.)) +
    geom_bar(stat="identity") +
    theme_classic(), 
  ggplot(fakeflu3, aes(x=chco.Date, y=as.vector.sim3.)) +
    geom_bar(stat="identity") +
    theme_classic()
  )

ggplot() +
  geom_bar(aes(x=fakeflu1$chco.Date[332:358], y=fakeflu1$as.vector.sim1.[332:358]), stat="identity") +
  geom_hline(aes(yintercept=4), linetype="dashed", 
             alpha=0.5, show.legend = FALSE) +
  geom_rect(aes(xmin = fakeflu1$chco.Date[340], 
                xmax = fakeflu1$chco.Date[349], 
                ymin = -.3, ymax = Inf), alpha = 0.2) +
  theme_classic() 


holder <- list()

holder[[1]] <- res$coefficients

holder[[2]] <- res$se

lower <- res$coefficients-res$se

upper <- res$coefficients+res$se

ar <- seq(from=lower[1], to=upper[1], length.out = 20)

beta <- seq(from=lower[3], to=upper[3], length.out = 20)

ar.beta <- expand.grid(ar, beta)

simulations <- list()

for (i in 1:nrow(ar.beta)){
  res$coefficients[1] <- ar.beta$Var1[i]
  res$coefficients[3] <- ar.beta$Var2[i]
  simulations[[i]] <- as.vector(simulate(res, nsim=1, seed=NULL,
           y.start=NULL, simplify=TRUE))
}

#get interesting ALERT parameters
i <- seq(1, length(simulations), 1)

#this is for when there is a firstMonth error. It will be filtered in the following steps.
filler <- (createALERT(season[[1]], firstMonth=1))

get_stats <- function (j) {
  print(j)
  result <- try(createALERT(season[[j]], firstMonth=1))
  if (class(result)=='try-error') {
    alert_stats <- filler$out[,1:3]
    alert_summaries <- (c(apply(alert_stats,
                                2, function (x) {max(x)*0}), parameters[j,]))
    success <- FALSE
  } else {
    alert_stats <- result$out[,1:3]
    alert_summaries <- (c(apply(alert_stats,
                                2, function (x) {median(x)}), parameters[j,]))
    success <- TRUE
  }
  alert_summaries['success'] <- success
  return(alert_summaries)
}

alertstats <- lapply(i, get_stats)

alertstatsposter <- list()


#for plotting

fakeflu1 <- data.frame(chco$Date, as.vector(simulations[[1]]))

fakeflu2 <- data.frame(chco$Date, as.vector(simulations[[2]]))

fakeflu3 <- data.frame(chco$Date, as.vector(simulations[[3]]))

grid.arrange(
  ggplot(fakeflu1, aes(x=chco.Date, y=as.vector.simulations..1...)) +
    geom_bar(stat="identity") +
    theme_classic(), 
  ggplot(fakeflu2, aes(x=chco.Date, y=as.vector.simulations..2...)) +
    geom_bar(stat="identity") +
    theme_classic(), 
  ggplot(fakeflu3, aes(x=chco.Date, y=as.vector.simulations..3...)) +
    geom_bar(stat="identity") +
    theme_classic()
)



















###get the other diseases

all <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/chco.csv", stringsAsFactors = F)


##alternatively, use the following.
#real data
##load the data
fluA <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/CHCO-fluA.csv", 
                 stringsAsFactors=F)
fluA$Date <- ymd(fluA$Date)
fluA <- arrange(fluA, Date)

fluB <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/CHCO-fluB.csv", 
                 stringsAsFactors=F)
fluB$Date <- ymd(fluB$Date)
fluB <- arrange(fluB, Date)

RSV <- read.csv("C:/Users/acbro0/Desktop/Projects/applied-ALERT-data/CHCO-RSV.csv", 
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

library(dplyr)
all <- arrange(all, by=Date)

all$state <- 0



holder <- list()

#specify the model for endemic and autoregressive components
f.end <- addSeason2formula(f = ~ 1+t, S=1, period=52)
model1 <- list(ar = list(f = ~ 1), end = list(f =f.end),
               family = "NegBinM", subset=2:469)

#RSVDisProg <- create.disProg(week = 1:nrow(all), 
#                            observed = all$Cases,
#                            state = all$state, start = c(2001, 35))

#for (i in 2:12){
#  print(i)
  
  RSVDisProg <- create.disProg(week = 1:nrow(all), 
                             observed = as.matrix(all[2]),
                             state = all$state, start = c(2001, 35))# convert to sts class
  rsv <- disProg2sts(RSVDisProg)# convert to sts class
#  # run model
  rsv_res <- hhh4(rsv, model1)
#  j <- as.integer(i-1)
  holder <- rsv_res$coefficients
#}

all <- data.frame(data.table::rbindlist(lapply(holder,as.list)))

colnames(all) <- "coefficients derived from CHCO dataset"

coefnames <- names(holder)
#rownames(all) <- bugnames



print(xtable(all, type="html"))



holder <- list()

for (i in 1:length(all)){
  holder[[i]] <- range(all[i])+(range(all[i])*0.1)
}

ranges <- t(data.frame(data.table::rbindlist(lapply(holder,as.list))))

colnames(ranges) <- coefnames

ranges <- data.frame(ranges)

holder <- list()

for (i in 1:length(ranges)){
  holder[[i]] <- seq(from=ranges[1,i], to=ranges[2,i], length.out=100)
}

names(holder) <- coefnames

coefs <- rsv_res$coefficients

simulations <- list()

#dates <- chco$Date



#this is for when there is a firstMonth error. It will be filtered in the following steps.
filler <- (createALERT(data, firstMonth=9))

get_stats <- function (j, params) {
  result <- try(createALERT(j, firstMonth=9))
  if (class(result)=='try-error') {
    alert_stats <- filler$out[,1:3]
    alert_summaries <- (c(apply(alert_stats,
                                2, function (x) {max(x)*0}), params))
    success <- FALSE
  } else {
    alert_stats <- result$out[,1:3]
    alert_summaries <- try((c(apply(alert_stats,
                                2, function (x) {median(x)}), params)))
    if (class(alert_summaries)=='try-error'){
      alert_stats <- filler$out[,1:3]
      alert_summaries <- (c(apply(alert_stats,
                                  2, function (x) {max(x)*0}), params))
      success <- FALSE
    } else {
    success <- TRUE
    }
  }
  alert_summaries['success'] <- success
  return(alert_summaries)
}


#choose the number of simulations to perform

snum <- 10

#create the empty data.frame to allocate memory.

num <- (length(holder)*length(holder[[1]])*snum)


alertstats2 <- data.frame(threshold=rep(NA, num),
                          median.dur=rep(NA, num),
                          median.pct.cases.captured=rep(NA, num),
                          parameter=rep(NA, num),
                          value=rep(NA, num),
                          success=rep(NA, num))

#fill the data.frame using the power of indexing

  for (j in 1:length(holder)){
    print(c(j))
    params <- c(names(holder)[i], holder[[i]][[j]])
    toycoef <- coefs
    toycoef[i] <- holder[[i]][[j]]
    res2 <- hhh4(rsv, model1)
    res2$coefficients <- toycoef
    simulation <- as.vector(simulate(res2, nsim=snum, seed=NULL,
                                           y.start=NULL, simplify=TRUE))
    minisim <- split(simulation, as.factor(rep(seq(1,snum,1), each=length(dates))))
    for(k in 1:snum){
      simset <- data.frame(minisim[[k]], dates)
      colnames(simset) <- c("Cases", "Date")
      alertstats2[j+(i-1)*length(holder[[i]])+
                    (length(holder)*length(holder[[i]])*(k-1)),] <- get_stats(simset, 
                                                           params=params) 
    }
  }
}

med.stats <- alertstats2 %>% mutate(median.pct.cases.captured=as.numeric(median.pct.cases.captured), 
                       threshold=as.numeric(threshold), 
                       median.dur=as.numeric(median.dur)) %>% 
  filter(success==TRUE) %>% group_by(parameter, value) %>%
  summarise(median=median(median.pct.cases.captured), 
            threshold=median(threshold),
            median.dur=median(median.dur)) %>%
  data.frame()

#it looks like ar.1, end.1, and end.t are going to be the most interesting.
#maybe overdispersion?

#plot this bit with each parameter adjusted while holding the others constant.

#I'm filtering the last row (overdispersion) because the threshold is crazy high
#and medians are either NA or 0.

med.stats <- med.stats[-60,]

ggplot(med.stats, group=parameter) +
  geom_point(aes(x=as.numeric(value), y=median, size=median.dur, color=threshold)) +
  scale_colour_gradient2(low = "blue", mid = "green", high = "yellow", 
                         midpoint = 1000, na.value = "grey50", guide = "colourbar") +
     facet_wrap(~parameter, ncol=3)+
  theme_classic() 

#I need some examples of datasets created with the sin and cos params
#as I don't exactly understand what they are doing

coefs2 <- coefs

coefs2[1] <- -10000

coefs2[5] <- 0

coefs2[6] <- 0

coefs2[2] <- 0

coefs2[3] <- 0

toycoef <- coefs2
toycoef[4] <- holder[[4]][[1]]
res2 <- hhh4(flu, model1)
res2$coefficients <- coefs2
simulation <- as.vector(simulate(res2, nsim=1, seed=NULL,
                                 y.start=NULL, simplify=TRUE))
simset <- data.frame(simulation, dates)
colnames(simset) <- c("Cases", "Date")

ggplot(simset, aes(x=Date, y=Cases)) +
  geom_line(stat="identity") +
  theme_classic()

coefs2[4] <- 0

toycoef <- coefs2
toycoef[5] <- holder[[5]][[10]]
res2 <- hhh4(flu, model1)
res2$coefficients <- toycoef
simulation <- as.vector(simulate(res2, nsim=1, seed=NULL,
                                 y.start=NULL, simplify=TRUE))
simset <- data.frame(simulation, dates)
colnames(simset) <- c("Cases", "Date")

ggplot(simset, aes(x=Date, y=Cases)) +
  geom_line(stat="identity") + 
  theme_classic()

#The sin term manages the width of the seasonal epidemic/peak height/baseline noise

#cos is the width of the sinus, Larger cos term should make the sinus more narrow. 
#smaller makes the sinus wider (less baseline noise)

#head(expand.grid(holder))

cor(all)

###################################
#use ALERT on each CHCO dataset###
###################################

head(all)

all$Date <- mdy(all$Date..Month.Year.) #FIXME

#RSV

rsv <- data.frame(all$Date, all$RSV)

colnames(rsv) <- c("Date", "Cases")

rsvALERT <- createALERT(rsv)[1]

#flu b
flub <- data.frame(all$Date, all$Flu.B)

colnames(flub) <- c("Date", "Cases")

flubALERT <- createALERT(flub)[1]

#HMPV
hmpv <- na.omit(data.frame(all$Date, all$HMPV))

colnames(hmpv) <- c("Date", "Cases")

hmpvALERT <- createALERT(hmpv)[1]

#Paraflu
pflu <- data.frame(all$Date, all$Paraflu)

colnames(pflu) <- c("Date", "Cases")

pfluALERT <- createALERT(pflu)[1]

#Adenovirus
advir <- data.frame(all$Date, all$Adenovirus)

colnames(advir) <- c("Date", "Cases")

advirALERT <- createALERT(advir)[1]

#Rhino

rhino <- data.frame(all$Date, all$Rhinovirus)

colnames(rhino) <- c("Date", "Cases")

rhinoALERT <- createALERT(rhino)[1]

#Coronavirus
corvir <- na.omit(data.frame(all$Date, all$Coronavirus))

colnames(corvir) <- c("Date", "Cases")

corvirALERT <- createALERT(corvir)[1]

#Pertussis
whoop <- na.omit(data.frame(all$Date, all$B..Pertussis))

colnames(whoop) <- c("Date", "Cases")

whoopALERT <- createALERT(whoop)[1]

#Enteroviruses
entero <- na.omit(data.frame(all$Date, all$Enterovirus))

colnames(entero) <- c("Date", "Cases")

enteroALERT <- createALERT(entero)[1]

#Diarrhea.Viruses
diavir <- na.omit(data.frame(all$Date, all$Diarrhea.Viruses))

colnames(diavir) <- c("Date", "Cases")

diavirALERT <- createALERT(diavir)[1]

grid.arrange(
  ggplot(rsv, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Respiratory Syncycial Virus"), 
  ggplot(flub, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Influenza B"),   
  ggplot(hmpv, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Human Metapneumovirus"), 
  ggplot(pflu, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Parainfluenza"),
  ggplot(advir, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Adenovirus"),
  ggplot(rhino, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Rhinovirus"), 
  ggplot(corvir, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Coronavirus"),
  ggplot(whoop, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Pertussis"),
  ggplot(entero, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Enterovirus"), 
  ggplot(diavir, aes(x=Date, y=Cases)) +
    geom_bar(stat="identity") +
    theme_classic()+
    ggtitle("Diarrhea Viruses")
)

#graph parameters against each other to look for interesting clusters.

# CUSUM

disProgObj <- create.disProg(week=chco_train$Date, 
                             observed= chco_train$Cases,
                             state=rep(1,210))

res <- algo.cusum(disProgObj, 
                  control = list(range = 1:210, 
                                 k=36, h = .2, 
                                 m=NULL,
                                 trans="none"))

chco_train$alarm <- res$alarm
chco_train$state <- c(rep(1, 40), rep (2, 50), rep (3, 50), rep (4, 70))

tab <- chco_train %>% group_by(alarm, state) %>% summarize(sum=sum(alarm))
tab[5:8,]$sum
median(tab[5:8,]$sum)

disProgObj <- create.disProg(week=chco$Date, 
                             observed= chco$Cases,
                             state=rep(1,414))

res <- algo.cusum(disProgObj, 
                  control = list(range = 1:414, 
                                 k=36, h = .2, 
                                 m=NULL,
                                 trans="none"))
chco$alarm <- res$alarm

chco$alertperiod <- ifelse(chco$alarm==1, TRUE, FALSE)


#assign the year
A <- rep(1, 42)
for(i in 2:7){
  A <- append(rep(i, 51), A)
}
A <- append(rep(8, nrow(chco)-length(A)), A)
chco$year <- as.factor(rev(A))

nrow(chco_train)
nrow(chco_test)

#chco1 <- chco
chco <- chco1
chco <- chco[196:nrow(chco),]
#chco <- chco[1:195,]

tot <- chco %>% group_by(year) %>% 
  summarise(total_cases=sum(Cases)) #total cases per year

alertcases <- chco %>% filter(chco$alertperiod==TRUE) %>% 
  group_by(year) %>% 
  summarise(alert_cases=sum(Cases)) #"ALERT" cases per year

calc <- full_join(tot, alertcases)

calc[is.na(calc)] <- 0

calc$perc_alert <- calc$alert_cases/calc$total_cases *100 #percent alert cases captured

median(calc$perc_alert) #median cases captured

alert_dura <- chco %>% filter(alertperiod==TRUE) %>% 
  group_by(year) %>% #duration
  summarise(duration=as.numeric(length(Date)))

calc <- unique(full_join(calc, alert_dura))

holder <- chco %>% group_by(year) %>% 
  summarise(Cases=max(Cases)) %>%
  data.frame()

holder <- left_join(holder, chco)

holder <- holder %>% group_by(year) %>% summarise(Date=first(Date))

holder <- left_join(holder, chco) %>% data.frame()

calc$peak_captured <- holder$alertperiod

peaky <- table(calc$peak_captured) %>% data.frame() 

if(nrow(peaky)==2){
  perc_peaks_captured <- peaky$Freq[2] / peaky$Freq[1]*100
} else {
  if (peaky$Var1==TRUE){
    perc_peaks_captured <- 100
  } else {
    perc_peaks_captured <- 0
  }
}

holder <- chco %>% group_by(year) %>% summarise(Cases=min(Cases)) %>%
  data.frame()

holder <- left_join(holder, chco)

#number of 0 weeks included
holder <- filter(holder, alertperiod==TRUE) %>% group_by(year) %>% 
  summarise(low_weeks_incl=length(alertperiod))

calc <- left_join(calc, holder) %>% data.frame()

calc$year <- NULL

calc[is.na(calc)] <- 0

cutoffs <- dates %>% arrange(-yearIdx)

calc <- cbind(cutoffs, calc)

calc$yearIdx <- NULL

#get the colnames for eval statistics

alertholder <- data.frame(alertholder$out)

snames <- colnames(alertholder)[-7]

stats_CUSUM <- t(data.frame(c(NA, 
                             round(median(calc$duration), 1),
                             round(median(calc$perc_alert), 1), 
                             round(min(calc$perc_alert), 1), 
                             round(max(calc$perc_alert), 1), 
                             round(perc_peaks_captured, 3), 
                             round(mean(calc$low_weeks_incl), 1))))

colnames(stats_CUSUM) <- snames

require(xtable)

print(xtable(stats_CUSUM), include.rownames=F)

str(stats_real)

