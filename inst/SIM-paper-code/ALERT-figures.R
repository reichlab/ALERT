#### figures for the applied ALERT paper
require(ALERT)
require(surveillance)
require(lubridate); year <- lubridate::year
require(ggplot2)
library(dplyr)
library(gridExtra)

source("inst/SIM-paper-code/ALERT-simulations-fxns.R")

#real data
##load the data
fluA <- read.csv("inst/SIM-paper-data/CHCO-fluA.csv", 
                 stringsAsFactors=F)
fluA$Date <- ymd(fluA$Date)
fluA <- arrange(fluA, Date)

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
fulldata <- data




###get the dates that CHCO actually used
dates <- read.csv("inst/SIM-paper-data/chco-trigger-dates.csv", stringsAsFactors = F)
dates$startDate <- ymd(dates$startDate)
dates$endDate <- ymd(dates$endDate)

#truncate to get years we have cutoffs for and cutoffs we have years for.
data <- filter(fulldata, Date >= head(dates$startDate, 1)-weeks(2))
#testingdata <- filter(fulldata, Date < head(dates$startDate, 1)-weeks(2))
dates <- filter(dates, endDate <= tail(data$Date, 1)+weeks(6))


############################################
#### ALERT DEMO FIGURE   ###################
############################################

rawdat <- ggplot(data) +
  geom_ribbon(stat="identity", aes(x=Date, ymax=Cases), ymin=0) +
 # scale_x_date(limits=c(as.Date("2001-07-07"), as.Date("2005-01-01")))+
  theme_classic()+
# ggtitle("Training years") +
  geom_hline(aes(yintercept=25), linetype="dashed", 
             alpha=0.5, show.legend = FALSE) +
  theme(axis.title.y=element_blank())
  
##zeros are dropped because that's what ALERT does
dens <- ggplot(data=filter(data, Cases>0), 
       aes(filter(data, Cases>0)$Cases)) +
  theme_classic()+
  geom_density(fill="grey") +
  coord_flip() +
  geom_vline(aes(xintercept=25), linetype="dashed", 
             alpha=0.5, show.legend = FALSE) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  labs(y = "Density") 

pdf(file="inst/SIM-paper-code/fig1-historical-data.pdf", width=7, height=4)
grid.arrange(rawdat, dens, ncol=2, left="Incidence")
dev.off()

#######################################################
###### ALERT THRESHOLD COMPARISON FIGURE  #############
#######################################################

##Run ALERT

alerttrain <- createALERT(data[1:193,], firstMonth=8, lag=0) #train on 2004-2008
true <- datetestALERT(data[1:193,], dates[1:4,1], dates[1:4,2])

trainstats <- full_join(data.frame(alerttrain[[1]]), data.frame(true[[1]]))

test <- thresholdtestALERT(data, firstMonth = 8, whichThreshold = 6)
ALERT_dates3 <- data.frame(test$details)
##convert dates to something human readable for plotting
ALERT_dates3$start <- as.Date(ALERT_dates3$start, origin="1970-01-01")
ALERT_dates3$end <- as.Date(ALERT_dates3$end, origin="1970-01-01")

test <- thresholdtestALERT(data, firstMonth = 8, whichThreshold = 10)
ALERT_dates12 <- data.frame(test$details)
##convert dates to something human readable for plotting
ALERT_dates12$start <- as.Date(ALERT_dates12$start, origin="1970-01-01")
ALERT_dates12$end <- as.Date(ALERT_dates12$end, origin="1970-01-01")

test <- thresholdtestALERT(data, firstMonth = 8, whichThreshold = 21)
ALERT_dates22 <- data.frame(test$details)
##convert dates to something human readable for plotting
ALERT_dates22$start <- as.Date(ALERT_dates22$start, origin="1970-01-01")
ALERT_dates22$end <- as.Date(ALERT_dates22$end, origin="1970-01-01")


###get the summary stats for the table

alerttest <- createALERT(data[194:nrow(data),], firstMonth=8, 
                          lag=0, allThresholds = TRUE) 

holder <- data.frame(alerttest[[1]]) %>% 
  filter(threshold==1 |
           threshold==2 |
           threshold==4 |
           threshold==6 |
           threshold==10 |
           threshold==21)


teststats <- full_join(holder, 
          data.frame(datetestALERT(data[194:nrow(data),], 
                                   dates[5:8,1], dates[5:8,2])[[1]]))

comparestats <- full_join(trainstats, teststats)

require(xtable)

print(xtable(comparestats), include.rownames=F)

datetestALERT(data, dates[,1], dates[,2])  ##full real data 
data.frame(createALERT(data, firstMonth=8, 
            lag=0, allThresholds = TRUE) [[1]]) %>%
  filter(threshold==6 | threshold==10 | threshold==21)

# get the CUSUM dates 

disProgObj <- create.disProg(week=data$Date, 
                             observed= data$Cases,
                             state=rep(1,414))

res <- algo.cusum(disProgObj, 
                  control = list(range = 1:414, 
                                 k=36, h = .2, 
                                 m=NULL,
                                 trans="none"))
data$alarm <- res$alarm

startDate1 <- data$Date[ifelse(data$alarm==1 & lag(data$alarm==0), TRUE, FALSE)==TRUE]
endDate1 <- data$Date[ifelse(data$alarm==0 & lag(data$alarm==1), TRUE, FALSE)==TRUE]

CUSUM_dates <- data.frame(startDate1, endDate1[2:length(endDate1)])

colnames(CUSUM_dates) <- c("startDate", "endDate")

ymin <- -2

threstrig <- ggplot() + #this is the ALERT dates
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14) ) +
 # ylim(-40, 150) +
  labs(y = "Incidence", x = "Date") +
  geom_ribbon(aes(ymax=data$Cases, x=data$Date),
              ymin=0,
           stat="identity") +
  geom_rect(aes(xmin = dates$startDate[1], 
                xmax = dates$endDate[1], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[2], 
                xmax = dates$endDate[2], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[3], 
                xmax = dates$endDate[3], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[4], 
                xmax = dates$endDate[4], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[5], 
                xmax = dates$endDate[5], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[6], 
                xmax = dates$endDate[6], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[7], 
                xmax = dates$endDate[7], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  geom_rect(aes(xmin = dates$startDate[8], 
                xmax = dates$endDate[8], 
                ymin = -27, ymax = 0), alpha = 0.2) +
  
  geom_rect(aes(xmin = ALERT_dates3$start[1], 
                xmax = ALERT_dates3$end[1], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates3$start[2], 
                xmax = ALERT_dates3$end[2], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates3$start[3], 
                xmax = ALERT_dates3$end[3], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates3$start[4], 
                xmax = ALERT_dates3$end[4], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates3$start[5], 
                xmax = ALERT_dates3$end[5], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates3$start[6], 
                xmax = ALERT_dates3$end[6], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates3$start[7], 
                xmax = ALERT_dates3$end[7], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates3$start[8], 
                xmax = ALERT_dates3$end[8], 
                ymin = -5, ymax = -10), alpha = 0.5) +
  
  geom_rect(aes(xmin = ALERT_dates12$start[1], 
                xmax = ALERT_dates12$end[1], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates12$start[2], 
                xmax = ALERT_dates12$end[2], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates12$start[3], 
                xmax = ALERT_dates12$end[3], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates12$start[4], 
                xmax = ALERT_dates12$end[4], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates12$start[5], 
                xmax = ALERT_dates12$end[5], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates12$start[6], 
                xmax = ALERT_dates12$end[6], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates12$start[7], 
                xmax = ALERT_dates12$end[7], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates12$start[8], 
                xmax = ALERT_dates12$end[8], 
                ymin = -11.5, ymax = -16), alpha = 0.5) +
  
  geom_rect(aes(xmin = ALERT_dates22$start[1], 
                xmax = ALERT_dates22$end[1], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates22$start[2], 
                xmax = ALERT_dates22$end[2], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates22$start[3], 
                xmax = ALERT_dates22$end[3], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates22$start[4], 
                xmax = ALERT_dates22$end[4], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates22$start[5], 
                xmax = ALERT_dates22$end[5], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates22$start[6], 
                xmax = ALERT_dates22$end[6], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates22$start[7], 
                xmax = ALERT_dates22$end[7], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  geom_rect(aes(xmin = ALERT_dates22$start[8], 
                xmax = ALERT_dates22$end[8], 
                ymin = -18, ymax = -23), alpha = 0.5) +
  #CUSUM
  geom_rect(aes(xmin = CUSUM_dates$start[1], 
                xmax = CUSUM_dates$end[1], 
                ymin = -27, ymax = -32), alpha = 0.5) +
  geom_rect(aes(xmin = CUSUM_dates$start[2], 
                xmax = CUSUM_dates$end[2], 
                ymin = -27, ymax = -32), alpha = 0.5) +
  geom_rect(aes(xmin = CUSUM_dates$start[3], 
                xmax = CUSUM_dates$end[3], 
                ymin = -27, ymax = -32), alpha = 0.5) +
  geom_rect(aes(xmin = CUSUM_dates$start[4], 
                xmax = CUSUM_dates$end[4], 
                ymin = -27, ymax = -32), alpha = 0.5) +
  geom_rect(aes(xmin = CUSUM_dates$start[5], 
                xmax = CUSUM_dates$end[5], 
                ymin = -27, ymax = -32), alpha = 0.5) +
  geom_rect(aes(xmin = CUSUM_dates$start[6], 
                xmax = CUSUM_dates$end[6], 
                ymin = -27, ymax = -32), alpha = 0.5) +
  geom_rect(aes(xmin = CUSUM_dates$start[7], 
                xmax = CUSUM_dates$end[7], 
                ymin = -27, ymax = -32), alpha = 0.5) +
  geom_rect(aes(xmin = CUSUM_dates$start[8], 
                xmax = CUSUM_dates$end[8], 
                ymin = -27, ymax = -32), alpha = 0.5) +
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.title.y=element_text()) +
  annotate("text", x = as.Date("2004-10-05"), y = -7, label = "6", size=3)+
  annotate("text", x = as.Date("2004-10-05"), y = -13.815, label = "10", size=3) +
annotate("text", x = as.Date("2004-10-05"), y = -21, label = "21", size=3) +
  annotate("text", x = as.Date("2004-09-10"), y = -29, label = "CUSUM", size=3) +
  geom_vline(xintercept=as.Date("2008-07-26"), linetype="dashed", 
                 alpha=0.5, show.legend = FALSE )

pdf(file="inst/SIM-paper-code/fig2-alert-periods.pdf", width=7, height=4)
threstrig
dev.off()



