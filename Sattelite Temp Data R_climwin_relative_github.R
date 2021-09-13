rm(list = ls())

require(readr)
require(lubridate)
require(plyr)
require(ggplot2)
require(EnvStats)


#upload remotely sensed temperature data
Beach_Haven <- read_csv("Beach Haven.csv")
False_Bay <- read_csv("False Bay.csv")
Fisherman_Bay <- read_csv("Fisherman Bay.csv")
Indian_Cove <- read_csv("Indian Cove.csv")
Mosquito_Pass <- read_csv("Mosquito Pass.csv")
North_Cove <- read_csv("North Cove.csv")
Picnic_Cove2 <- read_csv("Picnic Cove2.csv")
Ship_Harbor <- read_csv("Ship Harbor.csv")
Shoal_Bay <- read_csv("Shoal Bay.csv")

#merge datasets and organize so that it is compatible with climwin package
Sat_Temps<-rbind(Beach_Haven, False_Bay, Fisherman_Bay, Indian_Cove, Mosquito_Pass, North_Cove, Picnic_Cove2, Ship_Harbor, Shoal_Bay)
Sat_Temps$SiteID<-as.factor(ifelse(Sat_Temps$site=='Beach Haven', 'A', 
                         ifelse(Sat_Temps$site=='False Bay', 'B', 
                                ifelse(Sat_Temps=='Fisherman Bay','C', 
                                       ifelse(Sat_Temps=='Indian Cove', 'D', 
                                              ifelse(Sat_Temps=='Mosquito Pass', 'E', 
                                                     ifelse(Sat_Temps=='North Cove','F', 
                                                            ifelse(Sat_Temps=='Picnic Cove', 'G', 
                                                                   ifelse(Sat_Temps=='Ship Harbor', 'H', 
                                                                          ifelse(Sat_Temps=='Shoal Bay', 'I', 'Z'))))))))))
Sat_Temps$date<-as.Date(Sat_Temps$date,"%m/%d/%y")
Sat_Temps$year<-year(Sat_Temps$date)
Sat_Temps$date<-factor(format(Sat_Temps$date, "%d/%m/%Y"))
Sat_Temps<-Sat_Temps[Sat_Temps$SST<172,] #remove error codes
Sat_Temps<-Sat_Temps[!is.na(Sat_Temps$SST),]


#import eelgrass data
EGWD<- read_csv("EGWD 2013-2017_MEPS.csv")

#organize data for 'climwin' analysis
EGWD$date<-as.Date(EGWD$date, "%d/%m/%Y")
EGWD$julian<-yday(EGWD$date)
EGWD$date<-factor(format(EGWD$date, "%d/%m/%Y"))

EGWD$SiteID<-as.factor(ifelse(EGWD$site=='Beach Haven', 'A', 
                         ifelse(EGWD$site=='False Bay', 'B', 
                                ifelse(EGWD$site=='Fisherman Bay','C', 
                                       ifelse(EGWD$site=='Indian Cove', 'D', 
                                              ifelse(EGWD$site=='Mosquito Pass', 'E', 
                                                     ifelse(EGWD$site=='North Cove','F', 
                                                            ifelse(EGWD$site=='Picnic Cove', 'G', 
                                                                   ifelse(EGWD$site=='Ship Harbor', 'H', 
                                                                          ifelse(EGWD$site=='Shoal Bay', 'I', 'Z'))))))))))

EGWD<-EGWD[!is.na(EGWD$leafArea),]
EGWD$leafArea_dm2<-EGWD$leafArea/1000
EGWD$year<-as.factor(EGWD$year)
EGWD<-EGWD[EGWD$depth!='Middle',]
EGWD<-EGWD[EGWD$leafArea_dm2<20,]#remove a few extreme outliers


require(climwin)
require(AICcmodavg)
###First analyze EWD presence/ absence
#Run preliminary models to determine which predictors (aside from temperature) should be included, then use climwin to test temperature windows

model.a<-glmer(diseased~(1|site), data=EGWD, family="binomial")
model.b<-glmer(diseased~year+(1|site), data=EGWD, family="binomial")
model.c<-glmer(diseased~leafArea_dm2+(1|site), data=EGWD, family="binomial")
model.d<-glmer(diseased~depth+(1|site), data=EGWD, family="binomial")
model.e<-glmer(diseased~year+leafArea_dm2+(1|site), data=EGWD, family="binomial")
model.f<-glmer(diseased~depth+leafArea_dm2+(1|site), data=EGWD, family="binomial")
model.g<-glmer(diseased~year+depth+(1|site), data=EGWD, family="binomial")
model.h<-glmer(diseased~year+depth+leafArea_dm2+(1|site), data=EGWD, family="binomial")
AICc(model.a)
AICc(model.b)
AICc(model.c)
AICc(model.d)
AICc(model.e)
AICc(model.f)
AICc(model.g)
AICc(model.h)


#now include best model in climwin analysis
model<-glmer(diseased~1+year+leafArea_dm2+depth+(1|site), data=EGWD, family="binomial")

MassWin_Rela <- slidingwin(xvar = list(Temp = Sat_Temps$SST),
                      cdate = Sat_Temps$date,
                      bdate = EGWD$date,
                      baseline = model,
                      cohort = as.factor(EGWD$year),
                      cinterval = "month",
                      range = c(6, 0),
                      type = "relative",
                      cmissing='method2',
                      stat = "mean",
                      func = "lin", spatial = list(EGWD$SiteID, Sat_Temps$SiteID))


summary(MassWin_Rela[[1]]$BestModel)

#run 100 randomized models to compare result against
MassRand <- randwin(repeats = 100, xvar = list(Temp = Sat_Temps$SST),
                    cdate = Sat_Temps$date,
                    bdate = EGWD$date,
                    baseline = model,
                    cohort = as.factor(EGWD$year),
                    cinterval = "month",
                    range = c(6, 0),
                    type = "relative",
                    cmissing='method2',
                    stat = "mean",
                    func = "lin", spatial = list(EGWD$SiteID, Sat_Temps$SiteID))





head(MassWin_Rela[[1]]$Dataset)

MassWin_Rela[[1]]$BestModel



check<-(MassWin_Rela[[1]]$BestModelData)

pvalue(dataset = MassWin_Rela[[1]]$Dataset, datasetrand = MassRand[[1]], metric = "C", sample.size = 45)

MassOutput<-MassWin_Rela[[1]]$Dataset

MassOutput$WindowOpen<-MassOutput$WindowOpen+1

plotdelta(dataset = MassOutput)

p1<-ggplot(MassOutput, aes(x=WindowClose, y=WindowOpen, fill=deltaAICc))+geom_tile()+theme_bw()+ylab('Window Open (month)')+xlab('Window Close (month)')+
  scale_fill_continuous(low = "red", high = "grey")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.tag.position = c(0.1, 0.95),legend.position = c(0.8, 0.2))+
  labs(fill="Delta AICc", tag="a")

p1




#to calculate rsquared run best fit model alone
require(rsq)
mod<-glmer(yvar~climate+year+leafArea_dm2+depth+(1|site), data=check, family="binomial")

rsq(mod)

AICc(mod)

#and plot predicted results
newdata <- with(check, expand.grid(year=as.factor(seq(2013, 2017, 1)), site="Beach Haven", depth="Shallow", leafArea_dm2=seq(10,20, .1), climate=seq(11.5, 15.5, .1)))

fit<-predict(mod, type ="response", newdata,  se=TRUE)

newdata<-cbind(newdata, fit)
newdata<-newdata[newdata$depth=='Shallow',]


p<-ggplot(newdata, aes(x=climate, y=leafArea_dm2, fill=fit))+geom_tile()+scale_fill_gradient2(low = "blue", mid="yellow", high = "red", midpoint=0.5)
p+facet_grid(.~year)+theme_bw()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab(bquote('Leaf Area (' ~dm^2* ')'))+xlab(expression(paste("Average temperature (", degree,"C) for a 1 month window, beginning 1 month prior to sampling", )))+
  theme(panel.spacing.x = unit(4, "mm"))+labs(fill='EWD Prevalence')


####Analyze the Severity data
Severity_data<-EGWD[EGWD$severity>0,]#limit dataset to diseased leaves
Severity_data<-Severity_data[!is.na(Severity_data$severity),]
Severity_data<-Severity_data[Severity_data$year!='2013',]#remove years when severity was not sampled
Severity_data<-Severity_data[Severity_data$year!='2014',]
Severity_data<-Severity_data[!is.na(Severity_data$year),]

#Run modelswithout climate to determine which predictors to include
model.a<-lmer(severity~(1|site), data=Severity_data)
model.b<-lmer(severity~year+(1|site), data=Severity_data)
model.c<-lmer(severity~leafArea_dm2+(1|site), data=Severity_data)
model.d<-lmer(severity~depth+(1|site), data=Severity_data)
model.e<-lmer(severity~year+leafArea_dm2+(1|site), data=Severity_data)
model.f<-lmer(severity~depth+leafArea_dm2+(1|site), data=Severity_data)
model.g<-lmer(severity~depth+year+(1|site), data=Severity_data)
model.h<-lmer(severity~depth+year+leafArea_dm2+(1|site), data=Severity_data)
AICc(model.a)
AICc(model.b)
AICc(model.c)
AICc(model.d)
AICc(model.e)
AICc(model.f)
AICc(model.g)
AICc(model.h)

#models e has the lowest AICc. 

#relative time
require(lmerTest)
detach(package:lmerTest,unload=TRUE)#climwin models won't run if lmerTest is active

model<-lmer(severity~1+year+leafArea_dm2+(1|site), data=Severity_data)

MassWin_Rela <- slidingwin(xvar = list(Temp = Sat_Temps$SST),
                           cdate = Sat_Temps$date,
                           bdate = Severity_data$date,
                           baseline = model,
                           cohort = as.factor(Severity_data$year),
                           cinterval = "month",
                           range = c(6, 0),
                           type = "relative",
                           cmissing='method2',
                           stat = "mean",
                           func = "lin", spatial = list(Severity_data$SiteID, Sat_Temps$SiteID))

summary(MassWin_Rela[[1]]$BestModel)


MassRand <- randwin(repeats = 100, xvar = list(Temp = Sat_Temps$SST),
                    cdate = Sat_Temps$date,
                    bdate = Severity_data$date,
                    baseline = lm(severity~1+year+leafArea_dm2+site, data=Severity_data),
                    cohort = as.factor(Severity_data$year),
                    cinterval = "month",
                    range = c(6, 0),
                    type = "relative",
                    cmissing='method2',
                    stat = "mean",
                    func = "lin", spatial = list(Severity_data$SiteID, Sat_Temps$SiteID))

head(MassWin_Rela[[1]]$Dataset)

MassWin_Rela[[1]]$BestModel


pvalue(dataset = MassWin_Rela[[1]]$Dataset, datasetrand = MassRand[[1]], metric = "C", sample.size = 36) #p=0.0004

MassOutput<-MassWin_Rela[[1]]$Dataset
MassOutput$WindowOpen<-MassOutput$WindowOpen+1

p2<-ggplot(MassOutput, aes(x=WindowClose, y=WindowOpen, fill=deltaAICc))+geom_tile()+theme_bw()+ylab('Window Open (month)')+xlab('Window Close (month)')+
  scale_fill_continuous(low = "red", high = "grey")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.tag.position = c(0.1, 0.95),legend.position = c(0.8, 0.2))+
  labs(fill="Delta AICc", tag="b")

p2
plotdelta(dataset = MassOutput)


####plot predictions
check1<-(MassWin_Rela[[1]]$BestModelData)

mod1<-lm(yvar~climate+year+leafArea_dm2+site, data=check1)
require(rsq)
rsq(mod1)
AICc(mod1)
rsq(mod1)

newdata <- with(check1, expand.grid(year=as.factor(seq(2015, 2017, 1)), site="Beach Haven",  leafArea_dm2=seq(10,20, .1), climate=seq(9, 14, .1)))

fit<-predict(mod1, newdata, type="response", se=TRUE)

newdata<-cbind(newdata, fit$fit)
newdata<-cbind(newdata, fit$se.fit)



p<-ggplot(newdata, aes(x=climate, y=leafArea_dm2, fill=fit$fit))+geom_tile()
p+facet_grid(.~year)+theme_bw()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Predicted EWD Severity")+xlab(expression(paste("Average temperature (", degree,"C) for a 1 month period, beginning 3 months prior to sampling")))+
  theme(panel.spacing.x = unit(4, "mm"))+ scale_fill_gradient2(low = "blue", mid="yellow", high = "red", midpoint=0.04)+labs(fill='EWD severity')+ylab(bquote('Leaf Area'~(dm^2)))


#########analyze density data

#first make a data sset with summarized density data
SummaryEGWD<-ddply(EGWD, .(date, SiteID, year, depth), summarize, prevalence=sum(diseased)/length(diseased), MeanArea=mean(leafArea)/10000)

myvars <- c("date", "year", "site", "SiteID", "depth", "density", "leafArea")
EGWD_density<-EGWD[myvars]

EGWD_density<-EGWD_density[is.na(EGWD_density$density)==FALSE,]
Summarydensity<-ddply(EGWD_density, .(site, SiteID, year, depth), summarize, density=mean(density))
names(Summarydensity)[names(Summarydensity) == "Year"] <- "year"



EGWD<-EGWD[is.na(EGWD$diseased)==FALSE,]
EGWD<-EGWD[EGWD$julian>200,]
EGWD<-EGWD[EGWD$depth!="Middle",]#remove middle transect (only taken in a subset of surveys)
EGWD$severity<-ifelse(EGWD$severity>0, EGWD$severity, NA)

SummaryEGWD_density<-merge(SummaryEGWD, Summarydensity, by=c("year", "SiteID", "depth"))



#Run models without climate data to determine which predictors to include in climwin analysis
model.a<-lmer(density~(1|site), data=SummaryEGWD_density)
model.b<-lmer(density~year+(1|site), data=SummaryEGWD_density)
model.c<-lmer(density~depth+(1|site), data=SummaryEGWD_density)
model.d<-lmer(density~MeanArea+(1|site), data=SummaryEGWD_density)
model.e<-lmer(density~year+depth+(1|site), data=SummaryEGWD_density)
model.f<-lmer(density~year+MeanArea+(1|site), data=SummaryEGWD_density)
model.g<-lmer(density~depth+MeanArea+(1|site), data=SummaryEGWD_density)
model.h<-lmer(density~depth+year+MeanArea+(1|site), data=SummaryEGWD_density)
AICc(model.a)
AICc(model.b)
AICc(model.c)
AICc(model.d)
AICc(model.e)
AICc(model.f)
AICc(model.g)
AICc(model.h)# best fit model based on AICc, however depth and Leaf Area effect sizes have huge confidence intervals and very high p-values, unclear that they improve model fit


model=lmer(density~1+year+(1|site), data=SummaryEGWD_density)

MassWin_Rela <- slidingwin(xvar = list(Temp = Sat_Temps$SST),
                           cdate = Sat_Temps$date,
                           bdate = SummaryEGWD_density$date,
                           baseline = model, 
                           cohort = as.factor(SummaryEGWD_density$year),
                           cinterval = "month",
                           range = c(6, 0),
                           type = "relative",
                           cmissing='method2',
                           stat = "mean",
                           func = "lin", spatial = list(SummaryEGWD_density$SiteID, Sat_Temps$SiteID))


summary(MassWin_Rela[[1]]$BestModel)


MassRand <- randwin(repeats = 100,xvar = list(Temp = Sat_Temps$SST),
                    cdate = Sat_Temps$date,
                    bdate = SummaryEGWD_density$date,
                    baseline = model, 
                    cohort = as.factor(SummaryEGWD_density$year),
                    cinterval = "month",
                    range = c(6, 0),
                    type = "relative",
                    cmissing='method2',
                    stat = "mean",
                    func = "lin", spatial = list(SummaryEGWD_density$SiteID, Sat_Temps$SiteID))

head(MassWin_Rela[[1]]$Dataset)

summary(MassWin_Rela[[1]]$BestModel)

MassWin_Rela[[1]]$BestModel

pvalue(dataset = MassWin_Rela[[1]]$Dataset, datasetrand = MassRand[[1]], metric = "C", sample.size = 47) #p= 0.96, no evidence that the model is not by chance

MassOutput<-MassWin_Rela[[1]]$Dataset
MassOutput$WindowOpen<-MassOutput$WindowOpen+1

p3<-ggplot(MassOutput, aes(x=WindowClose, y=WindowOpen, fill=deltaAICc))+geom_tile()+theme_bw()+ylab('Window Open (month)')+xlab('Window Close (month)')+
  scale_fill_continuous(low = "red", high = "grey")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.tag.position = c(0.1, 0.95),legend.position = c(0.8, 0.2))+
  labs(fill="Delta AICc", tag="c")

p3

grid.arrange(p1, p2, nrow=1)
plotdelta(dataset = MassOutput)

check2<-(MassWin_Rela[[1]]$BestModelData)

mod1<-lmer(yvar~year+(1|site), data=check2)

rsq(mod1)
summary(mod1)


newdata <- with(check2, expand.grid(year=as.factor(seq(2013, 2017, 1)), site="Beach Haven",  climate=seq(11, 15, .5)))
newdata<-newdata[newdata$year!='2014',]
fit<-predict(mod1, newdata, type="response")

newdata<-cbind(newdata, fit)


p<-ggplot(newdata, aes(x=climate, y=fit))+geom_line(size=1)
p+facet_grid(.~year)+theme_bw()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab(bquote('Predicted shoots / '~(m^2)))+xlab(expression(paste("Average temperature (", degree,"C) for 2 month window beginning 2 months prior to sampling")))+
  theme(panel.spacing.x = unit(4, "mm"))+ scale_color_discrete(name =bquote('Leaf Area'~(dm^2)), l=40, c=55, guide = guide_legend(reverse=TRUE))







