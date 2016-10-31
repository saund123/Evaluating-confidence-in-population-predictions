#############################################################################
#############Data manipulation and model fitting code########################
##Evaluating confidence in climate-based predictions of population change in 
##a migratory species. Global Ecology & Biogeography 25:1000-1012
##Authors: Saunders, SP, Ries, L, Oberhauser, KS, Zipkin, EF
##Data sources are proprietary. Please contact authors for more information
##Note: only monarch observations were used from datasets of all butterfly
##species observations in Ohio and Illinois monitoring networks
#############################################################################

#load libraries
library(jagsUI)
library(reshape)

#set working directory to where input files are located on user computer
setwd("C:/Users/")

#read in data 
#SurveyData file consists of file with observations by row and the following columns:
#Latitude, longitude, Site ID, Year, Week, Sum of Monarch Counts, PDI by week, GDDs accumulated
#from weeks 10-28, average GDD, GDD differential (latter 4 columns calculated outside of R)
mon = read.csv("SurveyData.csv", header=TRUE,sep=',',na.strings=T)
mon[1:10,]
names(mon)
str(mon)

#Look at data (checking all years/weeks/sites present)
uyear = sort(unique(mon$Year))
uweek= sort(unique(mon$Week))
usite = sort(unique(mon$SiteID))

#Reshape the count data.
#NOTE: This assumes columns called Year, Week, SiteID and SumofMonarch count
#Does not need to be in any specific order
junk.melt=melt(mon,id.var=c("Year", "Week", "SiteID"), measure.var="SumOfMonarch.count")
monarchs=cast(junk.melt, SiteID ~ Week ~ Year, sum)
#This organizes into tables that are sums of counts according to Week by SiteID for each year

#Include only spring/summer (e.g., weeks 5-28)
monarchs1=monarchs[,1:24,]
#keep weeks 5-28

#####################
###Adding covariates
#####################
#WEEK
uweek1=uweek[1:24]
#standardize uweek and make quadratic term
mweek=mean(uweek1)
sdweek=sd(uweek1)
week1<-(uweek1-mweek)/sdweek
week2<-week1*week1

#YEAR
myear=mean(uyear)
sdyear=sd(uyear)
year1<-(uyear-myear)/sdyear
year2<-year1*year1

#SURVEYING EFFORT
junk.melt2=melt(mon,id.var=c("Year", "Week", "SiteID"), measure.var="SumOfTotalTime")
effort=cast(junk.melt2, SiteID ~ Week ~ Year, sum)
#this creates tables of week by site according to year with 'sum' being sum of total time
effort = effort[,1:24,] 
#this eliminates weeks 29-35 as above
effort=effort/60
#turns the variable into hours

#################
#HABITAT COVARIATES
################

#SiteEffects
#This file is organized according to sites as rows and the following columns:
#Latitude, longitude, Site ID, percent open, average GDD
coord=read.csv("siteEffects.csv", header=T,sep=',',na.strings=T)

#Site effects--constant across years
#% open habitat (standardized)
mopen=mean(coord$Open)
sdopen=sd(coord$Open)
open1<-(coord$Open-mopen)/sdopen

#Avg GDD
mavggdd=mean(coord$avgGDD)
sdavggdd=sd(coord$avgGDD)
avggdd1<-(coord$avgGDD-mavggdd)/sdavggdd

#SITE BY YEAR effects
#NOTE: sites are in numerical order in this file
#File is organized with each row as site-year combinations; thus, columns are:
#Site ID, Year, accumulated GDD from weeks 10 to 28, and average PDI accumulated through week 28
siteyear = read.csv("SiteYear.csv", header=T,sep=',',na.strings=T)

#organize site by year
#report appropriate GDD and PDSI values according to siteID and year
sitegdd = matrix(0,nrow=length(usite), ncol=length(uyear))
sitepalm = matrix(0,nrow=length(usite), ncol=length(uyear))
for (tt in 1:length(uyear)) {
  for (j in 1:length(usite)) {
    a=which(siteyear$SiteID == usite[j] & siteyear$Year == uyear[tt])
    sitegdd[j,tt] = siteyear$AccGDD[a]  
    sitepalm[j,tt] = siteyear$AvgPDSI[a] 
  }} 

#standardize site GDDs
msitegdd=mean(as.matrix(sitegdd))
sdsitegdd=sd(as.vector(sitegdd))
sitegdd1=as.matrix((sitegdd-msitegdd)/sdsitegdd)

#standardize site drought indeces
msitepalm=mean(as.matrix(sitepalm))
sdsitepalm=sd(as.vector(sitepalm))
sitepalm1=as.matrix((sitepalm-msitepalm)/sdsitepalm)

#YEAR effects (TX GDD and drought)
#NOTE: This file is first year of study (first row) to last year of study (bottom row) with columns:
#Year, spring GDD, spring precipitation
yeareffects=read.csv("YearEffects.csv", header=T,sep=',',na.strings=T)

spGDD=yeareffects$SprGDD
spPrec=yeareffects$SprPrecipTX

#standardize spring GDD
mspGDD=mean(as.matrix(spGDD))
sdspGDD=sd(as.vector(spGDD))
spGDD1=as.vector((spGDD-mspGDD)/sdspGDD)

#standardize spring Precip
mspPrec=mean(as.matrix(spPrec))
sdspPrec=sd(as.vector(spPrec))
spPrec1=as.vector((spPrec-mspPrec)/sdspPrec)

#SITE BY WEEK BY YEAR effects
#NOTE: survey file doesn't need to be in order by site for this loop since it reads by siteID
gdd=array(0, dim=dim(monarchs1))
gdddiff=array(0, dim=dim(monarchs1))
pdsi=array(0, dim=dim(monarchs1))
for (tt in 1:length(uyear)) {
  for (k in 1:length(uweek1)) { 
    for (j in 1:length(usite)) {
      a=which(mon$SiteID==usite[j] & mon$Year==uyear[tt] & mon$Week==uweek1[k])
      if (length(a)==1)   {
        gdd[j,k,tt]=mon$GDDwk10NEW[a]  
        gdddiff[j,k,tt]=mon$GDDdiffNEW[a]
        pdsi[j,k,tt]=mon$PDSIwk[a]  }
    }}   }

#standardize GDDdiff
#NOTE: adjust weeks below as needed
mgdddiff=mean(gdddiff[,6:24,], na.rm=T)
sdgdddiff=sd(as.vector(gdddiff[,6:24,]), na.rm=T)
gdddiff1<-(gdddiff[,6:24,]-mgdddiff)/sdgdddiff     

#standardize drought index
#NOTE: adjust weeks below as needed
mpdsi=mean(pdsi[,6:24,], na.rm=T)
sdpdsi=sd(as.vector(pdsi[,6:24,]), na.rm=T)
pdsi1<-(pdsi[,6:24,]-mpdsi)/sdpdsi 

#use summer weeks
#NOTE: adjust weeks below as needed
uweek2=uweek1[6:24]
week11<-(uweek2-mean(uweek2))/sd(uweek2)

###############################################
### Specify model in BUGS language#############
##############################################

sink("GEBmodel.jags")
cat("
    model {

    #Priors
    a1  ~ dnorm(0,0.01)
    a2  ~ dnorm(0,0.01)
    a3  ~ dnorm(0,0.01)
    a4  ~ dnorm(0,0.01)
    a5  ~ dnorm(0,0.01)
    a6  ~ dnorm(0,0.01)
    a7  ~ dnorm(0,0.01)
    a8  ~ dnorm(0,0.01)
    a9  ~ dnorm(0,0.01)
    a10  ~ dnorm(0,0.01)
    a11  ~ dnorm(0,0.01)
    a12  ~ dnorm(0,0.01)
    a13  ~ dnorm(0,0.01)
    a14  ~ dnorm(0,0.01)
    a15  ~ dnorm(0,0.01)
    a16  ~ dnorm(0,0.01)
    a17  ~ dnorm(0,0.01)
    
    rprec ~ dgamma(0.01,0.01)
    
    #Negative binomial regression model
    for (j in 1:usite){
    for (k in 1:uweek){
    for (t in 1:uyear){

    y[j,k,t]  ~ dnegbin(p[j,k,t],rprec)

    p[j,k,t] <- rprec/(rprec + phi[j,k,t])

    log(phi[j,k,t]) <- log(effort[j,k,t])+ a1 + a2*week1[k] + a3*spPrec1[t] + a4*spPrec2[t] +
    a5*spGDD1[t] + a6*spGDD2[t] + a7*spPrec1[t]*week1[k] + a8*spGDD1[t]*week1[k] + 
    a9*gdddiff1[j,k,t] + a10*avggdd1[j] + a11*avggdd2[j] + a12*gdddiff1[j,k,t]*week1[k] +                                                
    a13*gdddiff1[j,k,t]*avggdd1[j]*week1[k] + a14*sitepalm1[j,t]+a15*sitepalm2[j,t] + 
    a16*sitepalm1[j,t]*week1[k] + a17*open1[j]

    }
    }
    }

    #Posterior predictive check for first 8 years as example (survey data included)
    for (j in 1:usite){
    for (t in 1:8){

    EPRes[j,t] <- (p[j,19,t]*rprec)/(1-rprec) #Expected count
    Denom[j,t] <- 1-rprec
    VarPRes[j,t] <- (p[j,19,t]*rprec)/pow(Denom[j,t],2) #Variance

    #Pearson's residual
    PRes[j,t] <- (y[j,19,t]-EPRes[j,t])/sqrt(VarPRes[j,t]+0.5)
    sq[j,t]<-pow(PRes[j,t],2)
    
    #Generation of ideal data set from estimated parameters (y.new)
    y.new[j,t] ~ dnegbin(p[j,19,t],rprec)
    PRes.new[j,t] <- (y.new[j,t]-EPRes[j,t])/sqrt(VarPRes[j,t]+0.5)
    sq.new[j,t]<-pow(PRes.new[j,t],2)

    } }

    #Calculating Bayesian p-value for first 8 years (to determine model fit)
    fit <- sum(sq[,])
    fit.new <- sum(sq.new[,])

    #Posterior predictive check for last 8 years (survey data omitted)
    for (j in 1:usite){
    for (t in 1:8){
    yy[j,t] ~ dnegbin(p[j,19,t+8],rprec)
    
    EPRes1[j,t] <- (p[j,19,t+8]*rprec)/(1-rprec) #Expected count
    Denom1[j,t] <- 1-rprec
    VarPRes1[j,t] <- (p[j,19,t+8]*rprec)/pow(Denom1[j,t],2) #Variance

    #Pearson's residual
    PRes1[j,t] <- (yy[j,t]-EPRes1[j,t])/sqrt(VarPRes1[j,t]+0.5)
    sq1[j,t]<-pow(PRes1[j,t],2)
    
    #Generation of 'ideal' data set using observed counts (y.new1)
    y.new1[j,t] ~ dnegbin(p[j,19,t+8],rprec)
    PRes.new1[j,t] <- (y.new1[j,t]-EPRes1[j,t])/sqrt(VarPRes1[j,t]+0.5)
    sq.new1[j,t]<-pow(PRes.new1[j,t],2)

    } }

    #Calculating Bayesian p-values for each of last 8 years (years to predict)
    for (t in 1:8) {
    fit.year[t] <- sum(sq1[,t])
    fit.new.year[t] <- sum(sq.new1[,t])

    }

    #Calculating Bayesian p-value for last 8 years combined (years to predict)
    fit1 <- sum(sq1[,])
    fit.new1 <- sum(sq.new1[,])
  
    }

    ", fill=TRUE)
sink()

#PREP BUGS DATA
bugsdata<-list(uyear=length(uyear), usite=length(usite), uweek=length(uweek2),
               y=monarchs1[,6:24,], week1=week11, effort=effort[,6:24,], 
               gdddiff1=gdddiff1, sitepalm1=sitepalm1,
               spGDD1=spGDD1, open1=open1, spPrec1=spPrec1, 
               sitepalm2=sitepalm1*sitepalm1, spGDD2=spGDD1*spGDD1,
               spPrec2=spPrec1*spPrec1, avggdd1=avggdd1, avggdd2=avggdd1*avggdd1) 

#inits for NEG BINOM model
inits<-function(){
  list(rprec=rgamma(1,1))
}

#parameters to monitor
parameters<-c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 
              'a10', 'a11', 'a12', 'a13', 'a14', 'a15', 'a16', 'a17', 'fit1', 'fit.new1')

###################################################################
##RUN BUGS MODEL
###################################################################

library(jagsUI)
output<-jags(data = bugsdata,
              inits = inits,
              parameters.to.save = parameters,
              model.file = 'GEBmodel.jags',
              n.chains = 3,
              n.adapt = 100,
              n.iter = 4000,
              n.burnin = 1000,
              n.thin = 3,
              #parallel=TRUE
)

print(output,digits=4)

whiskerplot(output,parameters=c("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15","a16","a17"))

###########Evaluation of fit####################################################
mean(output$sims.list$fit.new1>output$sims.list$fit1) #this should be near 0.5
mean(output$mean$fit1)/mean(output$mean$fit.new1) #this should be near 1
pp.check(output,actual="fit1",new="fit.new1")

####################################################################
####Replace last year of data with NA to determine predictive ability
####################################################################

#NOTE: change this according to number of years to omit as NAs
monarchs6=monarchs1
monarchs6[,,16]<-NA

#######################################################################################
### Specify model in BUGS language
#The following model is for predicting last year of data; adjust as needed depending on what
#year is of interest to predict
#######################################################################################

sink("GEBmodel.predict.jags")
cat("
    model {

    a1  ~ dnorm(0,0.01)
    a2  ~ dnorm(0,0.01)
    a3  ~ dnorm(0,0.01)
    a4  ~ dnorm(0,0.01)
    a5  ~ dnorm(0,0.01)
    a6  ~ dnorm(0,0.01)
    a7  ~ dnorm(0,0.01)
    a8  ~ dnorm(0,0.01)
    a9  ~ dnorm(0,0.01)
    a10  ~ dnorm(0,0.01)
    a11  ~ dnorm(0,0.01)
    a12  ~ dnorm(0,0.01)
    a13  ~ dnorm(0,0.01)
    a14  ~ dnorm(0,0.01)
    a15  ~ dnorm(0,0.01)
    a16  ~ dnorm(0,0.01)
    a17  ~ dnorm(0,0.01)
    
    rprec ~ dgamma(0.01,0.01)
    
    for (j in 1:usite){
    for (k in 1:uweek){
    for (t in 1:uyear){
    
    y[j,k,t]  ~ dnegbin(p[j,k,t],rprec)
    
    p[j,k,t] <- rprec/(rprec + phi[j,k,t])
    
    log(phi[j,k,t]) <- log(effort[j,k,t]) + a1 + a2*week1[k] + a3*spPrec1[t] + a4*spPrec2[t] + 
    a7*spPrec1[t]*week1[k] + a5*spGDD1[t] + a6*spGDD2[t] + a8*spGDD1[t]*week1[k] + 
    a9*gdddiff1[j,k,t] + a10*avggdd1[j] + a11*avggdd2[j] + a12*gdddiff1[j,k,t]*week1[k] +                               
    a13*gdddiff1[j,k,t]*avggdd1[j]*week1[k] + a14*sitepalm1[j,t] + a15*sitepalm2[j,t] + 
    a16*sitepalm1[j,t]*week1[k] + a17*open1[j]
    
    }
    }
    }
    
    for (j in 1:usite){
    for (t in 1:15){
    
    EPRes[j,t] <- (p[j,19,t]*rprec)/(1-rprec)
    Denom[j,t] <- 1-rprec
    VarPRes[j,t] <- (p[j,19,t]*rprec)/pow(Denom[j,t],2)
    PRes[j,t] <- (y[j,19,t]-EPRes[j,t])/sqrt(VarPRes[j,t]+0.5)
    sq[j,t]<-pow(PRes[j,t],2)
    
    y.new[j,t] ~ dnegbin(p[j,19,t],rprec)
    PRes.new[j,t] <- (y.new[j,t]-EPRes[j,t])/sqrt(VarPRes[j,t]+0.5)
    sq.new[j,t]<-pow(PRes.new[j,t],2)
    } 
    }
    
    fit <- sum(sq[,])
    fit.new <- sum(sq.new[,])
    
    for (j in 1:usite){
    
    yy[j] ~ dnegbin(p[j,19,16],rprec)
    
    EPRes1[j] <- (p[j,19,16]*rprec)/(1-rprec)
    Denom1[j] <- 1-rprec
    VarPRes1[j] <- (p[j,19,16]*rprec)/pow(Denom1[j],2)
    PRes1[j] <- (yy[j]-EPRes1[j])/sqrt(VarPRes1[j]+0.5)
    sq1[j]<-pow(PRes1[j],2)
    
    y.new1[j] ~ dnegbin(p[j,19,16],rprec)
    PRes.new1[j] <- (y.new1[j]-EPRes1[j])/sqrt(VarPRes1[j]+0.5)
    sq.new1[j]<-pow(PRes.new1[j],2)
    
    }
    
    fit1 <- sum(sq1[])
    fit.new1 <- sum(sq.new1[])

    }

    ", fill=TRUE)
sink()

#PREP WINBUGS DATA
#change years in yy according to those to omit
#change y dataset to monarchs2, monarchs3, etc. depending on what year omitted in dataset
bugsdata<-list(uyear=length(uyear), usite=length(usite), uweek=length(uweek2),
               yy=(monarchs1[,24,16]), 
               y=monarchs6[,6:24,],
               week1=week11, effort=effort[,6:24,], 
               gdddiff1=gdddiff1, sitepalm1=sitepalm1,
               spGDD1=spGDD1, open1=open1, spPrec1=spPrec1, 
               sitepalm2=sitepalm1*sitepalm1, spGDD2=spGDD1*spGDD1,
               spPrec2=spPrec1*spPrec1, avggdd1=avggdd1, avggdd2=avggdd1*avggdd1) 

#inits for NEG BINOM model
inits<-function(){
  list(rprec=rgamma(1,1))
}

#add new parameters to monitor
parameters<-c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 
              'a10', 'a11', 'a12', 'a13', 'a14', 'a15', 'a16', 'a17','fit','fit.new','fit1','fit.new1') 

###################################################################
##RUN BUGS MODEL
###################################################################

#NOTE: change model input file according to years to omit
out.predict<-jags(data = bugsdata,
                    inits = inits,
                    parameters.to.save = parameters,
                    model.file = 'GEBmodel.predict.jags',
                    n.chains = 3,
                    n.adapt = 100,
                    n.iter = 4000,
                    n.burnin = 1000,
                    n.thin = 3,
                    #parallel=TRUE
)

#checking fit of 15 years
mean(out.predict$sims.list$fit.new>out.predict$sims.list$fit) 
#checking fit of given year
mean(out.predict$sims.list$fit.new1>out.predict$sims.list$fit1)
