#Test PK examples script for PKADVAN functions
rm(list=ls(all=TRUE))
graphics.off()

library(PKADVAN)
library(plyr)
library(doBy)
library(ggplot2)

#--------------------------------------------------------------
#Customize ggplot2 theme - R 2.15.3
theme_bw2 <- theme_set(theme_bw(base_size = 20))
theme_bw2 <- theme_update(plot.margin = unit(c(1.5,1.5,3,1.5), "lines"),
                          axis.title.x=element_text(size = 18, vjust = 0),
                          axis.title.y=element_text(size = 18, vjust = 0, angle = 90),
                          strip.text.x = element_text(size=12, margin = margin(3, 0, 3, 0)),
                          strip.text.y=element_text(size = 16, angle = 90))

#--------------------------------------------------------------
#Confidence intervals - from function utility

CI90lo <- function(x) quantile(x, probs=0.05,na.rm=T)
CI90hi <- function(x) quantile(x, probs=0.95,na.rm=T)

CI95lo <- function(x) quantile(x, probs=0.025,na.rm=T)
CI95hi <- function(x) quantile(x, probs=0.975,na.rm=T)

####################################################
## Three simple Steps for performing simulations:
# 1.  Generate a simulation data frame with the individual PK parameters and include any covariate effects on the PK parameters.
# 2.  Run the "PKADVAN" function of the selected model: the "PKADVAN" function will return the simulated amounts in each compartment of the PK system and IPREDs for the central compartment.
# 3.  Add residul un explained variability to the IPREDs.

######################################################################
########      Basic Pharmacokinetic models    ########################
######################################################################
# 1.  Iv bolus models (1,2, and 3 compartment models)
# 2.  IV infusion models (1,2, and 3 compartment models)
# 3.  First order absorption models(1,2, and 3 compartment models)

#====================
# IV bolus Models
#====================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(seq(0,72,4))
tlast <- 96
#set number of subjects
nsub <- 1000
ID  <- 1:nsub
#Make dataframe: CLCR: is creatinine clearance
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimes))),"AMT"=0,"MDV"=0,"CLCR"=90)
df$CLCR[df$ID >= 0.5*nsub] <- 70

doserows <- subset(df, TIME %in% dosetimes)

#Dose = 100 mg, Dose 1  at time 0
doserows$AMT <- 100

#Add back dose information
df <- rbind(df,doserows)
df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
df <- subset(df, (TIME==0 & AMT==0)==FALSE)      # remove the row that has a TIME=0 and AMT=0

#--------------------------------------------
# 1 compartment-IV Bolus via PKADVAN package
#--------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Define population PK parameters for 1-compartment model
CLpop <- 2       # clearance
Vpop  <- 10      # central volume of distribution

#Modify df for ADVAN calculation and include any covariate effects on the PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V  <- Vpop

#------------------------
# Apply PKADVAN function
#-----------------------
simdf <- ddply(dfadvan, .(ID), OneCompIVbolus)

#Add residual unexplained variability (within subject variability)
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)

head(simdf)
tail(simdf)

plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#Check processing time. Depends on PK sampling duration!
#system.time(ddply(dfadvan, .(ID), OneCompIVbolus))

#--------------------------------------------
# 2 compartment-IV Bolus via PKADVAN package
#--------------------------------------------
#Set population PK parameters for 2-compartment model
CLpop <- 2       # clearance
V1pop <- 10      # central volume of distribution
Qpop  <- 1       # inter-compartmental clearance
V2pop <- 25      # peripheral volume of distribution

#Modify df for ADVAN calculation and include any covariate effects on the PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*(dfadvan$CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1 <- V1pop
dfadvan$Q  <- Qpop
dfadvan$V2 <- V2pop

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), TwoCompIVbolus)
head(simdf)

#system.time(ddply(dfadvan, .(ID), TwoCompIVbolus))

#--------------------------------------------
# 3 compartment-IV Bolus via PKADVAN package
#--------------------------------------------
#Set population PK parameters for 1-compartment model
CLpop  <- 2          # clearance
V1pop  <- 10         # central volume of distribution
Q12pop <- 0.5        # inter-compartmental clearance (1)
V2pop  <- 30         # peripheral volume of distribution (1)
Q13pop <- 0.3        # inter-compartmental clearance (2)
V3pop  <- 40         # peripheral volume of distribution (2)

#Modify df for ADVAN calculation and include any covariate effects on the PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL  <- CLpop*(dfadvan$CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1  <- V1pop
dfadvan$Q2 <- Q12pop
dfadvan$V2  <- V2pop
dfadvan$Q3 <- Q13pop
dfadvan$V3  <- V3pop

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), ThreeCompIVbolus)
head(simdf)

#?processing time
system.time(ddply(dfadvan, .(ID), ThreeCompIVbolus))

#==================================
# First-order absorption Models
#==================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(seq(0,48,12))
tlast <- 72

#set number of subjects
nsub <- 1000
ID <- 1:nsub

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimes))),"AMT"=0,"MDV"=0,"CLCR"=90)
df$CLCR[df$ID <= 0.5*nsub] <- 70

doserows <- subset(df, TIME%in%dosetimes)

#Dose = 100 mg. It can be any arbitrary dose
doserows$AMT <- 100
doserows$MDV <- 1

#Add back dose information
df <- rbind(df,doserows)
df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#---------------------------------------------------------------
# 1 compartment-first order absorption via PKADVAN package
#---------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVF1 <- rnorm(nsub, mean = 0, sd = 0.2)  #BSV on F

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Define population PK parameters for 1-compartment model
CLpop <- 2       # clearance
Vpop  <- 10      # central volume of distribution
KApop <- 0.4     # first-order absorption rate constant
F1pop <- 1       # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V  <- Vpop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop*exp(BSVF1)

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), OneCompFirstOrderAbs)
head(simdf)

#Add residual error model
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)


plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#?processing time. Depends on PK sampling duration
#system.time(ddply(dfadvan, .(ID), OneCompFirstOrderAbs))


#----------------------------------------------------------------
# 2 compartment-first order absorption  via PKADVAN package
#----------------------------------------------------------------
#Set population PK parameters for 2-compartment model
CLpop <- 2       # clearance
V2pop  <- 10     # central volume of distribution
Qpop  <- 1       # inter-compartmental clearance
V3pop <- 25      # peripheral volume of distribution
KApop <- 0.4     # first-order absorption rate constant
F1pop <- 1       # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df

#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*(dfadvan$CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V2 <- V2pop
dfadvan$Q  <- Qpop
dfadvan$V3 <- V3pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), TwoCompFirstOrderAbs)
head(simdf)

#system.time(ddply(dfadvan, .(ID), TwoCompFirstOrderAbs))

#---------------------------------------------------------------
# 3 compartment-first order absorption via PKADVAN package
#---------------------------------------------------------------
#Set population PK parameters for 3-compartment model
CLpop <- 2          # clearance
V2pop <- 10         # central volume of distribution
Q3pop <- 0.5        # inter-compartmental clearance (1)
V3pop <- 30         # peripheral volume of distribution (1)
Q4pop <- 0.3        # inter-compartmental clearance (2)
V4pop <- 40         # peripheral volume of distribution (2)

KApop <- 0.4        # first-order absorption rate constant
F1pop <- 0.5          # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*(dfadvan$CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V2 <- V2pop
dfadvan$Q3 <- Q3pop
dfadvan$V3 <- V3pop
dfadvan$Q4 <- Q4pop
dfadvan$V4 <- V4pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), ThreeCompFirstOrderAbs)
head(simdf)

#system.time(ddply(dfadvan, .(ID), ThreeCompFirstOrderAbs))


#=======================================================================================
#                                  IV infusion
#=======================================================================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
#number of subjects
nsub <- 1000
ID <- 1:nsub

tlast <- 72
#Set dose records
dosetimessim <- c(0,20.5)

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimessim))),"AMT"=0,"RATE"=0,"CLCR"=90)

df$CLCR[df$ID>=0.5*nsub] <- 70

#subset doserows
doserowssim <- subset(df, TIME%in%dosetimessim)

#Dose & RATE can be arbitrary
doserowssim$AMT  <- 100
doserowssim$RATE <- 4

#Add back dose information
df <- rbind(df,doserowssim)
df <- df[order(df$TIME,-df$AMT),]       # arrange df by TIME with the AMT from high to low
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#----------------------------------------------------
# 1 compartment-IV infusion via PKADVAN package
#----------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVV  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.09)

#Define population PK parameters for 1-compartment model
CLpop <- 2       # clearance
Vpop  <- 10      # central volume of distribution

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V  <- Vpop*exp(BSVV)

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), OneCompIVinfusion)
head(simdf)

#Add residual unexplained variability (within subject variability)
#For example; additive error model
simdf$DV <- simdf$IPRED+PROP

#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#?provcessing time
#system.time(simdf <- ddply(dfadvan, .(ID), OneCompIVinfusion))


#--------------------------------------------------
# 2 compartment-IV infusion via PKADVAN package
#--------------------------------------------------
#Set population PK parameters for 1-compartment model
CLpop <- 2       # clearance
V1pop <- 10      # central volume of distribution
Qpop  <- 1       # inter-compartmental clearance
V2pop <- 25      # peripheral volume of distribution

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*(dfadvan$CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1 <- V1pop
dfadvan$Q  <- Qpop
dfadvan$V2 <- V2pop

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), TwoCompIVinfusion)

#system.time(simdf <- ddply(dfadvan, .(ID), TwoCompIVinfusion))

#-------------------------------------------------
# 3 compartment-IV infusion via PKADVAN package
#-------------------------------------------------
#Set population PK parameters for 1-compartment model
CLpop  <- 2          # clearance
V1pop  <- 10         # central volume of distribution
Q12pop <- 0.5        # inter-compartmental clearance (1)
V2pop  <- 30         # peripheral volume of distribution (1)
Q13pop <- 0.3        # inter-compartmental clearance (2)
V3pop  <- 40         # peripheral volume of distribution (2)

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL  <- CLpop*(dfadvan$CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1  <- V1pop
dfadvan$Q2 <- Q12pop
dfadvan$V2  <- V2pop
dfadvan$Q3 <- Q13pop
dfadvan$V3  <- V3pop

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), ThreeCompIVinfusion)
#system.time(simdf <- ddply(dfadvan, .(ID), ThreeCompIVinfusion))

#=======================================================================================
#                                 Transit Absorption models Models
#=======================================================================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(0)
tlast <- 36

#set number of subjects
nsub <- 1000
ID <- 1:nsub

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimes))),"AMT"=0,"MDV"=0,"CLCR"=90)
df$CLCR[df$ID <= 0.5*nsub] <- 70

doserows <- subset(df, TIME%in%dosetimes)

#Dose = 100 mg. It can be any arbitrary dose
doserows$AMT <- 100000 #ng
doserows$MDV <- 1

#Add back dose information
df <- rbind(df,doserows)
df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#----------------------------------------------------------------
# 1 compartment-1-transit  via PKADVAN package
#----------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVV2  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.2)

#Set population PK parameters for 2-compartment model
CLpop  <- 129       # clearance
V2pop  <- 861       # central volume of distribution
KTRpop <- 0.95      # first-order absorption rate constant
F1pop <- 0.75       # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
dfadvan$V2 <- V2pop*exp(BSVV2)
dfadvan$KTR <- KTRpop
dfadvan$F1 <- F1pop

#Apply the ADVAN function
simdf1 <- ddply(dfadvan, .(ID), OneCompOneTransit)
simdf1$TRANSIT <- "1-transit"

#Add residual unexplained variability (within subject variability)
#For example; additive error model
simdf1$DV <- simdf1$IPRED*(1+PROP)

#----------------------------------------------------------------
# 1 compartment-2-transit  via PKADVAN package
#----------------------------------------------------------------
simdf2 <- ddply(dfadvan, .(ID), OneCompTwoTransit)
simdf2$TRANSIT <- "2-transit"
#----------------------------------------------------------------
# 1 compartment-3-transit  via PKADVAN package
#----------------------------------------------------------------
simdf3 <- ddply(dfadvan, .(ID), OneCompThreeTransit)
simdf3$TRANSIT <- "3-transit"
#----------------------------------------------------------------
# 1 compartment-4-transit  via PKADVAN package
#----------------------------------------------------------------
simdf4 <- ddply(dfadvan, .(ID), OneCompFourTransit)
simdf4$TRANSIT <- "4-transit"

#Plot all!
simdf1 <- subset(simdf1, select = c(ID,TIME,IPRED,TRANSIT))
simdf2 <- subset(simdf2, select = c(ID,TIME,IPRED,TRANSIT))
simdf3 <- subset(simdf3, select = c(ID,TIME,IPRED,TRANSIT))
simdf4 <- subset(simdf4, select = c(ID,TIME,IPRED,TRANSIT))

simdfall <- rbind(simdf1,simdf2,simdf3,simdf4)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             ))
plotobj <- ggplot(data=simdfall)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED,colour=TRANSIT),fun.y=median, geom="line",  size=1)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

plotobj <- plotobj + facet_wrap(~TRANSIT)
plotobj

#system.time(simdf <- ddply(dfadvan, .(ID), OneCompOneTransit))
#system.time(simdf <- ddply(dfadvan, .(ID), OneCompTwoTransit))
#system.time(simdf <- ddply(dfadvan, .(ID), OneCompThreeTransit))
#system.time(simdf <- ddply(dfadvan, .(ID), OneCompFourTransit))

#----------------------------------------------------------------
# 2 compartment-1-transit  via PKADVAN package
#----------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVV2  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.2)

#Set population PK parameters for 2-compartment model
CLpop  <- 129       # clearance
V2pop  <- 861     # central volume of distribution
Qpop  <- 153       # inter-compartmental clearance
V3pop <- 2340      # peripheral volume of distribution
KTRpop <- 2.05     # first-order absorption rate constant
F1pop <- 0.75       # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
dfadvan$V2 <- V2pop*exp(BSVV2)
dfadvan$Q  <- Qpop
dfadvan$V3 <- V3pop
dfadvan$KTR <- KTRpop
dfadvan$F1 <- F1pop

#Apply PKADVAN functions
simdf1 <- ddply(dfadvan, .(ID), TwoCompOneTransit)
simdf1$TRANSIT <- "1-transit"
head(simdf)

#------------------------
# 2 compartment-2-transit
#------------------------
simdf2 <- ddply(dfadvan, .(ID), TwoCompTwoTransit)
simdf2$TRANSIT <- "2-transit"
head(simdf)

#------------------------
# 2 compartment-3-transit
#------------------------
simdf3 <- ddply(dfadvan, .(ID), TwoCompThreeTransit)
simdf3$TRANSIT <- "3-transit"
head(simdf)

#------------------------
# 2 compartment-4-transit
#------------------------
simdf4 <- ddply(dfadvan, .(ID), TwoCompFourTransit)
simdf4$TRANSIT <- "4-transit"
head(simdf)

#Plot all!
simdf1 <- subset(simdf1, select = c(ID,TIME,IPRED,TRANSIT))
simdf2 <- subset(simdf2, select = c(ID,TIME,IPRED,TRANSIT))
simdf3 <- subset(simdf3, select = c(ID,TIME,IPRED,TRANSIT))
simdf4 <- subset(simdf4, select = c(ID,TIME,IPRED,TRANSIT))

simdfall <- rbind(simdf1,simdf2,simdf3,simdf4)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
))
plotobj <- ggplot(data=simdfall)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED,colour=TRANSIT),fun.y=median, geom="line",  size=1)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

plotobj <- plotobj + facet_wrap(~TRANSIT)
plotobj

#system.time(simdf <- ddply(dfadvan, .(ID), TwoCompOneTransit))
#system.time(simdf <- ddply(dfadvan, .(ID), TwoCompTwoTransit))
#system.time(simdf <- ddply(dfadvan, .(ID), TwoCompThreeTransit))
#ystem.time(simdf <- ddply(dfadvan, .(ID), TwoCompFourTransit))


#===================================================================================
#     First-order absorption first-order formation metabolite models
#===================================================================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(seq(0,24,12))
tlast <- 48

#set number of subjects
nsub <- 1000
ID <- 1:nsub

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimes))),"AMT"=0,"MDV"=0,"CLCR"=90)
df$CLCR[df$ID <= 0.5*nsub] <- 70

doserows <- subset(df, TIME%in%dosetimes)

#Dose = 100 mg. It can be any arbitrary dose
doserows$AMT <- 500
doserows$MDV <- 1

#Add back dose information
df <- rbind(df,doserows)
df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#----------------------------------------------------------------------------------------------------
# 3 compartment-first order absorption with 1-compartment metabolite model via PKADVAN package
#----------------------------------------------------------------------------------------------------
# This example is an allustration on the possibility to extend ADVA-style equations to parent-metabolite models
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
#Parent BSV
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-parent
BSVV2  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-parent

#Metabolite BSV
BSVCLM <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-metabolite


#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.2) #parent
PROPMET <- rnorm(nsub, mean = 0, sd = 0.15) #metabolite

#Set population PK parameters for 3-compartment parent model
CLpop <- 2          # clearance
V2pop <- 10         # central volume of distribution
Q3pop <- 0.5        # inter-compartmental clearance (1)
V3pop <- 30         # peripheral volume of distribution (1)
Q4pop <- 0.3        # inter-compartmental clearance (2)
V4pop <- 40         # peripheral volume of distribution (2)
KApop <- 0.4        # first-order absorption rate constant
F1pop <- 1          # bioavailability

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 3        # Clearance of the metabolite
VpopMet  <- 10       # Volume of distribution of th metabolite

#Set fraction of parent drug converted into the metabolite
FR = 0.85

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
dfadvan$FR <- FR
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
dfadvan$V2 <- V2pop*exp(BSVV2)
dfadvan$Q3 <- Q3pop
dfadvan$V3 <- V3pop
dfadvan$Q4 <- Q4pop
dfadvan$V4 <- V4pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop
dfadvan$CLM <- CLpopMet**exp(BSVCLM)*(dfadvan$CLCR/100) #creatinine clearance (CLCR) added as a covariate on CLM
dfadvan$VM <- VpopMet

#Apply PKADVAN function to each ID in df
simdf <- ddply(dfadvan, .(ID), ThreeCompFirstOrderAbsOneCompMetab)
head(simdf)
tail(simdf)

#Add residual unexplained variability (within subject variability)
#For example; proportional error model
simdf$DVP <- simdf$IPREDP*(1+PROP)
simdf$DVM <- simdf$IPREDM*(1+PROPMET)


#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),geom="ribbon",fill="red", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),geom="ribbon", fill="red",fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj


#===========================================================================================
#                     infusion first-order formation Metabolite models
#===========================================================================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
#number of subjects
nsub <- 1000
ID <- 1:nsub

tlast <- 72
#Set dose records
dosetimessim <- c(0,20.5)

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimessim))),"AMT"=0,"RATE"=0,"CLCR"=90)

df$CLCR[df$ID>=0.5*nsub] <- 70

#subset doserows
doserowssim <- subset(df, TIME%in%dosetimessim)

#Dose & RATE can be arbitrary
doserowssim$AMT  <- 100
doserowssim$RATE <- 4

#Add back dose information
df <- rbind(df,doserowssim)
df <- df[order(df$TIME,-df$AMT),]       # arrange df by TIME with the AMT from high to low
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#---------------------------------------------------------------
# 3 compartment-IV infusion with 1-compartment metabolite model via PKADVAN package
#---------------------------------------------------------------
# This example is an allustration on the possibility to extend ADVA-style equations to parent-metabolite models
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
#Parent BSV
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-parent
BSVV  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-parent

#Metabolite BSV
BSVCLM <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-metabolite


#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.2) #parent
PROPMET <- rnorm(nsub, mean = 0, sd = 0.15) #metabolite

# Set population PK parameters for parent 3-compartment model
CLpop  <- 2          # clearance
V1pop  <- 10         # central volume of distribution
Q2pop <- 0.5        # inter-compartmental clearance (1)
V2pop  <- 30         # peripheral volume of distribution (1)
Q3pop <- 0.3        # inter-compartmental clearance (2)
V3pop  <- 40         # peripheral volume of distribution (2)

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 3        # Clearance of the metabolite
VpopMet  <- 20       # Volume of distribution of th metabolite

#Set fraction of parent drug converted into the metabolite
FR = 0.75

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
dfadvan$FR <- FR
#Calculate group parameter values including any covariate effects
dfadvan$CL  <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1  <- V1pop#*exp(BSVV)
dfadvan$Q2 <-  Q2pop
dfadvan$V2  <- V2pop
dfadvan$Q3 <-  Q3pop
dfadvan$V3  <- V3pop
dfadvan$CLM <- CLpopMet*exp(BSVCLM)*(dfadvan$CLCR/100)
dfadvan$VM <- VpopMet

#Apply the ADVAN function for each ID in df
simdf <- ddply(dfadvan, .(ID), ThreeCompIVinfusionOneCompMetab)
head(simdf)
tail(simdf)

#Add residual unexplained variability (within subject variability)
#For example; proportional error model
simdf$DVP <- simdf$IPREDP*(1+PROP)
simdf$DVM <- simdf$IPREDM*(1+PROPMET)


#DV
plotobj <- NULL
titletext <- expression(atop('Three compartment multiple IV infusion administration',
                             atop(italic("with discrete changes in creatinine clearance covariate values on parent systemic clearance"),
                                  italic("One compartment first-order formation for the metabolite"))))

plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),geom="ribbon",fill="red", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#IPRED
plotobj <- NULL
titletext <- expression(atop('Three compartment multiple IV infusion administration',
                             atop(italic("with discrete changes in creatinine clearance covariate values on parent systemic clearance"),
                                  italic("One compartment first-order formation for the metabolite"))))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),geom="ribbon", fill="red",fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

