#WCOP2016 simulation PK example using PKADVAN package
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
# 2.  Run the "PKADVAN" function for the selected model: the "PKADVAN" function will return the simulated amounts in each compartment of the PK system and IPREDs for the central compartment.
# 3.  Add residul un-explained variability to the IPREDs.

#########################################################################################
########      One compartment IV bolus with 1 comp metabolite    ########################
#########################################################################################

#================================
# Generate simulation data frame
#================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(seq(0,48,24))
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
doserows$MDV <- 1

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
BSVV  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Define population PK parameters for 1-compartment model
CLpop <- 0.2       # clearance
Vpop  <- 5      # central volume of distribution

#Modify df for ADVAN calculation and include any covariate effects on the PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V  <- Vpop*exp(BSVV)

#-----------------------
# Apply PKADVAN function
#-----------------------
simdf <- ddply(dfadvan, .(ID), OneCompIVbolus)

#Add residual unexplained variability (within subject variability)
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)

head(simdf)
tail(simdf)

#plotting
#Subset MDV==0
simdf <- subset(simdf, MDV==0)

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

#Check processing time. Depends on PK sampling duration!
system.time(ddply(dfadvan, .(ID), OneCompIVbolus))

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

system.time(ddply(dfadvan, .(ID), TwoCompIVbolus))

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

#================================
# First-order absorption Models
#================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(seq(0,48,12))
tlast <- 72

#Now define finer sample times for after a dose to capture Cmax
doseseq <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10)

#Use the outer product but with addition to expand this doseseq for all dosetimes
PKtimes <- outer(dosetimes,doseseq,"+")

#set number of subjects
nsub <- 1000
ID <- 1:nsub

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),PKtimes))),"AMT"=0,"MDV"=0,"CLCR"=90)
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
CLpop <- 0.5       # clearance
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

#plotting
#subset data for MDV ==0
simdf <- subset(simdf, MDV==0)

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

#?processing time. Depends on PK sampling duration
system.time(ddply(dfadvan, .(ID), OneCompFirstOrderAbs))

#----------------------------------------------------------------
# 2 compartment-first order absorption  via PKADVAN package
#----------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVF1 <- rnorm(nsub, mean = 0, sd = 0.2)  #BSV on F

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Set population PK parameters for 2-compartment model
CLpop  <- 0.5       # clearance
V2pop  <- 10     # central volume of distribution
Qpop  <- 1       # inter-compartmental clearance
V3pop <- 25      # peripheral volume of distribution
KApop <- 0.4     # first-order absorption rate constant
F1pop <- 1       # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df

#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V2 <- V2pop
dfadvan$Q  <- Qpop
dfadvan$V3 <- V3pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop*exp(BSVF1)

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), TwoCompFirstOrderAbs)
head(simdf)

#Add residual error model
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)

#plotting
#subset data for MDV ==0
simdf <- subset(simdf, MDV==0)

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

#?processing time
system.time(ddply(dfadvan, .(ID), TwoCompFirstOrderAbs))

#---------------------------------------------------------------
# 3 compartment-first order absorption via PKADVAN package
#---------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVF1 <- rnorm(nsub, mean = 0, sd = 0.2)  #BSV on F

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Set population PK parameters for 3-compartment model
CLpop <- 0.5        # clearance
V2pop <- 20          # central volume of distribution
Q3pop <- 0.5        # inter-compartmental clearance (1)
V3pop <- 30         # peripheral volume of distribution (1)
Q4pop <- 0.3        # inter-compartmental clearance (2)
V4pop <- 40         # peripheral volume of distribution (2)

KApop <- 0.4        # first-order absorption rate constant
F1pop <- 0.5          # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V2 <- V2pop
dfadvan$Q3 <- Q3pop
dfadvan$V3 <- V3pop
dfadvan$Q4 <- Q4pop
dfadvan$V4 <- V4pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop*exp(BSVF1)

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), ThreeCompFirstOrderAbs)
head(simdf)

#Add residual error model
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)

#plotting
#subset data for MDV ==0
simdf <- subset(simdf, MDV==0)

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

#?processing time
system.time(ddply(dfadvan, .(ID), ThreeCompFirstOrderAbs))

#=============================
# IV infusion Basic Models
#=============================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
#number of subjects
nsub <- 1000
ID <- 1:nsub

tlast <- 72
#Set dose records
dosetimessim <- c(0,20.5)

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(0.5,seq(0,tlast,1),dosetimessim))),"AMT"=0,MDV=0,"RATE"=0,"CLCR"=90)

df$CLCR[df$ID>=0.5*nsub] <- 70

#subset doserows
doserowssim <- subset(df, TIME%in%dosetimessim)

#Dose & RATE can be arbitrary
doserowssim$AMT  <- 100
doserowssim$RATE <- 4
doserowssim$MDV  <- 1

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

#plotting
#subset data for MDV==0
simdf <- subset(simdf, MDV==0)

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

#?provcessing time
system.time(simdf <- ddply(dfadvan, .(ID), OneCompIVinfusion))

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

#=======================================================================
#   First-order absorption first-order metabolite formation models
#=======================================================================
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
CLpop <- 0.5          # clearance
V2pop <- 20         # central volume of distribution
Q3pop <- 0.5        # inter-compartmental clearance (1)
V3pop <- 30         # peripheral volume of distribution (1)
Q4pop <- 0.3        # inter-compartmental clearance (2)
V4pop <- 40         # peripheral volume of distribution (2)
KApop <- 0.4        # first-order absorption rate constant
F1pop <- 1          # bioavailability

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 1        # Clearance of the metabolite
VpopMet  <- 15       # Volume of distribution of th metabolite

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
nsub <- 10000
ID <- 1:nsub

tlast <- 72
#Set dose records
dosetimessim <- c(0,20.5)

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimessim))),"AMT"=0,"RATE"=0,"CLCR"=90,"MDV"=0)

df$CLCR[df$ID>=0.5*nsub] <- 70

#subset doserows
doserowssim <- subset(df, TIME%in%dosetimessim)

#Dose & RATE can be arbitrary
doserowssim$AMT  <- 100
doserowssim$RATE <- 4
doserowssim$MDV  <- 1

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
BSVVM  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-metabolite

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
dfadvan$V1  <- V1pop*exp(BSVV)
dfadvan$Q2 <-  Q2pop
dfadvan$V2  <- V2pop
dfadvan$Q3 <-  Q3pop
dfadvan$V3  <- V3pop
dfadvan$CLM <- CLpopMet*exp(BSVCLM)*(dfadvan$CLCR/100)
dfadvan$VM <- VpopMet*exp(BSVVM)

#Apply the ADVAN function for each ID in df
simdf <- ddply(dfadvan, .(ID), ThreeCompIVinfusionOneCompMetab)
head(simdf)
tail(simdf)

#Add residual unexplained variability (within subject variability)
#For example; proportional error model
simdf$DVP <- simdf$IPREDP*(1+PROP)
simdf$DVM <- simdf$IPREDM*(1+PROPMET)

#plotting
#Subset MDV==0
simdf <- subset(simdf, MDV==0)

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

