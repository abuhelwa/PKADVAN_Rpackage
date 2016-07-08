#' @title Process simulations using 3-compartment IV-infusion model with 1-compartment first-order formation metabolite model.
#' @description \code{ThreeCompIVinfusionOneCompMetab} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments.
#' The data frame should have the following columns: \code{ID, TIME, AMT, RATE, CL, V1, Q2, V2, Q3, V3, CLM, VM, FR}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{RATE}:\tab is the infusion rate\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V1}:\tab is the central volume of distribution\cr
#' \code{Q2}:\tab is the inter-compartmental clearance (1)\cr
#' \code{V2}:\tab is the peripheral volume of distribution (1)\cr
#' \code{Q3}:\tab is the inter-compartmental clearance (2)\cr
#' \code{V3}:\tab is the peripheral volume of distribution (2)\cr
#' \code{CLM}:\tab is the central clearance of the metabolite\cr
#' \code{VM}:\tab is the  central volume of distribution of the metabolite\cr
#' \code{FR}:\tab is the  fraction of the parent drug converted into metabolite\cr
#' }
#'
#' @usage ThreeCompIVinfusionOneCompMetab(inputDataFrame)
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for: \code{ID, TIME, AMT, RATE, CL, V1, Q2, V2, Q3, V3, CLM, VM, FR}.\cr
#' @return The function calculates the parent and metabolite amounts in the respective compartments of the pharmacokinetic model.
#' This includes the amount in the parent central (\code{A1} & individual predicted parent concentrations, \code{IPREDP})
#' parent two peripheral \code{(A2, A3)}, and metabolite central (\code{AM} & individual predicted metabolite concentrations, \code{IPREDM}) compartments.
#' The function returns the output added to the \code{inputDataFrame}
#'
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),  the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), ThreeCompIVinfusionOneCompMetab)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{ThreeCompFirstOrderAbsOneCompMetab}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @export

#------------------------------------------------------------------------------------------------------------------------
# 3 compartment-IV infusion with 1 compartment first-order metabolite formation via ADVAN-style equations: RCppfunction
#-----------------------------------------------------------------------------------------------------------------------
ThreeCompIVinfusionOneCompMetab <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV,RATE, CL, V1, Q2, V2, Q3, V3, CLM, VM, FR,
#Returns a dataframe with populated columns for A1, A2, A3, AM, IPREDP, IPREDM

  #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
  k10 <- k12 <- k21 <- k20 <- k13 <- k31 <- k30 <- kmf <- kme <- NULL
  RATEALL <- TIME <- CLpop1 <- CLpop2 <- NULL

  #Sampling Times
  sampletimes <- inputDataFrame$TIME

  #Process infusion doses
  inputDataFrame <- ProcessInfusionDoses(inputDataFrame)

  #Calculate micro-rate constants-Parent
  FR = inputDataFrame$FR[1]
  inputDataFrame$CLpop1 <- inputDataFrame$CL*(1-FR)   # Clearance of the parent drug to outside the body
  inputDataFrame$CLpop2 <- inputDataFrame$CL*FR       # Clearance of the parent drug into the metabolite compartment

  #Calculate rate constants
  inputDataFrame$k10 <- inputDataFrame$CL/inputDataFrame$V1
  inputDataFrame$k12 <- inputDataFrame$Q2/inputDataFrame$V1
  inputDataFrame$k21 <- inputDataFrame$Q2/inputDataFrame$V2
  inputDataFrame$k20 <- 0
  inputDataFrame$k13 <- inputDataFrame$Q3/inputDataFrame$V1
  inputDataFrame$k31 <- inputDataFrame$Q3/inputDataFrame$V3
  inputDataFrame$k30 <- 0

  #Calculate micro-rate constants-Metabolite
  inputDataFrame$kmf <- inputDataFrame$CLpop2/inputDataFrame$V1  #Rate constant for metabolite formation
  inputDataFrame$kme <- inputDataFrame$CLM/inputDataFrame$VM	 #Rate constant for metabolite elimination

  #Add columns for amounts for marginals speed gain!
  inputDataFrame$A1 <- 0
  inputDataFrame$A2 <- 0
  inputDataFrame$A3 <- 0
  inputDataFrame$AM <- 0

  #set initial values in the compartments
  inputDataFrame$A1[inputDataFrame$TIME==0] <- 0                   # Amount in the central compartment at time zero.
  inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                   # Amount in the 1st peripheral compartment at time zero.
  inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                   # Amount in the 2nd peripheral compartment at time zero.
  inputDataFrame$AM[inputDataFrame$TIME==0] <- 0                   # Amount in the metabolite compartment at time zero.

  ThreeCompIVinfusionOneCompMetabCpp( inputDataFrame )

  #Remove end infusion time points
  inputDataFrame <- subset(inputDataFrame, (TIME%in%sampletimes))

  #Calculate IPRED for the central compartment
  inputDataFrame$IPREDP <- inputDataFrame$A1/inputDataFrame$V1  #Concentration of parent
  inputDataFrame$IPREDM <- inputDataFrame$AM/inputDataFrame$VM  #Concentration of metabolite

  #subset extra columns
  inputDataFrame <- subset(inputDataFrame, select=-c(k10,k12,k21,k20,k13,k31,k30,RATEALL,CLpop1,CLpop2,kmf,kme))

  #Return output
  inputDataFrame
}
