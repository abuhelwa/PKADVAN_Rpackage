#' Process simulations using 1-compartment first-order absorption model with 1-compartment first-order formation metabolite model.
#'
#' @description    \code{OneCompFirstOrderAbsOneCompMetab} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments.
#' The data frame should have the following columns: \code{ID, TIME, AMT, F1, KA, CL, V, CLM, VM, FR}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KA}:\tab is the absorption rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance of parent drug\cr
#' \code{V}:\tab is the central volume of distribution of parent drug\cr
#' \code{CLM}:\tab is the central clearance of the metabolite\cr
#' \code{VM}:\tab is the    central volume of distribution of the metabolite\cr
#' \code{FR}:\tab is the    fraction of the parent drug converted into metabolite\cr
#' }
#'
#' @usage OneCompFirstOrderAbsOneCompMetab(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, MDV, F1, KA, CL, V, CLM, VM, FR}.
#'
#' @return The function calculates the parent and metabolite amounts in the respective compartments of the pharmacokinetic model.
#' This includes the amount in the absorption \code{(A1)}, parent central (\code{A2} & individual predicted parent concentrations, \code{IPREDP})
#' and metabolite central (\code{AM} & individual predicted metabolite concentrations, \code{IPREDM}) compartments.
#' The function returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompFirstOrderAbsOneCompMetab)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompFirstOrderAbsOneCompMetab}}, \code{\link{ThreeCompFirstOrderAbsOneCompMetab}}
#' @seealso \code{\link{OneCompIVbolusOneCompMetab}}, \code{\link{TwoCompIVbolusOneCompMetab}}, \code{\link{ThreeCompIVbolusOneCompMetab}}
#' @seealso \code{\link{OneCompIVinfusionOneCompMetab}}, \code{\link{TwoCompIVinfusionOneCompMetab}}, \code{\link{ThreeCompIVinfusionOneCompMetab}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#-----------------------------------------------------------------------------------------------------------------
# First-order absorption: 1 compartment parent with 1 compartment first-order metabolite formation: RCpp function
#-----------------------------------------------------------------------------------------------------------------
OneCompFirstOrderAbsOneCompMetab <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT, MDV, CL, V, CLM, VM, FR, KA & F1
    #Returns a dataframe with populated columns for A1, A2, AM, IPREDP, IPREDM

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- kmf <- kme <- NULL
    CLpop1 <- CLpop2 <- NULL

    #Calculate micro-rate constants
    FR = inputDataFrame$FR[1]
    inputDataFrame$CLpop1 <- inputDataFrame$CL*(1-FR)	# Clearance of the parent drug to outside the body
    inputDataFrame$CLpop2 <- inputDataFrame$CL*FR		# Clearance of the parent drug into the metabolite compartment

    #Calculate micro-rate constants-Parent
    inputDataFrame$k20    <- inputDataFrame$CLpop1/inputDataFrame$V

    #Calculate micro-rate constants-Metabolite
    inputDataFrame$kmf <- inputDataFrame$CLpop2/inputDataFrame$V    #Rate constant for metabolite formation
    inputDataFrame$kme <- inputDataFrame$CLM/inputDataFrame$VM	 #Rate constant for metabolite elimination

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]	#Amount of parent in the absorption compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0																	#Amount of parent in the central compartment at time zero.
    inputDataFrame$AM[inputDataFrame$TIME==0] <- 0                                     								#Amount of metabolite in the metabolite compartment at time zero.

    OneCompFirstOrderAbsOneCompMetabCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPREDP <- inputDataFrame$A2/inputDataFrame$V     #Concentration of parent
    inputDataFrame$IPREDM <- inputDataFrame$AM/inputDataFrame$VM    #Concentration of metabolite

    #subset unnecessary columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,CLpop1,CLpop2,kmf,kme))

    #Return output
    inputDataFrame
}

#' Process simulations using 2-compartment first-order absorption model with 1-compartment first-order formation metabolite model.
#'
#' @description    \code{TwoCompFirstOrderAbsOneCompMetab} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments.
#' The data frame should have the following columns: \code{ID, TIME, AMT, F1, KA, CL, V2, Q3, V3, CLM, VM, FR}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KA}:\tab is the absorption rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance of parent drug\cr
#' \code{V2}:\tab is the central volume of distribution of parent drug\cr
#' \code{Q3}:\tab    is the inter-compartmental clearance of parent drug\cr
#' \code{V3}:\tab is the peripheral volume of distribution of parent drug\cr
#' \code{CLM}:\tab is the central clearance of the metabolite\cr
#' \code{VM}:\tab is the    central volume of distribution of the metabolite\cr
#' \code{FR}:\tab is the    fraction of the parent drug converted into metabolite\cr
#' }
#'
#' @usage TwoCompFirstOrderAbsOneCompMetab(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, MDV, F1, KA, CL, V2, Q3, V3, CLM, VM, FR}.
#'
#' @return The function calculates the parent and metabolite amounts in the respective compartments of the pharmacokinetic model.
#' This includes the amount in the absorption \code{(A1)}, parent central (\code{A2} & individual predicted parent concentrations, \code{IPREDP})
#' parent peripheral \code{(A3)}, and metabolite central (\code{AM} & individual predicted metabolite concentrations, \code{IPREDM}) compartments.
#' The function returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompFirstOrderAbsOneCompMetab)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompFirstOrderAbsOneCompMetab}}, \code{\link{ThreeCompFirstOrderAbsOneCompMetab}}
#' @seealso \code{\link{OneCompIVbolusOneCompMetab}}, \code{\link{TwoCompIVbolusOneCompMetab}}, \code{\link{ThreeCompIVbolusOneCompMetab}}
#' @seealso \code{\link{OneCompIVinfusionOneCompMetab}}, \code{\link{TwoCompIVinfusionOneCompMetab}}, \code{\link{ThreeCompIVinfusionOneCompMetab}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#-----------------------------------------------------------------------------------------------------------------
# First-order absorption: 2 compartment parent with 1 compartment first-order metabolite formation: RCpp function
#-----------------------------------------------------------------------------------------------------------------
TwoCompFirstOrderAbsOneCompMetab <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT, MDV, CL, V2, Q3, V3, CLM, VM, FR, KA & F1
    #Returns a dataframe with populated columns for A1, A2, A3, AM, IPREDP, IPREDM

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- kmf <- kme <- NULL
    CLpop1 <- CLpop2 <- NULL

    #Calculate micro-rate constants
    FR = inputDataFrame$FR[1]
    inputDataFrame$CLpop1 <- inputDataFrame$CL*(1-FR)     # Clearance of the parent drug to outside the body
    inputDataFrame$CLpop2 <- inputDataFrame$CL*FR             # Clearance of the parent drug into the metabolite compartment

    #Calculate micro-rate constants-Parent
    inputDataFrame$k20 <- inputDataFrame$CLpop1/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #Calculate micro-rate constants-Metabolite
    inputDataFrame$kmf <- inputDataFrame$CLpop2/inputDataFrame$V2    #Rate constant for metabolite formation
    inputDataFrame$kme <- inputDataFrame$CLM/inputDataFrame$VM	 #Rate constant for metabolite elimination

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]	#Amount of parent in the absorption compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                             										#Amount of parent in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0																                                #Amount of parent in the peripheral compartment at time zero.
    inputDataFrame$AM[inputDataFrame$TIME==0] <- 0                                     												    #Amount of metabolite in the metabolite compartment at time zero.

    TwoCompFirstOrderAbsOneCompMetabCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPREDP <- inputDataFrame$A2/inputDataFrame$V2    #Concentration of parent
    inputDataFrame$IPREDM <- inputDataFrame$AM/inputDataFrame$VM    #Concentration of metabolite

    #subset unnecessary columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30,CLpop1,CLpop2,kmf,kme))

    #Return output
    inputDataFrame
}

#' Process simulations using 3-compartment first-order absorption model with 1-compartment first-order formation metabolite model.
#'
#' @description    \code{ThreeCompFirstOrderAbsOneCompMetab} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments.
#' The data frame should have the following columns: \code{ID, TIME, AMT, MDV, F1, KA, CL, V2, Q3, V3, Q4, V4, CLM, VM, FR}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KA}:\tab is the absorption rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance of parent drug\cr
#' \code{V2}:\tab is the central volume of distribution of parent drug\cr
#' \code{Q3}:\tab    is the inter-compartmental clearance of parent drug (1)\cr
#' \code{V3}:\tab is the peripheral volume of distribution of parent drug (1)\cr
#' \code{Q4}:\tab is the inter-compartmental clearance of parent drug (2)\cr
#' \code{V4}:\tab is the peripheral volume of distribution of parent drug(2)\cr
#' \code{CLM}:\tab is the central clearance of the metabolite\cr
#' \code{VM}:\tab is the    central volume of distribution of the metabolite\cr
#' \code{FR}:\tab is the    fraction of the parent drug converted into metabolite\cr
#' }
#'
#' @usage ThreeCompFirstOrderAbsOneCompMetab(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, MDV, F1, KA, CL, V2, Q3, V3, Q4, V4, CLM, VM, FR}.
#'
#' @return The function calculates the parent and metabolite amounts in the respective compartments of the pharmacokinetic model.
#' This includes the amount in the absorption \code{(A1)}, parent central (\code{A2} & individual predicted parent concentrations, \code{IPREDP})
#' parent two peripheral \code{(A3, A4)}, and metabolite central (\code{AM} & individual predicted metabolite concentrations, \code{IPREDM}) compartments.
#' The function returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), ThreeCompFirstOrderAbsOneCompMetab)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompFirstOrderAbsOneCompMetab}}, \code{\link{TwoCompFirstOrderAbsOneCompMetab}}
#' @seealso \code{\link{OneCompIVbolusOneCompMetab}}, \code{\link{TwoCompIVbolusOneCompMetab}}, \code{\link{ThreeCompIVbolusOneCompMetab}}
#' @seealso \code{\link{OneCompIVinfusionOneCompMetab}}, \code{\link{TwoCompIVinfusionOneCompMetab}}, \code{\link{ThreeCompIVinfusionOneCompMetab}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#-----------------------------------------------------------------------------------------------------------------
# First-order absorption: 3 compartment parent with 1 compartment first-order metabolite formation: RCpp function
#-----------------------------------------------------------------------------------------------------------------
ThreeCompFirstOrderAbsOneCompMetab <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT, MDV, CL, V2, Q3, V3, Q4, V4, CLM, VM, FR, KA & F1
    #Returns a dataframe with populated columns for A1, A2, A3, A4, AM, IPREDP, IPREDM

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- k24 <- k42 <- k40 <- kmf <- kme <- NULL
    CLpop1 <- CLpop2 <- NULL

    #Calculate micro-rate constants
    FR = inputDataFrame$FR[1]
    inputDataFrame$CLpop1 <- inputDataFrame$CL*(1-FR)     # Clearance of the parent drug to outside the body
    inputDataFrame$CLpop2 <- inputDataFrame$CL*FR             # Clearance of the parent drug into the metabolite compartment

    #Calculate micro-rate constants-Parent
    inputDataFrame$k20 <- inputDataFrame$CLpop1/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q3/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q3/inputDataFrame$V3
    inputDataFrame$k30 <- 0
    inputDataFrame$k24 <- inputDataFrame$Q4/inputDataFrame$V2
    inputDataFrame$k42 <- inputDataFrame$Q4/inputDataFrame$V4
    inputDataFrame$k40 <- 0

    #Calculate micro-rate constants-Metabolite
    inputDataFrame$kmf <- inputDataFrame$CLpop2/inputDataFrame$V2    #Rate constant for metabolite formation
    inputDataFrame$kme <- inputDataFrame$CLM/inputDataFrame$VM	 #Rate constant for metabolite elimination

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]        # Amount in the absorption compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount of parent in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount of parent in the 1st peripheral compartment at time zero.
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount of parent in the 2nd peripheral compartment at time zero.
    inputDataFrame$AM[inputDataFrame$TIME==0] <- 0                                     # Amount of metabolite in the metabolite compartment at time zero.

    ThreeCompFirstOrderAbsOneCompMetabCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPREDP <- inputDataFrame$A2/inputDataFrame$V2    #Concentration of parent
    inputDataFrame$IPREDM <- inputDataFrame$AM/inputDataFrame$VM    #Concentration of metabolite

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30,k24,k42,k40,CLpop1,CLpop2,kmf,kme))

    #Return output
    inputDataFrame
}
