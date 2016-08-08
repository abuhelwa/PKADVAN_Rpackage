#' PKADVAN: Pharmacokinetic Simulations Using the ADVAN-style Analytical Solutions of Common Pharmacokinetic Models.
#'
#' @description The \code{PKADVAN} is an R-package for simulating pharmacokinetic data from various population pharmacokinetic models.
#' The \code{PKADVAN} package presents the ADVAN-style analytical solutions functions for 26 different pharmacokinetic models, including the basic models published
#' by \emph{Abuhelwa et al.} (1), for the 1, 2, and 3 compartments IV bolus, infusion, and first-order absoprtion models.
#' The ADVAN-style analytical solutions functions simulate the time-course of drug amounts in the respective compartments of a pharmacokinetic system and
#' they \emph{ADVANCE} the solution of the model from one time point to the next, allowing for any dose or time-changing covariates to be accounted for.
#'
#' The \code{PKADVAN} functions presented herein were written using an integrated R/C++ functions for significantly faster computational speed as compared to the R-only coded functions.
#'
#' A main application for the \code{PKADVAN} package is for simulation from stochastic population pharmacokinetic
#' models coded in R and incorporated into reactive Shiny Web applications. The speed of the integrated R/C++
#' analytical solutions is an important factor as they allow for simulating larger populations without compromising the reactivity benefits of the Shiny application.
#'
#' All the \code{PKADVAN} functions were successfully validated against the commercially available population pharmacokinetic modelling software; NONMEM (2).
#'
#' @details \tabular{ll}{
#' Package:\tab \code{PKADVAN}\cr
#' Type:\tab Package\cr
#' Version:\tab 0.1.0\cr
#' Date:\tab 05-08-2016\cr
#' License:\tab GPL (>= 2)\cr
#' }
#'
#' The \code{PKADVAN} functions of a pharmacokinetic model accept a NONMEM-style data frame and calculate drug
#' amounts in the respective compartments and the individual predicted concentrations in the central compartment of a pharmacokinetic system.
#'
#' The data frame should have the following columns: \code{ID, TIME, AMT}, in addition to the individual pharmacokinetic parameters of the respective pharmacokinetic model
#' (e.g. \code{CL, V, Q}). Please see help pages for more details!
#'
#' @section PKADVAN Functions:
#' The PKADVAN functions are capable of simulating arbitrary dosing regemins and can account for time-changing covariate
#' structures; however, covariate effects on respective parameters must be calculated prior porocessing
#' simulations. Please see examples for more details!
#'
#' @seealso The following links provides more details on how to use the \code{PKADVAN} functions for processing
#' simulations from various pharmacokinetic models.
#'
#' @seealso Basic Pharmacokinetic Models:
#' @seealso \code{\link{OneCompIVbolus}}, \code{\link{TwoCompIVbolus}}, \code{\link{ThreeCompIVbolus}}
#' @seealso \code{\link{OneCompIVinfusion}}, \code{\link{TwoCompIVinfusion}}, \code{\link{ThreeCompIVinfusion}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#'
#' @seealso Transit Absorption Models:
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompTwoTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompThreeTransit}}, \code{\link{TwoCompFourTransit}}
#'
#' @seealso First-order Formation Metabolite Models:
#' @seealso \code{\link{OneCompIVbolusOneCompMetab}}, \code{\link{TwoCompIVbolusOneCompMetab}}, \code{\link{ThreeCompIVbolusOneCompMetab}}
#' @seealso \code{\link{OneCompIVinfusionOneCompMetab}}, \code{\link{TwoCompIVinfusionOneCompMetab}}, \code{\link{ThreeCompIVinfusionOneCompMetab}}
#' @seealso \code{\link{OneCompFirstOrderAbsOneCompMetab}}, \code{\link{TwoCompFirstOrderAbsOneCompMetab}}, \code{\link{ThreeCompFirstOrderAbsOneCompMetab}}
#'
#'
#' @references (1) Abuhelwa, A.Y., D.J. Foster, and R.N. Upton, \emph{ADVAN-style analytical solutions for common pharmacokinetic models}. Journal of pharmacological and toxicological methods, 2015. \bold{73}: p. 42-48.
#'
#' (2) Beal, S., et al., \emph{NONMEM User's Guides, Part V}. (1989-2009), Icon Development Solutions, Ellicott City, MD, USA, 2009.

#' @name PKADVAN-package
#' @author Ahmad Abuhelwa, David Foster, Richard Upton\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}

NULL


