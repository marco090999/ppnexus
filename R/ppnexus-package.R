#' ppnexus: Simulation, Modelling and Diagnostics for Multidimensional
#' Point Processes
#'
#' Tools for simulation, estimation, diagnostics, and exploratory analysis of
#' multidimensional spatial and spatio-temporal point processes, including
#' Poisson, marked Poisson, functional-covariate, log-Gaussian Cox, and ETAS
#' models.
#'
#' @keywords internal
#' @useDynLib ppnexus, .registration = TRUE, .fixes = "C_"
#' @importFrom graphics plot
#' @importFrom stats predict setNames simulate
#' @importFrom utils tail
"_PACKAGE"
