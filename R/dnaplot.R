#' Densityplot for DNA Profiles
#' 
#' A wrapper for [flowViz::densityplot()] with defaults for
#' DNA profiles. 
#' 
#' @md
#' @param ... Arguments to be handed to [flowViz::densityplot()]
#' @param data `flowFrame` or `flowSet` with data
#' @param adjust Numeric factor to adjust the bandwith in the kernel
#'   density function [stats::density()], decreased from 1 to 0.2
#' @param n The number of equally spaced points at which the density is to
#'   be estimated, increased from 50 to 512
#' @param plot If `TRUE`, lattice object will be plotted as well as 
#'   being returned as an `invisible` object
#'
#' @details
#' This is a convenience function with parameters more suited for
#' typical DNA profiles. Specifically, the optional arguments to
#' [stats::density()], which are provided in `darg`, are changed from
#' `n = 50` and `adjust = 1` to `n = 512` and `adjust = 0.2`.
#'
#' The call effectively executes the following:
#' ```
#'   densityplot(..., dargs = list(adjust = adjust, n = n))
#' ```
#' 
#' @return
#' Optionally plot the data as a side-effect and invisibly
#' return the lattice plot.
#' 
#' @examples 
#'  fs <- readSet(system.file("extdata", "RPE_synch/", package = "flowExtra"))
#'  dnaplot(~ FL2.A, fs) # unadjusted values
#'
#' @import flowViz
#'
#' @export
#'
dnaplot <- function(..., adjust = 0.2, n = 512, plot = TRUE)
{
	dots <- list(...)
	darg <- list(adjust = adjust, n = n)
	argList <- c(dots, darg = list(darg))
	obj <- do.call(flowViz::densityplot, args = argList)
	if (plot == TRUE)
		plot(obj)
	invisible(obj)
}
