#' Densityplot for DNA Profiles
#' 
#' This is a wrapper for [flowViz::densityplot()] with defaults for
#' DNA profiles. 
#' 
#' @md
#' @param ... Arguments to be handed to [flowViz::densityplot()]
#' @param data `flowFrame` or `flowSet` with data
#' @param adjust Numeric factor to adjust the bandwith in the kernel
#'   density function [stats::density()], decreased from 1 to 0.2
#' @param n The number of equally spaced points at which the density is to
#'   be estimated, increased from 50 to 512
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
#' @return Invisibly return the lattice object.
#' 
#' @import flowViz
#'
#' @export
#'
dnaplot <- function(..., adjust = 0.2, n = 512)
{
	dots <- list(...)
	darg <- list(adjust = adjust, n = n)
	argList <- c(dots, darg = list(darg))
	obj <- do.call(flowViz::densityplot, args = argList)
	plot(obj)
	invisible(obj)
}
