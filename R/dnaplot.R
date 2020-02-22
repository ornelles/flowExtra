#' Densityplot for DNA Profiles
#' 
#' A wrapper for [flowViz::densityplot()] with defaults for
#' DNA profiles. 
#' 
#' @md
#' @param ... Arguments to be handed to [flowViz::densityplot()] *or*
#'   a `flowFrame` or `flowSet`
#' @param darg List of arguments passed to the kernel density function
#'   [stats::density()] to increase the resolution of the plot. In
#'   addition to `bw`, these can include `adjust`, `n`. The option
#'   `na.rm = TRUE` will be included in all cases
#' @param chan An option character scalar identifying the data to be
#'   plotted if `...` is a `flowFrame` or `flowSet`
#' @param plot If `TRUE`, lattice object will be plotted as well as 
#'   being returned as an `invisible` object
#'
#' @details
#' This is a convenience function with parameters more suited for
#' typical DNA profiles by using the bandwidth algorithm of Sheather
#' & Jones in `bw.SJ` with the adjustment in `adjust` and using more
#' points (512 rather than 50) to determine the density. These additional
#' parameters are passed to [stats::density()], which in `darg`.
#'
#' The call effectively executes the following:
#' ```
#'   densityplot(..., dargs = dargs)
#' ```
#' ... **except** for the *lazy* case of calling it with first argument a
#' `flowSet` or `flowFrame`. In this case the function guesses
#' the most likely `chan` parameter (in order) from `FL2.A, FL3.A, FL2.H`
#' or `FL3.H`. Alternatively, the channel name can be passed in `chan` and
#' the call now becomes:
#' ```
#'   densityplot(~ chan, ..., dargs = dargs)
#' ```
#' 
#' @return
#' Optionally plot the data as a side-effect and invisibly
#' return the lattice plot.
#' 
#' @examples 
#' # load and filter example data with 16 flowFrames
#'   fs <- readSet(system.file("extdata", "synch/", package = "flowExtra"))
#'   fs <- Subset(fs, boundaryFilter("FL2.A"))
#' # default density parameters from flowViz 
#'   dnaplot(~ FL2.A, fs, xlim = c(50, 500), darg = NULL,
#'     main = "Default parameters for densityplot")
#' # revised density parameters 
#'   dnaplot(~ FL2.A, fs, xlim = c(50, 500),
#'     main = "Adjusted parameters for 'dnaplot'")
#' # lazy use if "FL2.A" is present
#'   dnaplot(fs)
#'
#' @importFrom flowViz densityplot
#'
#' @export
#'
dnaplot <- function(..., darg = list(bw = "SJ", n = 512, adjust = 1),
	chan = NULL, plot = TRUE)
{
# extract first argument
	dots <- list(...)

# special case of lazily treating this like an S4 function...
	if (is(dots[[1]], "flowFrame") | is(dots[[1]], "flowSet")) {
		x <- dots[[1]]
		if (is.null(chan)) {
			choices <- c("FL2.A", "FL3.A", "FL2.H", "FL3.H")
			sel <- pmatch(colnames(x), choices)
			if (all(is.na(sel)))
				stop("need one of: ", paste(choices, collapse = ", "), " in data")
			else	
				chan <- choices[min(sel, na.rm = TRUE)]
		}
		else if (!chan %in% colnames(x))
			stop ('"', deparse(substitute(chan)), '" not found')
		form <- as.formula(paste("~", chan))
		dots <- c(form, dots)
	}
# ensure that 'na.rm = TRUE' is in 'darg'
	if (!"na.rm" %in% names(darg))
		darg <- c(darg, na.rm = TRUE)
# make call
	argList <- c(dots, darg = list(darg))
	obj <- do.call(flowViz::densityplot, args = argList)
# plot object unless told not to...
	if (plot == TRUE)
		plot(obj)
	invisible(obj)
}
