#' Densityplot for DNA Profiles
#' 
#' A wrapper for [flowViz::densityplot()] with defaults for
#' DNA profiles. 
#' 
#' @md
#' @param ... Arguments to be handed to [flowViz::densityplot()] *or*
#'   a `flowFrame` or `flowSet`
#' @param chan An option character scalar identifying the data to be
#'   plotted if `...` is a `flowFrame` or `flowSet`
#' @param bw,n,adjust Parameters passed to the kernel density function
#'   [stats::density()] to increase the resolution of the plot with
#'   default values of `bw = "SJ"`, `n = 512` and `adjust = 1`. The option
#'   `na.rm = TRUE` will be included in all cases
#'
#' @details
#' This is a convenience function with parameters more suited for
#' typical DNA profiles by using the bandwidth algorithm of Sheather
#' & Jones in `bw.SJ` with the adjustment in `adjust` and using more
#' points (512 rather than 50) to determine the density. 
#'
#' The call effectively executes the following:
#' ```
#' densityplot(..., dargs = list(bw = bw, n = n, adjust = adjust, na.rm = TRUE))
#' ```
#' **Except** for the *lazy* case of calling it with first argument a
#' `flowSet` or `flowFrame`. In this case the function guesses
#' the most likely `chan` parameter (in order) from `FL2.A, FL3.A, FL2.H`
#' or `FL3.H`. Alternatively, the channel name can be passed in `chan` and
#' the call now becomes:
#' ```
#'   densityplot(~ chan, ..., dargs = list(bw = bw, n = n, adjust = adjust, na.rm = TRUE))
#' ```
#' 
#' @return
#' Plot the data and invisibly return the lattice object.
#' 
#' @examples 
#' # load and filter example data with 16 flowFrames
#'   fs <- readSet(system.file("extdata", "synch/", package = "flowExtra"))
#'   fs <- Subset(fs, boundaryFilter("FL2.A"))
#' # default density parameters from flowViz 
#'   densityplot(~ FL2.A, fs, xlim = c(50, 500),
#'     main = "Default parameters for densityplot")
#' # revised density parameters 
#'   dnaplot(~ FL2.A, fs, xlim = c(50, 500),
#'     main = "Refined parameters for 'dnaplot'")
#' # lazy use if "FL2.A" is present
#'   dnaplot(fs)
#'
#' @importFrom flowViz densityplot
#'
#' @export
#'
dnaplot <- function(..., chan = NULL, bw = "SJ", n = 512, adjust = 1)
{
# extract arguments that are not passed to density()
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

# extract three arguments for density()
	darg <- list(bw = bw, n = n, adjust = adjust, na.rm = TRUE)

# make call
	argList <- c(dots, darg = list(darg))
	obj <- do.call(flowViz::densityplot, args = argList)

# plot and return invisible object
	plot(obj)
	invisible(obj)
}
