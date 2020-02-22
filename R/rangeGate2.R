#' Improved search for separation between positive and negative populations
#'
#' This function replaces [flowStats::rangeGate()] which never seemed to 
#' work for me. It does not perform all the functions of the original
#' but will return a mores sensible `rectangleGate` object more
#' often than [flowStats::rangeGate()].
#'
#' @md
#' @param fs A [`flowCore::flowset`] or [`flowCore::flowFrame`] object
#' @param stain A `character` scalar (also known as `chan`) of the flow
#'  parameter to analyze
#' @param from,to Starting and ending values to search for the break
#' @param cutoff Peaks and valleys smaller than this fraction of the range
#'   will be ignored, default value of 0.05
#' @param adjust Adjustment to bandwidth used by [stats::density()]
#' @param half.range If `TRUE` and a single population is found, the
#'   normal distribution will be fit using the left half of the population
#' @param sd In the case of a single population, this is the multiplier
#'  for the standard deviation to estimate the upper edge of the 
#'  negative population, default value of 2.5
#' @param positive If `TRUE`, gate the positive population
#' @param method Either `"minimum"` or `"left"` to choose the break as the
#'  minimum between two peaks or to treat the left-most population
#'  as a Gaussian population of the negative population
#' @param plot If `TRUE`, show a plot of the data and breakpoint
#' @param legend If `TRUE`, include a legend in the plot (if `plot == TRUE`)
#' @param filterId The name assigned to the filter (`character` value)
#' @param ... Other arguments to the original function are accepted
#'  such as `borderQuant`, `absolute`, and `refLine` but are discarded 
#' 
#' @details
#'
#' This is a replacement function for [flowStats::rangeGate] that may be
#' more robust than the original function. 
#' 
#' This function seeks the most likely breakpoint between two
#' populations to differentiate between the positive and negative
#' population. The
#' "negative" population is assumed to have a lower values and would be
#' located on the left of a histogram.
#' 
#' The default method (`method = "minimum"`) attempts to identify the
#' lowest minimum between the regions of greatest density or peaks on
#' a densityplot. The search can be limited to values between `from` and
#' `to`. If more than four high density regions exist *or* if only a
#' single high density region is detected, the breakpoint will be
#' determined from the properties of the leftmost population where this
#' population is assumed to follow a normal distribution. This approach
#' also can be specified directly by `method = "left"` as described
#' below.
#'
#' With `method = "left"`, the maximum of the left most (minimum)
#' population is determined from a kernel density estimate. If
#' `half.range = TRUE`, the *left* half of the population will be used to
#' determine the Gaussian distribution in order to estimate the standard
#' deviation of the population. With `method = "left"`, the value returned
#' is the position of the peak + `sd` times the standard deviation of the
#' distribution.
#' 
#' A diagnostic plot will be generated with base graphics if `plot = TRUE`.
#' In this case, the default argument `legend = TRUE` will add
#' an informative legend to the plot showing the data and breakpoint. The
#' plotting option can be useful to iteratively adjust the search parameters
#' such as `adjust`, `from`, `to` and `sd`. 
#'
#' @seealso [flowStats::rangeGate()]
#'
#' @return
#' An object of class [`flowCore::rectangleGate`]
#' 
#' @import
#' flowCore
#'
#' @importFrom MASS fitdistr
#' 
#' @examples
#' # Read and clean up synchronized cell data
#'   fs <- readSet(system.file("extdata", "synch", package = "flowExtra"))
#'   fs <- Subset(fs, linearGate(fs, "FL2.A", "FL2.H"))
#' # breakpoint by default (method = "minimum")
#'   rangeGate2(fs[[8]], "FL2.A", plot = TRUE)
#' # breakpoint by left population (method = "left")
#'   rangeGate2(fs[[8]], "FL2.A", method = "left", plot = TRUE)
#' # puzzling choice with original rangeGate() function
#'   flowStats::rangeGate(fs[[8]], "FL2.A", plot = TRUE)
#' # example of multiple possible breakpoints
#'   rangeGate2(fs[[2]], "FL2.A", plot = TRUE)
#' # adjust by limiting search range with 'to'
#'   rangeGate2(fs[[2]], "FL2.A", plot = TRUE, to = 275)
#'
#' @export
#' 
rangeGate2 <- function(fs, stain, from = NULL, to = NULL, cutoff = 0.05,
	adjust = 1, half.range = TRUE, sd = 2.5, positive = TRUE,
	method = c("minimum", "left"), plot = FALSE, legend = TRUE,
	filterId = "defaultRectangleGate", ...)
{
# check and process arguments
	flowCore:::checkClass(fs, c("flowFrame", "flowSet"))
	if (is(fs, "flowFrame"))
		fs <- as(fs, "flowSet")
	flowCore:::checkClass(stain, "character", 1)
	if (!stain %in% colnames(fs))
		stop("'", stain, "' is not a valid parameter in this flowSet")
	ilim <- range(fs[[1]], type = "instrument")[, stain]
	dlim <- fsApply(fs[, stain], range, type = "data")
	dlim <- as.matrix(do.call(cbind, dlim))
	dlim <- apply(dlim, 1, range)[c(1, 4)] # min of min and max of max
	from <- if(is.null(from)) dlim[1] else max(from[1], dlim[1])
	to <- if(is.null(to)) dlim[2] else min(to[1], dlim[2])
	flowCore:::checkClass(cutoff, "numeric", 1)
	if (cutoff < 0 | cutoff > 1)
		stop("'cutoff' must be between 0 and 1")
	flowCore:::checkClass(adjust, "numeric", 1)
	flowCore:::checkClass(half.range, "logical", 1)
	flowCore:::checkClass(sd, "numeric", 1)
	flowCore:::checkClass(positive, "logical", 1)
	flowCore:::checkClass(plot, "logical", 1)
	flowCore:::checkClass(legend, "logical", 1)
	method <- match.arg(method)

# extract data between 'from' and 'to'
	dat <- fsApply(fs, exprs)[, stain]
	x <- dat[dat >= from & dat <= to]

# find peaks and valleys
	d <- stats::density(x, adjust = adjust)
	zc <- flowExtra:::zero.cross(d$x, d$y)

# remove artifactual adjacent peaks or valleys
	keep <- c(diff(zc$sign) != 0, TRUE)	# discard first duplicates
	zc$x <- zc$x[keep]
	zc$y <- zc$y[keep]
	zc$sign <- zc$sign[keep]

# discard those peaks below the adjusted cutoff
	adjcut <- cutoff * diff(range(d$y)) + min(d$y)
	sel <- which(d$x %in% zc$x[zc$sign == 1])
	xp <- d$x[sel]
	yp <- zc$model$y[zc$model$x %in% xp]
	sel <- yp > adjcut
	xp <- xp[sel]
	yp <- yp[sel]
	
# if 0 peaks, fail
	if (length(xp) == 0)
		stop("unable to identify any peaks above cutoff")

# if 2 to 4 peaks, use lowest minimum between peaks
	if (method == "minimum" && length(xp) >= 2 && length(xp) <= 4) {
		sel <- zc$x > min(xp) & zc$x < max(xp) & zc$sign == -1 # select min
		yp.loc <- zc$y[sel] # all minimum values between xp[1] and xp[2]
		loc <- zc$x[sel][which.min(yp.loc)] # x-value from lowest valley
	}

# else select the left-most peak from kernel density estimate
	else {
		xmid <- min(xp)
		if (half.range == TRUE) {
			xl <- x[x <= xmid]
			xx <- c(xl, 2*xmid  - xl)
		}
		else
			xx <- x
		if (!length(xx))
			stop("unable to determine left-most peak")
		fit <- MASS::fitdistr(xx, "normal")
		loc <- fit$est[1] + sd * fit$est[2]
	}

# warn about too many peaks
	if (length(xp) > 3)
		warning(length(xp), " peaks were found with adjust = ", adjust,
			". Consider increasing this value.")

# show plot with base graphics
	if (plot == TRUE) {
	# need flat caps for fat line
		opar <- par(lend = 3)
	# assemble legend components
		leg.txt <- c(sprintf("breakpoint (%0.4g)", loc),
			ifelse(length(xp) > 1, "peaks", "peak"), "gated region",
			"search region")
		if (method == "left" | length(xp) > 4 | length(xp) == 1)
			leg.title <- paste('"Left" method, sd =', round(sd, 1))
		else
			leg.title <- '"Minimum" method'
	# assemble plot
		main.txt <- paste("Breakpoint for ", stain)
		d <- density(dat, adjust = adjust)
	# polygon to highlight gated population
		xv <- if (positive) d$x[d$x > loc] else d$x[d$x < loc]
		yv <- if (positive) d$y[d$x > loc] else d$y[d$x < loc]
	# plot empty plot first, then polygon, then overlay with lines
		plot(d, main = main.txt, xlim = ilim, type = "n")
		polygon(c(min(xv), xv, max(xv)), c(0, yv, 0),
			col = "gray90", border = NA)
		lines(d)
		rug(xp, col = "blue")
		abline(v = loc, col = "red", lwd = 2)
		ymid <- max(d$y)/4
		arrows(min(x), ymid, max(x), ymid, length = 0.2, angle = 90,
			lty = 2, code = 3)
	# add legend if requested
		if (legend == TRUE) 
			legend("topright", legend = leg.txt, title = leg.title,
				inset = 0.05, bg = "white", box.col = "white",
				col = c("red", "blue", "gray90", "black"),
				lty = c(1, 1, 1, 2), lwd = c(2, 1, 10, 1))
		par(opar)
	}

# return rectangleGate with channel limits defined by instrument
	if (positive)
		bounds <- list(c(loc, ilim[2]))
	else
		bounds <- list(c(ilim[1], loc))
	names(bounds) <- stain
	rectangleGate(bounds, filterId = filterId)
}
