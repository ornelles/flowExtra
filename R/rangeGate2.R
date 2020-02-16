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
#' @param method One of `"min"` or `"left"` to choose the break as the
#'  minimum between two peaks or to treat the left-most population
#'  as a Gaussian population of negative values
#' @param plot If `TRUE`, show a plot of the data and breakpoint
#' @param filterId The name assigned to the filter (`character` value)
#' @param ... Other arguments to the original function are accepted
#'  such as `borderQuant`, `absolute`, and `refLine` but are ignored 
#' 
#' @details
#' 
#' If `method = "minimum"`, the breakpoint is sought first as the minimum
#' between two high density regions between `from` and `to`. If more than
#' two high density regions exist or if only a single high
#' density region be detected, the breakpoint will be 
#' determined from the properties of the leftmost population 
#' where this population is assumed to follows a normal distribution.
#' This approach can be specified directly by `method = "left"`.
#' The maximum of this population is determined from a kernel density
#' estimate. If `half.range = TRUE`, the **left** half of the population
#' will be fit to the presumed Gaussian distribution. With a single
#' peak, the value returned is the position of
#' the peak + `sd` times the standard deviation of the distribution.
#'
#' @seealso [flowStats::rangeGate()]
#'
#' @return
#' An object of class ['flowCore::rectangleGate']
#' 
#' @import
#' flowCore
#'
#' @importFrom MASS fitdistr
#' 
#' @examples
#'
#' @export
#' 
rangeGate2 <- function(fs, stain, from = NULL, to = NULL, cutoff = 0.05,
	adjust = 1, half.range = TRUE, sd = 2.5, positive = TRUE,
	method = c("minimum", "left"), plot = FALSE,
	filterId = "defaultRectangleGate", ...)
{
# check and process arguments
	flowCore:::checkClass(fs, c("flowFrame", "flowSet"))
	if (is(fs, "flowFrame"))
		fs <- as(fs, "flowSet")
	flowCore:::checkClass(stain, "character", 1)
	if (!stain %in% colnames(fs))
		stop("'", stain, "' is not a valid parameter in this flowSet")
	lims <- fsApply(ws[, stain], range, type = "data")
	lims <- as.matrix(do.call(cbind, lims))
	lims <- apply(lims, 1, range)[c(1, 4)] # min of min and max of max
	from <- if(is.null(from)) lims[1] else max(from[1], lims[1])
	to <- if(is.null(to)) lims[2] else min(to[1], lims[2])
	flowCore:::checkClass(cutoff, "numeric", 1)
	if (cutoff < 0 | cutoff > 1)
		stop("'cutoff' must be between 0 and 1")
	flowCore:::checkClass(adjust, "numeric", 1)
	flowCore:::checkClass(half.range, "logical", 1)
	flowCore:::checkClass(sd, "numeric", 1)
	flowCore:::checkClass(positive, "logical", 1)
	flowCore:::checkClass(plot, "logical", 1)
	method <- match.arg(method)

# extract and adjust data between 'from' and 'to'
	dat <- fsApply(fs, exprs)[, stain]
	x <- dat[dat >= from & dat <= to]

# find peaks
	d <- stats::density(x, adjust = adjust)
	zc <- flowExtra:::zero.cross(d$x, d$y)

# remove adjacent peaks or valley artifacts
	keep <- c(diff(zc$sign) != 0, TRUE)	# discard first duplicates
	zc$x <- zc$x[keep]
	zc$y <- zc$y[keep]
	zc$sign <- zc$sign[keep]

# extract peaks and discard those below adjusted cutoff
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

# if 2 peaks, use lowest minimum between peaks
	if (method == "minimum" && length(xp) == 2) {
		sel <- zc$x > xp[1] & zc$x < xp[2] & zc$sign == -1 # selects minimums
		yp.loc <- zc$y[sel] # all minimum values between xp[1] and xp[2]
		loc <- zc$x[sel][which.min(yp.loc)] # x-value from lowest peak
	}

# select left-most peak from kernel density estimate
	else {
		xmid <- min(xp)
		if (half.range == TRUE) {
			xl <- x[x <= xmid]
			xx <- c(xl, 2*xmid  - xl)
		}
		else
			xx <- x
		fit <- MASS::fitdistr(xx, "normal")
		loc <- fit$est[1] + sd * fit$est[2]
	}

# warn about too many peaks
	if (length(xp) > 3)
		warning(length(xp), " peaks were found with adjust = ", adjust, ".\n",
			"Consider increasing the value of 'adjust'.")

# show plot?
	if (plot == TRUE) {
		opar <- par(lend = 3)
		main.txt <- paste("Breakpoint for ", stain)
		leg.txt <- c(sprintf("breakpoint (%0.4g)", loc), "gated data", "included data")
		if (method == "left" | length(xp) != 2)
			leg.title <- paste('"Left" method, sd =', round(sd, 1))
		else
			leg.title <- '"Minimum" method'
		d <- density(dat, adjust = adjust)
		plot(d, main = main.txt, xlim = lims)
		xv <- if (positive) d$x[d$x > loc] else d$x[d$x < loc]
		yv <- if (positive) d$y[d$x > loc] else d$y[d$x < loc]
		polygon(c(min(xv), xv, max(xv)), c(0, yv, 0), col = "gray90")
		abline(v = loc, col = 2, lwd = 2)
		ymid <- max(d$y)/3
		arrows(min(x), ymid, max(x), ymid, length = 0.2, angle = 90,
			lty = 3, code = 3)
		legend("topright", legend = leg.txt, title = leg.title, inset = 0.05,
			bg = "white", box.col = "white", col = c("red", "gray90", "black"),
			lty = c(1, 1, 3), lwd = c(2, 10, 1))
		par(opar)
	}

# return rectangleGate with channel limits
	if (positive)
		bounds <- list(c(loc, lims[2]))
	else
		bounds <- list(c(lims[1], loc))
	names(bounds) <- stain
	rectangleGate(bounds, filterId = filterId)
}
