#' Fit Gate for Single Cells
#' 
#' Fit a rectangular polygon gate for single cells in FL2.A x FL2.H
#' 
#' @details
#' This function fits values in \code{xchan} and \code{ychan} to a robust
#' linear regression with \code{\link[MASS]{lqs}}. Note that default names
#' for the x & y channels are "FL2.A" and "FL2.H" and \strong{not} "FL2-A"
#' and "FL2-H". If \code{zero.intercept} is \code{TRUE}, the fit is forced
#' through origin. The regression will be limited to values within the
#' 2.5 to 97.5 percentile if \code{xrange} and \code{yrange} is \code{NULL}.
#' A different range can be specified with these arguments.
#' 
#' The argument \code{x} can be a matrix from \code{\link[Biobase]{exprs}},
#' a \code{flowFrame}, or a \code{flowSet}. A recantgular polygon gate (or
#' list of gates) with a width and height of \code{dx} and \code{dy} will
#' be returned. If \code{limit.gate} is \code{TRUE}, the returned gate will
#' also be limited by the values in \code{xrange} and \code{yrange}.
#' Otherwise, the gate will be limited by the range of values in each
#' flowFrame up to limits set by \code{cutoff} and the maximum range of
#' values in the x channel (FL2.A) limited to \code{cutoff * maxRange} and
#' \code{(1 - cutoff) * maxRange}.
#' 
#' @param x Data to fit as either a matrix from \code{\link[Biobase]{exprs}},
#'   a \code{flowFrame}, or a \code{flowSet}
#' @param xchan,ychan Character strings to define single cells, default
#'   values of \code{"FL2.A"} and \code{"FL2.H"}
#' @param zero.intercept If \code{TRUE} (default), force regression through
#'   origin
#' @param dx,dy Width and height of rectangular gate, default value of 0.06.
#'   If only \code{dx} is specified, \code{dy} will be assigned the same value 
#' @param xrange,yrange If \code{NULL}, regression will be limited to
#'   values within the 2.5 - 97.5 percentile, otherwise in the specified
#'   range of values
#' @param cutoff Numeric value between 0 and 1 specifying the fraction of
#'   low and high values to ignore in the regression, default value of 0.02
#' @param limit.gate Logical value (default of \code{FALSE}) to limit the
#'   returned gate by the values in \code{xrange} and \code{yrange}
#' @param filterId Character string with default value of \code{"singlets"}
#' @param groupFilterId Character string with default value of \code{
#'   "linearGateList"}
#' @return
#' A rectangular \code{polygonGate} or list of these gates.
#' 
#' @import flowCore
#' @importFrom MASS lqs
#'
#' @export
#'

linearGate <- function(x, xchan = "FL2.A", ychan = "FL2.H",
		zero.intercept = TRUE, dx = 0.06, dy = dx, xrange = NULL,
		yrange = NULL, cutoff = 0.02, limit.gate = FALSE,
		filterId = "singlets", groupFilterId = "linearGateList")
{
# argument checks
	if(!(is(x, "flowFrame") | is(x, "flowSet") | is(x, "matrix")))
		stop("require either a flowFrame, flowSet or matrix")
	nms <- colnames(x)
	if (!xchan %in% nms)
		stop('"', xchan, '" not found in ', deparse(substitute(x)))
	if (!ychan %in% nms)
		stop('"', ychan, '" not found in ', deparse(substitute(x)))

# set maxRange for x values
	if (is(x, "flowFrame"))
		maxRange <- range(x)[2, xchan]
	else if (is(x, "flowSet"))
		maxRange <- range(x[[1]])[2, xchan]
	else if (is(x, "matrix"))
		maxRange <- max(x, na.rm = TRUE)
	else
		stop("Whoa! Should not have arrived here...")

# working function expects matrix with columns named xchan, ychan
	.linear.Gate <- function(mat, xchan, ychan, dx, dy, xrange, yrange,
			limit.gate, cutoff, maxRange, filterId)
	{
		dat <- data.frame(mat[, c(xchan, ychan)])
		names(dat) <- c("x", "y")
	
		if (is.null(xrange))
			xrange <- quantile(dat$x, c(0.025, 0.975), na.rm = TRUE)
		if (is.null(yrange))
			yrange <- quantile(dat$y, c(0.025, 0.975), na.rm = TRUE)
		dat <- subset(dat, x > min(xrange) & x < max(xrange))
		dat <- subset(dat, y > min(yrange) & y < max(yrange))

		if (zero.intercept)
			fm <- try(MASS::lqs(y ~ x + 0, data = dat), silent = TRUE)
		else
			fm <- try(MASS::lqs(y ~ x, data = dat), silent = TRUE)

		if (class(fm) != "try-error") { # continue processing
			if (limit.gate)
				xp <- xrange
			else {
				xp <- range(mat[, xchan], na.rm = TRUE)
				xp[1] <- max(xp[1], round(cutoff * maxRange))
				xp[2] <- min(xp[2], round((1 - cutoff) * maxRange))
			}
			yp <- predict(fm, newdata = data.frame(x = xp))
	
			xinc <- 0.5 * dx * diff(range(dat$x, na.rm = TRUE))
			yinc <- 0.5 * dy * diff(range(dat$y, na.rm = TRUE))
	
			x.vertices <- c(xp[1] - xinc, xp[1] + xinc, xp[2] + xinc, xp[2] - xinc)
			y.vertices <- c(yp[1] + yinc, yp[1] - yinc, yp[2] - yinc, yp[2] + yinc)
			filterIdValue <- filterId
		}
		else { # lqs failed, return zero-width gate
			x.vertices <- c(xrange[1], xrange[1], xrange[2], xrange[2])
			y.vertices <- c(yrange[1], yrange[1], yrange[2], yrange[2])
			filerIdValue <- "unfiltered"
		}
		bd <- cbind(x.vertices, y.vertices)
		colnames(bd) <- c(xchan, ychan)

		return(polygonGate(.gate = bd, filterId = filterIdValue))
	}

# dispatch function on nature of argument 'x'
	if (class(x) == "matrix")
		ret <- .linear.Gate(x, xchan, ychan, dx, dy, xrange, yrange,
						limit.gate, cutoff, maxRange, filterId)
	else if (class(x) == "flowFrame")
		ret <- .linear.Gate(exprs(x), xchan, ychan, dx, dy, xrange, yrange,
						limit.gate, cutoff, maxRange, filterId)
	else if (class(x) == "flowSet") {
		ret <- fsApply(x, .linear.Gate, xchan, ychan, dx, dy, xrange, yrange,
						limit.gate, cutoff, maxRange, filterId, use.exprs = TRUE)
		ret <- filterList(ret, filterId = groupFilterId)
	}
	else
		stop("cannot process argument of class \"", class(x), "\"")

	return(ret)
}
