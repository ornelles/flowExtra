#' Create Gate for Single Cells
#' 
#' Fit a narrow rectangular gate to select single cells in a 2-D
#' plot of area x height
#' 
#' @details#' This function creates a data-driven gate to select single cells from
#' \code{area} and \code{height} values collected on a linear scale. A 
#' robust linear regression is fit to the data with \code{\link[MASS]{lqs)}.
#' This regression is used to create a narrow rectangular gate 
#' centered about the single cells. If \code{zero.intercept = TRUE}, the 
#' regression is forced through origin. By default, the regression is 
#' applied to values within the 2nd to 98th percentile of data in each
#' \code{flowFrame}. Different percentiles can be specified with the 
#' arguments \code{xquant} and \code{yquant}.
#' 
#' The resulting \code{polygonGate(s)} will be limited by a cutoff 
#' factor determined from \code{gate.limit} and the range of data 
#' in each channel (\code{maxRange}). If \code{gate.limit} is a numeric 
#' vector of values less or equal to 1, the gate will be limited to
#' \code{maxRange*min(gate.limit)} to \code{maxRange*max(gate.limit)}.
#' If \code{gate.limit} is single numeric value less than 1, the limits
#' will be \code{maxRange * gate.limit} to \code{maxRange*(1 - gate.limit)}.
#' If \code{gate.limit} is larger than 1, is must be a numeric value of  
#' length two that will be interepretted as the lower and upper limits of
#' the gate. If \code{gate.limit = TRUE}, \code{gate.limit} will be set to
#' \code{c(0.02, 0.98). If \code{gate.limit = FALSE}, the gate will span 
#' the full range of data.
#' 
#' Since this is typically used for cell cycle analysis with DNA content,
#' the default parameters are \code{xchan = "FL2.A"} and
#' \code{ychan = "FL2.H"}
#'  
#' @param x Either a matrix from \code{\link[Biobase]{exprs}}, a
#'   \code{flowFrame}, or a \code{flowSet} with the data
#' @param xchan,ychan Character values identifying the data to search for
#'   single cells with default values of \code{"FL2.A"} and \code{"FL2.H"}
#' @param zero.intercept If \code{TRUE} (default), the regression used to
#'   identify the singles will be forced through the origin
#' @param dx,dy The width and height of rectangular gate, default value of
#'   0.06. If only \code{dx} is specified, \code{dy} will be assigned the
#'   same value 
#' @param xquant A numeric value of length two to specify the upper and
#'   lower percentile range of the data in \code{xchan} to be used for the
#'   regression. If \code{NULL}, a default value of \code{c(0.02, 0.98)}
#'   will be used
#' @param gate.limit Numeric vector of length 1 or 2 \emph{or} a 
#'   \code{logical} value to specify limits on the returned gate. See
#'   Details for more information
#' @param filterId Character string with default value of \code{"singlets"}
#' @param groupFilterId Character string with default value of
#'   \code{"singletsGateList"}
#'
#' @return
#' A rectangular \code{polygonGate} or list of these gates.
#' 
#' @import flowCore
#' @importFrom MASS lqs
#'
#' @export
#'
singletGate <- function(x, xchan = "FL2.A", ychan = "FL2.H",
		zero.intercept = TRUE, dx = 0.06, dy = dx, xquant = NULL,
		gate.limit = TRUE,
		filterId = "singlets", groupFilterId = "singletsGateList")
{
# default quantile range
	QUANT_RANGE <- c(0.02, 0.98)

# argument checks
	if(!(is(x, "flowFrame") | is(x, "flowSet") | is(x, "matrix")))
		stop("require either a flowFrame, flowSet or matrix")
	nms <- colnames(x)
	if (!xchan %in% nms)
		stop('"', xchan, '" not found in ', deparse(substitute(x)))
	if (!ychan %in% nms)
		stop('"', ychan, '" not found in ', deparse(substitute(x)))
	if (is.null(xquant))
		xquant <- QUANT_RANGE
	if (any(xquant) > 1 | any(xquant) < 0)
		stop("'xquant' must be within [0,1]")

# process gate.limit to have two values
	if (identical(gate.limit, TRUE))
		gate.limit <- QUANT_RANGE
	else if (identical(gate.limit, FALSE))
		gate.limit <- c(-Inf, Inf)
	else if (is.numeric(gate.limit) && all(gate.limit <= 1))
		if (length(gate.limit) == 1)
			gate.limit <- c(gate.limit, 1 - gate.limit)
		else
			gate.limit <- range(gate.limit)
	else if (is.numeric(gate.limit) && all(gate.limit > 1))
		if (length(gate.limit) != 2)
			stop("values of gate.limit > 1 must be a vector of length 2")
		else
			gate.limit <- range(gate.limit)
	else
		stop("unable to use values in 'gate.limit'")

# get range of data in xchan
	if (is(x, "flowFrame"))
		xRange <- range(range(x[, xchan], type = "data"))
	else if (is(x, "flowSet"))
		xRange <- range(fsApply(x[, xchan], range, type = "data"))
	else if (is(x, "matrix"))
		xRange <- range(x[, xchan], na.rm = TRUE)
	else
		stop("Whoa! Should not have arrived here...")

# working function expects matrix with columns names in xchan, ychan
# xquant holds quantiles to limit data for regression
# xRange has the range of data for xchan in each frame
# gate.limit is two values, each less than 1 or each greater than 1
#
	.singletGate <- function(mat, xchan, ychan, dx, dy, xquant,
			xRange, gate.limit, filterId)
	{
	# pre-pocess data into data.frame
		dat <- data.frame(mat[, c(xchan, ychan)])
		names(dat) <- c("x", "y")

	# subset data for regression (use 
		xlimits <- quantile(dat$x, xquant, na.rm = TRUE)
		ylimits <- quantile(dat$y, yquant, na.rm = TRUE)
		dat <- subset(dat, x > min(xlimits) & x < max(xlimits))
		dat <- subset(dat, y > min(ylimits) & y < max(ylimits))

		if (zero.intercept)
			fm <- try(MASS::lqs(y ~ x + 0, data = dat), silent = TRUE)
		else { # no zero-intercept, much slower fit
			if (nrow(dat) > 1600)
				idx <- sample(seq_len(nrow(dat)), 800) # to speed up regression
			else
				idx <- seq_len(nrow(dat))
			fm <- try(MASS::lqs(y ~ x, data = dat, subset = idx), silent = TRUE)
		}
		if (class(fm) == "try-error") { # lqs failed, return zero-width gate 
			x.vertices <- c(xlimits[1], xlimits[1], xlimits[2], xlimits[2])
			y.vertices <- c(ylimits[1], ylimits[1], ylimits[2], ylimits[2])
			filerIdValue <- "unfiltered"
		}
##### r * diff(range(a)) + 1) + min(a)
		else { # proceed by setting limits from gate.limit
			xp <- xquant
			if (identical(gate.limit, FALSE)
				xp <- c(1, 1)
			else if (identical(gate.limit, TRUE))
				xp <- xquant
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
		bd <- cbind(x.vertices, y.vertices)
		colnames(bd) <- c(xchan, ychan)

		return(polygonGate(.gate = bd, filterId = filterIdValue))
	}

# dispatch function on nature of argument 'x'
	if (class(x) == "matrix")
		ret <- .singletGate(x, xchan, ychan, dx, dy, xquant, yquant,
						gate.limit, cutoff, maxRange, filterId)
	else if (class(x) == "flowFrame")
		ret <- .singletGate(exprs(x), xchan, ychan, dx, dy, xquant, yquant,
						gate.limit, cutoff, maxRange, filterId)
	else if (class(x) == "flowSet") {
		ret <- fsApply(x, .singletGate, xchan, ychan, dx, dy, xquant, yquant,
						gate.limit, cutoff, maxRange, filterId, use.exprs = TRUE)
		ret <- filterList(ret, filterId = groupFilterId)
	}
	else
		stop("cannot process argument of class \"", class(x), "\"")

	return(ret)
}
