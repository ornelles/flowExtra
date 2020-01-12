#' Create Gate for Single Cells
#'
#' Fit a narrow rectangular gate to select single cells in a 2-D
#' plot of \code{area} x \code{height}
#'
#' @details
#' This function creates a data-driven gate to select single cells from
#' \code{area} and \code{height} values collected on a linear scale. A
#' robust linear regression is fit to the data with
#' \code{\link[MASS]{lqs}}. This regression is used to create a narrow
#' rectangular gate centered about the single cells. If
#' \code{zero.intercept = TRUE}, the regression is forced through origin.
#' By default, the regression is applied to values within the 2.5 to 97.5
#' percentiles of data in channel \code{xchan}. Other percentiles can be
#' specified with the argument \code{xRange}.
#' 
#' The resulting \code{polygonGate(s)} will be limited to the range of
#' values specified by \code{gRange}. If \code{gRange = TRUE}, the gate
#' will be limited to the 2.5 to 97.5 percentile of all data in
#' channel \code{xchan}. If \code{gRange = FALSE}, the gate will
#' span the full range of data. If \code{gRange} is a numeric vector of
#' values less than or equal to 1, \code{gRange} will be treated as
#' the percentiles to limit the gate. If \code{gRange} is a numeric
#' vector with values greater than 1, \code{gRange} will be
#' interpreted as the lower and upper limits of the gate in the \code{xchan}
#' dimension.
#' 
#' Since this is typically used for cell cycle analysis with DNA
#' content, the default parameters are \code{xchan = "FL2.A"} and
#' \code{ychan = "FL2.H"}
#'
#' @param x Either a matrix from \code{\link[Biobase]{exprs}}, a
#'   \code{\link[flowCore]{flowFrame}}, or a \code{\link[flowCore]{flowSet}}
#' @param xchan,ychan Character strings identifying the data to search for
#'   single cells with default values of \code{"FL2.A"} and \code{"FL2.H"}
#' @param zero.intercept If \code{TRUE} (the default), force the regression
#'   used to identify the singlets through the origin
#' @param width The width (narrow dimension) of the rotated rectangular gate
#'   expressed as \emph{either} as a fraction of the \strong{instrument data 
#'   range} for the \code{xchan} parameter (if less than 1) \emph{or}
#'   the absolute width of the gate in the same units as \code{xchan}, if
#'   greater than 1
#' @param xRange A numeric value of length two to specify the upper and
#'   lower percentile of the data in \code{xchan} to be used for the
#'   regression. If \code{NULL}, a default value of \code{c(0.025, 0.975)}
#'   will be used
#' @param gRange A \code{logical} value \emph{or} a numeric vector of length
#'   2 to specify limits on the dimension of the returned gate. The limits
#'   can be specified by quantiles or absolute values. See the Details
#'   section for more information
#' @param filterId Character string with default value of \code{"singlets"}
#'   for successful fits. An unsuccessful fit is labeled
#'   \code{"unfiltered"}
#' @param groupFilterId Character string with default value of
#'   \code{"singletsGateList"}
#'
#' @return
#' A rectangular \code{polygonGate} or list of these gates.
#'
#' @examples
#' # Read synchronized cell data
#'   fs <- readSet(system.file("extdata", "synch/", package = "flowExtra"))
#'   lg <- linearGate(fs) # default settings
#'   xyplot(FL2.H ~ FL2.A, fs[1:4], filter = lg[1:4], stats = TRUE)
#'   lg2 <- linearGate(fs, width = 0.02, gRange = c(75, 750))
#'   xyplot(FL2.H ~ FL2.A, fs[1:4], filter = lg2[1:4], stats = TRUE)
#'
#' @importFrom flowCore polygonGate fsApply filterList
#' @importFrom MASS lqs
#'
#' @export
#'
linearGate <- function(x, xchan = "FL2.A", ychan = "FL2.H",
		zero.intercept = TRUE, width = 0.05, xRange = NULL,
		gRange = TRUE, filterId = "singlets",
		groupFilterId = "singletsGateList")
{
# default quantile range
	QUANTILES <- c(0.025, 0.975)

# argument checks
	if(!(is(x, "flowFrame") | is(x, "flowSet") | is(x, "matrix")))
		stop("require either a flowFrame, flowSet or matrix")
	nms <- colnames(x)
	if (!xchan %in% nms)
		stop('"', xchan, '" not found in ', deparse(substitute(x)))
	if (!ychan %in% nms)
		stop('"', ychan, '" not found in ', deparse(substitute(x)))

# convert xRange to a vector of two quantiles 
	if (identical(xRange, NULL))
		xRange <- QUANTILES
	else if (identical(xRange, TRUE))
		xRange <- QUANTILES
	else if (length(xRange) > 1 && all(xRange <= 1))
		xRange <- range(xRange)
	else
		stop("'xRange' must be NULL or a vector of values less than 1")

# convert gRange to a vector of two values 
	if (identical(gRange, TRUE))
		gRange <- QUANTILES
	else if (identical(gRange, FALSE))
		gRange <- c(0, 1)
	else if (length(gRange) > 1 && all(gRange >= 1))
		gRange <- range(gRange)
	else if (length(gRange) > 1 && all(gRange <= 1))
		gRange <- range(gRange)
	else
		stop("unable to use values in 'gRange'")

# get maximum range of actual data in xchan
	if (is(x, "flowFrame"))
		dataRange <- range(range(x[, xchan], type = "data"))
	else if (is(x, "flowSet"))
		dataRange <- range(fsApply(x[, xchan], range, type = "data"))
	else if (is(x, "matrix"))
		dataRange <- range(x[, xchan], na.rm = TRUE)
	else
		stop("Whoa! Should not have arrived here...")

# adjust width to absolute scale
	if (is.numeric(width) && width[1] <= 1)
		width <- abs(diff(dataRange)) * width[1]
	else if (is.numeric(width) && width[1] > 1 && width[1] <= max(dataRange))
		width <- width[1]
	else
		stop("unable to use value in 'width'")

# working function expects matrix with columns names in xchan, ychan
# xRange holds quantiles to limit data for regression
# gRange holds either quantiles or absolute values to limit gate
# dataRange holds the actual range of values for xchan
# width holds width of rectangle in units for xchan
#
	.linearGate <- function(mat, xchan, ychan, width,
			xRange, gRange, dataRange, filterId)
	{
	# collect data as data.frame of two variables
		dat <- data.frame(mat[, c(xchan, ychan)])
		names(dat) <- c("x", "y")

	# determine limits of x values for regression
		xlimits <- quantile(dat$x, xRange, na.rm = TRUE)
		dat <- subset(dat, x > min(xlimits) & x < max(xlimits))

	# perform robust regression
		if (zero.intercept)
			fm <- try(MASS::lqs(y ~ x + 0, data = dat), silent = TRUE)
		else { # no zero-intercept, much slower fit
			if (nrow(dat) > 1600)
				idx <- sample(seq_len(nrow(dat)), 800) # to speed up regression
			else
				idx <- seq_len(nrow(dat))
			fm <- try(MASS::lqs(y ~ x, data = dat, subset = idx), silent = TRUE)
		}
	# with successful regression
		if (class(fm) != "try-error") {
			if (all(gRange >= 1))
				xp <- gRange
			else
				xp <- gRange * diff(dataRange) + min(dataRange)
			yp <- predict(fm, newdata = data.frame(x = xp))
	
	# calculate vertices of rotated rectangle of width 'w'
			theta <- atan2(diff(yp), diff(xp))
			dx <- (width/2) * sin(theta)
			dy <- (width/2) * cos(theta)

			x.vertices <- c(xp[1] - dx, xp[1] + dx, xp[2] + dx, xp[2] - dx)
			y.vertices <- c(yp[1] + dy, yp[1] - dy, yp[2] - dy, yp[2] + dy)
			filterIdValue <- filterId
		}
		else { # lqs failed, return zero-width gate
			x.vertices <- c(min(dat$x), min(dat$x), max(dat$x), max(dat$x))
			y.vertices <- c(min(dat$y), min(dat$y), max(dat$y), max(dat$y))
			filerIdValue <- "unfiltered"
		}
	# assemble and return gate
		bd <- cbind(x.vertices, y.vertices)
		colnames(bd) <- c(xchan, ychan)
		return(polygonGate(.gate = bd, filterId = filterIdValue))
	}

# dispatch function on nature of argument 'x'
	if (class(x) == "matrix")
		ret <- .linearGate(x, xchan, ychan, width, xRange, gRange,
			dataRange, filterId)
	else if (class(x) == "flowFrame")
		ret <- .linearGate(exprs(x), xchan, ychan, width, xRange, gRange,
			dataRange, filterId)
	else if (class(x) == "flowSet") {
		ret <- fsApply(x, .linearGate, xchan, ychan, width, xRange, gRange,
			dataRange, filterId, use.exprs = TRUE)
		ret <- filterList(ret, filterId = groupFilterId)
	}
	else
		stop("cannot process argument of class \"", class(x), "\"")

	return(ret)
}