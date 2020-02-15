#' Determine Extreme Values
#' 
#' Determine minimum and maximum values from a smoothed curve through x-y points.
#' 
#' @param x,y x- and y-values of the curve to search for extreme values
#' @param first if \code{TRUE}, return only the first extreme value
#' @param spar smoothing parameter for \code{smooth.spline} between 0 and 1
#' 
#' @return
#' An \code{invisible} list identifying each extreme value found in the smoothed
#' curve. The components of the list include \code{x} and \code{y} for the
#' x- and y-values, \code{sign} which is +1 for a maximum and -1 for a minimum,
#' and the smoothing model returned by \code{smooth.spline()} in \code{model}.
#' 
#' @details
#' Multiple \code{x} values are collapsed into a single value and the 
#' corresponding \code{y} value treated as the mean of multiple values.
#' The resulting \code{x} and \code{y}
#' values are smoothed by \code{smooth.spline()} with the smoothing parameter
#' \code{spar}. If \code{spar} is \code{NULL}, the value will be determined by 
#' the algorithm and returned in \code{model}. This is the same code used in
#' the \code{ribofrag} package although neither are exported.
#' 
#' Neither \code{x} nor \code{y} can have missing values.
#' 
#' The endpoints of the \code{x, y} points are included in the search for extreme
#' values.
#' 
#' The function is called "\code{zero.cross}" because the point at which the first
#' derivative crosses zero is used to determine the extreme. The fit and choice of
#' \code{spar} can be evaluated by \code{plot(x, y); lines(zero.cross(x, y)$model)}.
#' 
#' @examples
#'   set.seed(123)
#'   y_ex <- 100*cumsum(rnorm(250))
#'   x_ex <- seq_along(y_ex)
#'   ans <- zero.cross(x_ex, y_ex, spar = 0.5)
#'   as.data.frame(ans[1:3])
#'   plot(x_ex, y_ex)
#'   lines(ans$model)
#'   abline(v = ans$x, col = ifelse(ans$sign > 0, "black", NA))
#' 
#'   print(zero.cross(1:10, 1:10))
#' 
zero.cross <- function(x, y, first = FALSE, spar = NULL) {
	keep <- complete.cases(x, y)	# smooth.spline can't handle missing data
	x <- x[keep]
	y <- y[keep]
	if (length(y) <5)
		stop("requires at least 4 points")
	if (anyDuplicated(x)) {			# updated
		warning("duplicated x observations averaged")
		y <- as.numeric(tapply(y, x, mean))
		x <- unique(x)
	}
	fms <- smooth.spline(x, y, spar = spar)
	fmp <- predict(fms, deriv=1)
	sign1 <- sign(fmp$y[-length(fmp$y)])
	sign2 <- sign(fmp$y[-1])
	zc.index <- which(sign1 * sign2 == -1)
	
	if(length(zc.index) == 0) {	# no curvature, collect absolute min and max
		zc.index[1] <- which.min(y)
		zc.index[2] <- which.max(y)
		sign1[zc.index] <- c(-1,1)
	}
	zc.index <- sort(zc.index)
	if (first)
		zc.index <- zc.index[1]
	xx <- x[zc.index]
	yy <- y[zc.index]
	sign <- sign1[zc.index]
	invisible(list(x=xx, y=yy, sign=sign, model=fms))
}
