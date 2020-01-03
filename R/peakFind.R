#' Find G1 and G2 peaks
#' 
#' Find G1 and G2 (and other) peaks in FL2.A values
#' 
#' @md
#' @details
#' This function identifies the probable location of the G1 and G2 peaks
#' in the `"FL2.A"` channel (or the channel
#' specified by the argument `chan`) in a `flowFrame` or `flowSet`. Areas
#' of high density are identified by the function defined in `curveFilter`. If 
#' this is `NULL`, the [flowStats::curv1Filter()] function will be used 
#' with optional arguments (`bwFac` or `gridsize`) passed in `...`.
#' 
#' Peak searching is limited to 0.5 to 99.5\% quantiles (default) or to the 
#' values specified in `probs`. If no more than two peaks are found, they will 
#' be labeled `"G1"` and `"G2"`. A single peak will be labeled as `"G1"` with
#' an `NA` for the second peak, labeled as `"<G2>"`. If more peaks are found,
#' they will be labeled `"peak1"`, `"peak2"`, `"peak3"`, etc.

#' @param x A `flowFrame` or `flowSet` with the parameter named in `chan`
#' @param curveFilter Optional function that returns a `multipleFilterResult`
#'   defining the regions of 1-D curvature in `chan`. If `NULL`, the
#'   [flowStats::curv1Filter()] will be used
#' @param chan Name of the data to be evaluated in `x`, default of `"FL2.A"`
#' @param searchFun Optional function to determine the position of each
#'   peak in defined by the argument `curvFilter`. If `NULL`, the
#'   [stats::density()] function will be used with the bandwith adjusted
#'   by the option `adj = 2`
#' @param range.search Optional numeric vector of length 2 defining the
#'   lower and upper limits to accept peaks. If `NULL` the search will be
#'   limited to the values defined by the quantiles specified in `probs`
#' @param probs Numeric vector of length 2 defining the lower and upper
#'   quantiles for valid peak values, default value of 0.5 to 0.95\%
#'  (`probs = c(0.005, 0.995)`)
#' @param ... Additional arguments to be passed to
#'   [flowStats::curv1Filter()]
#' 
#' @return
#' An 1-D or 2-D array with rownames obtained from [flowCore::identifier()]
#' and column names determined as described above.
#' 
#' @import flowCore flowStats
#' 
#' @export
#' 
peakFind <- function(x, curveFilter = NULL, chan = "FL2.A",
		searchFun = NULL, range.search = NULL, probs = c(0.005, 0.995), ...)
{
	# assign values to arguments
	if (is.null(curveFilter))
		curveFilter <- curv1Filter(chan, ...)
	if (is.null(searchFun)) # kernel density with 2x default bandwidth
		searchFun <- function(x, adj = 2) {
			kd <- stats::density(x, adj = adj, na.rm = TRUE)
			kd$x[which.max(kd$y)]	# return first 'x' at maximum
		}
# working function to perform peak finding on a flowFrame
	.peakFind <- function(x, curveFilter, chan, searchFun, range.search)
	{
	# apply filter to define high density regions in one dimension (curv1Filter)
		res <- filter(x, curveFilter)
		peaks <- tapply(exprs(x)[, chan], res@subSet, searchFun)[-1]
	# assign value to 'range.search' if needed
		if (is.null(range.search))
			range.search <- quantile(exprs(x)[, chan], probs = probs, na.rm = TRUE)
	# drop NA values, drop peaks outside of 'range.search'
		peaks <- peaks[!is.na(peaks)]
		sel <- peaks > range.search[1] & peaks < range.search[2]
		peaks <- peaks[peaks > range.search[1] & peaks < range.search[2]]
	# assign best guess names to peaks
		if (length(peaks) == 2) # assume G1 and G2 peak found
			names(peaks) <- c("G1", "G2")
		else if (length(peaks) == 1) { # assume only G1 peak found
			peaks[2] <- NA
			names(peaks) <- c("G1", "<G2>")
			warning(identifier(x), ": one peak found", call. = FALSE)
		}
		else {
			warning(identifier(x), ": ", length(peaks), " peaks found", call. = FALSE)
			names(peaks) <- paste("peak", 1:length(peaks), sep="")
		}
		return(peaks)	# return x-position of peaks
	}

# dispatch working function according to arguments
	if (class(x) == "flowFrame") {
		ans <- t(.peakFind(x, curveFilter, chan, searchFun, range.search))
		rownames(ans) <- identifier(x)
	}
	else if (class(x) == "flowSet") {
		ans <- fsApply(x, .peakFind, curveFilter, chan, searchFun, range.search,
				simplify = FALSE)
		ans.range <- range(sapply(ans, length), na.rm = TRUE)
		cnames <- names(ans[[which(sapply(ans, length) == ans.range[2])[1]]])
		ans <- lapply(ans, function(x) c(x, rep(NA, ans.range[2] - length(x))))
		ans <- t(as.data.frame(ans))
		colnames(ans) <- cnames
	}
	else
		stop("unable to operate on argument of class ", class(x))
	return(ans) # either 1-D or 2-D array
}
