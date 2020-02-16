#' Find G1 and G2/M peaks
#' 
#' Find G1 and G2/M (and other) peaks in linear values
#' 
#' @md
#' @details
#' This function identifies the probable location of the G1 and G2/M peaks
#' in the channel identified by `chan` in either a `flowFrame` or `flowSet`
#' or "raw" values as described below.
#' Areas of local high density are first identified by
#' [flowStats::curv1Filter()] with arguments in `bwFac` and `gridsize`.
#' The values are the default values for by `curv1Filter()`. 
#' The peak location within each local area of high density is
#' identified from a kernel density estimate [stats::density()] using
#' arguments passed in the list `darg`. Peak positions are returned as
#' an array of two or more columns.
#' 
#' If `range.search` is `NULL`, peaks are accepted in the range of data
#' spanning 5th to 95th percentile or as specified by the percentile values
#' in `probs`. For well-behaved data, `range.search` can be set
#' to `c(50, 500)` to exclude apoptotic and polyploid populations.
#'
#' If exactly two peaks are found, they will be labeled `"G1"` and `"G2"`.
#' If a mixture of one and two peaks are found, the single peak will be
#' treated as the `"G1"` peak and the `"G2"` peak will be assigned an `NA`
#' value. If only one peak is found or if more than two peaks
#' are found, the peaks will be labeled as `"peak1"`, `"peak2"`, `"peak3"`,
#' etc. A warning message will report finding other than 2 peaks.
#'
#' The argument `x` can also be a `matrix` or `numeric` vector of
#' "raw" values which will be
#' coerced into a matrix. If the name of a column (or vector) is `chan`,
#' that parameter will be analyzed. If the `matrix` or `vector` is not
#' named, the first column will be analyzed.  
#'
#' @param x A `flowFrame` or `flowSet` with the parameter named in `chan`
#'   *or* a `numeric` vector or matrix (see Details)
#' @param chan The name of the data to be evaluated in `x`
#' @param darg List of arguments passed to the kernel density estimating
#'   function ([stats::density()]) used to find peaks
#' @param range.search An optional numeric vector of length 2 setting the
#'   range of values in `chan` to search for peaks. The default value of `NULL` 
#'   limits the search to the range defined by the quantiles specified
#'   in `probs`
#' @param probs Numeric vector of length 2 defining the lower and upper
#'   quantiles to set the value for `range.search` if not provided
#' @param bwFac,gridsize Numeric values of length 1 and 2 used by
#'   [flowStats::curv1Filter()] to identify regions of 1-D curvature in
#'   the data (`bwFac` increased from default of 1.2 to 2)
#' 
#' @return
#' An 2-D array with row names obtained from [flowCore::identifier()]
#' and column names determined as described above.
#' 
#' @examples 
#'  fs <- readSet(system.file("extdata", "synch/", package = "flowExtra"))
#'  peakFind(fs) # find all by default
#'  peakFind(fs, range.search = c(50, 500)) # find probable G1 and G2
#'
#' @importFrom flowCore filter exprs fsApply
#' @importFrom flowStats curv1Filter
#' 
#' @export
#' 
peakFind <- function(x, chan = "FL2.A", darg = list(bw = "nrd0", n = 512),
		range.search = NULL, probs = c(0.05, 0.95), bwFac = 2,
		gridsize = rep(401, 2))
{
	# argument check
		if (is(x, "flowSet") | is(x, "flowFrame")) {
			if (!chan %in% colnames(x))
				stop('"', chan, '" is not in ', deparse(substitute(x)))
		}
		else if (is(x, "numeric")) { # assembly dummy flowFrame
			x <- as.matrix(x)
			sel <- which(colnames(x) == chan)[1]
			if (!is.na(sel))
				x <- x[, sel]
			else {
				x <- x[, 1, drop = FALSE]
				colnames(x) <- chan <- "raw"
			}
			x <- flowFrame(x)
			identifier(x) <- "raw"
		}
		else
			stop("'", deparse(substitute(x)), "' must be a flowFrame, a flowSet,
				a numeric vector, or a matrix")

# working function to perform peak finding on a flowFrame
	.peakFind <- function(x, chan, range.search, darg)
	{
	# identify high-density regions with curv1Filter
		c1f <- curv1Filter(list(chan), bwFac = bwFac, gridsize = gridsize)
		res <- flowCore::filter(x, c1f)
	# helper function to pick peak position
		.localMax <- function(x) {
			kd <- do.call(stats::density, c(list(x), darg))
			kd$x[which.max(kd$y)] }
	# find the peak within each high density region
		peaks <- tapply(exprs(x)[, chan], res@subSet, .localMax)
		peaks <- peaks[-1] # drop "rest"
	# assign value to 'range.search' if needed
		if (is.null(range.search))
			range.search <- quantile(exprs(x)[, chan], probs = probs, na.rm = TRUE)
	# drop NA values, drop peaks outside of 'range.search'
		peaks <- peaks[!is.na(peaks)]
		sel <- peaks > range.search[1] & peaks < range.search[2]
		peaks <- peaks[peaks > range.search[1] & peaks < range.search[2]]
	# add NA value for only one peak
		if (length(peaks) == 1) { # add NA for G2 peak
			peaks[2] <- NA
			warning(identifier(x), ": one peak found", call. = FALSE)
		}
		else if (length(peaks) > 2) # report
			warning(identifier(x), ": ", length(peaks), " peaks found", call. = FALSE)
		return(signif(peaks, 5))	# return x-position of peaks
	}

# dispatch working function according to arguments
	if (class(x) == "flowFrame") {
		ans <- t(.peakFind(x, chan, range.search, darg))
		rownames(ans) <- identifier(x)
	}
	else if (class(x) == "flowSet") {
		ans <- fsApply(x, .peakFind, chan, range.search, darg, simplify = FALSE)
		n <- max(lengths(ans), na.rm = TRUE)
		ans <- t(sapply(ans, function(x) c(x, rep(NA, n))[1:n]))
	}
	else
		stop("SHOULDN'T BE HERE with argument of class ", class(x))

# clean up return value
	if (ncol(ans) == 2 && all(is.na(ans[, 2])))
			ans <- ans[, 1, drop = FALSE]
	if (ncol(ans) == 2 && all(ans[, 2]/ans[, 1] >= 1.9, na.rm = TRUE))
		colnames(ans) <- c("G1", "G2")
	else
		colnames(ans) <- paste0("peak", seq_len(ncol(ans)))
	return(ans) # array of 2 or more dimensions
}
