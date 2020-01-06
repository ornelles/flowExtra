#' Translate and Scale Linear Data
#'
#' Scale linear DNA content to have identical "G1" peak values
#'
#' @md
#' @details
#' This function scales and translates linear DNA profile
#' data to align peaks corresponding to a "G1" population
#' of cells. Data will be adjusted by scaling the linear data in order
#' to place the "G1" peak at the value in `g1`.
#'
#' The position of the actual "G1" population (and "G2/M" population)
#' is provided in the optional argument `peaks`. This can either be a
#' vector of "G1" peak positions or a matrix of "G1" and "G2/M" peak
#' positions. A larger matrix will be accepted, such as one produced
#' with [peakFind()], but only the first two columns will be used and
#' treated as "G1" and "G2" positions, respectively. Missing values in
#' "G1" that have a "G2/M" value will be replaced by a value derived
#' from `ratio`. If the argument `peaks` is missing, the function
#' [peakFind()] will be called with arguments `range.search` and `...`
#' to identify the peaks.
#'
#'
#' Transformed values will be trimmed to the range specified by
#' `limits`, uses the instrument range for `chan` if `limits = NULL`.
#' No trimming will occur if `limits = FALSE`.
#'
#' @param fs A `flowSet` with data to be transformed in channel in `chan`
#' @param peaks An optional matrix of "G1" and "G2/M" peak values *or* a
#'   vector of "G1" values. If missing, [peakFind()] to identify the peaks
#'   within the range specified by `range.search`
#' @param g1 Value of the rescaled G1 peak (`200`)
#' @param ratio G2 to G1 ratio to place missing peaks "G1" peaks given
#'   a corresponding "G2/M" peak, default value of `1.90`
#' @param chan Character vector to specify the channel to normalize (`"FL2.A"`)
#' @param limits Optional numeric vector of length 2 specifying the limits
#'   for the transformed data. If `NULL`, the instrument range for data in
#'   `chan` will be used. If `FALSE` the range will be `-Inf` to `+Inf`.
#' @param range.search A numeric vector of length 2 specifying the range
#'   of values in the original data in which peaks will be accepted. The
#'   default value of `c(50, 500)` is chosen to ignore apoptotic values and
#'   polyploid peaks in typical data
#' @param ... Additional arguments passed to [peakFind()] and then to
#'   [flowStats::curv1Filter()]
#' 
#' @return
#' A transformed and scaled `flowSet`.
#' 
#' @import flowCore flowStats
#' 
#' @export
#' 
peakNormalize <- function(fs, peaks, g1 = 200, ratio = 1.92 chan = "FL2.A",
		limits = NULL, range.search = c(50, 500), ...)
{
# argument checks accept only 'flowSet'
	if (class(fs) != "flowSet")
		stop("flowSet required")
# check for channel
	if (!chan %in% colnames(fs))
		stop('"', chan, '" is not found in ', deparse(substitute(fs)))
# extract G1 and G2 peaks if not provided in 'peaks'
	if (missing(peaks)) {
		cat("Identifying peaks from", range.search[1], "to", range.search[2], "\n")
		flush.console()
		peaks <- fsApply(fs, flowExtra::peakFind, chan = chan,
			range.search = range.search, ...)
	}
# allow peaks to be a vector of G1 peaks or a "single" G1 peak
	if (is(peaks, "vector"))
		if (length(peaks) > 1)
			peaks <- cbind(peaks, NA)
		else
			peaks <- cbind(rep(peaks, length(fs)), NA)
	peaks <- peaks[, 1:2] # use first two peaks

# adjust any invisible G1 peaks with a G2 counterpart
	missingG1 <- is.na(peaks[,1])
	peaks[missingG1, 1] <- peaks[missingG1, 2] / ratio

# assign limits to transformed values
	if (is.null(limits) || identical(limits, TRUE))
		lim <- range(range(fs[[1]][, chan], type = "instrument"))
	else if (identical(limits, FALSE))
		lim <- c(-Inf, Inf)
	else if (is.numeric(limits) && length(limits) == 2)
		lim <- limits
	else 
		stop("unable to use value for 'limits'")

# create copy of flowSet in order to change values 'in situ'
	x <- Subset(fs, TRUE)
	sel <- which(!is.na(peaks[, 1])) # act only on non-missing G1 peaks
	k <- peaks[,1]/g1

# scale data to align G1 peak with 'g1'
	for (i in sel) {
		v <- exprs(x[[i]])[, chan]
		v <- a * k # scale
		v <- pmax(v, lim[1]) # trim
		v <- pmin(v, lim[2])
		exprs(x[[i]])[, chan] <- v # replace with transformed values
	}
	return(x)
}
