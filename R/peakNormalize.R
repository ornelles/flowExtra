#' Translate and Scale G1 and G2 Peaks
#'
#' Normalize "G1" and "G2" peaks in DNA profile data
#'
#' @md
#' @details
#' This function attempts to translate and scale linear DNA profile
#' data to align peaks corresponding to a "G1" and "G2/M" population
#' of cells. Data will be adjusted by translating the data in order
#' to place the "G1" peak at the value in `g1`. If `scale = TRUE`,
#' data will be scaled to place the "G2/M" peak at the value
#' in `g2`. If `g2` is missing, the value for `g2` is given by `ratio *
#' g1`. Although this ratio should theoretically be 2.0, it is
#' empirically closer to 1.9.
#'
#' The position of the actual peaks is specified in the optional argument
#' `peaks`. If provided, this can either be a vector of "G1" peak
#' positions or a matrix of "G1" and "G2/M" peak positions. A larger matrix
#' will be accepted, such as one produced with [peakFind()], but only the
#' first two columns will be used and treated as "G1" and "G2" positions,
#' respectively. Missing values in either "G1" or "G2"
#' will be replaced by a value derived from `ratio`. If the argument `
#' peaks` is `NULL`, the function [peakFind()] will be called with
#' arguments `range.search` and `...` to identify the peaks.
#'
#' If the `logical` argument `scale = FALSE`, a complete set of "G1"
#' values must be provided and data will **not** be scaled. The values
#' will only be translated to place the "G1" peak at value specified
#' in the argument `g1`. Note that this **should not** be set to `TRUE`
#' for data with *bona fide* peaks that are not true "G1" or "G2" peaks
#' such those found in synchronized populations of cells as in the
#' example.
#'
#' Transformed values will be trimmed to the range specified by
#' `limits`, which has a default value of `c(0, 1023)`. The full range of
#' original values will be used if `limits = NULL` and no trimming will
#' occur if `limits = FALSE`.
#'
#' @param fs A `flowSet` with data to be transformed in channel in `chan`
#' @param peaks An optional matrix of "G1" and "G2/M" peak values *or* a
#'   vector of "G1" values. The default value of `NULL` will use
#'	 [peakFind()] to identify the peaks within the range specified by
#'   `range.search`
#' @param g1 Value of the rescaled G1 peak (`200`)
#' @param g2 Value of the rescaled G2 peak (`ratio * g1`)
#' @param ratio G2 to G1 ratio, used to scale data and place missing peaks
#'   (`1.90`)
#' @param chan Character vector to specify the channel to normalize
#'   (`"FL2.A"`)
#' @param limits Numeric vector of length 2 specifying the limits for the
#'   transformed data (default value of `c(0, 1023)`). If `NULL`, the range
#'   will be the original limits of the data. If `FALSE` the range will be
#'   `-Inf` to `+Inf`.
#' @param range.search A numeric vector of length 2 specifying the range
#'   of values in the original data in which peaks will be accepted. The
#'   default value of `c(50, 500)` is chosen to ignore apoptotic values and
#'   polyploid peaks in typical data
#' @param scale A `logical` value to scale the data such that the "G1" and
#'   "G2/M" peaks occur at the values in `g1` and `g2`, respectively.
#'   See Details for more information
#' @param ... Additional arguments passed to [peakFind()] and then to
#'   [flowStats::curv1Filter()]
#' 
#' @return
#' A transformed and (optionally) scaled `flowSet`.
#' 
#' @import flowCore flowStats
#' 
#' @export
#' 
peakNormalize <- function(fs, peaks = NULL, g1 = 200, g2 = g1 * ratio,
		ratio = 1.90, chan = "FL2.A", limits = c(0, 1023),
		range.search = c(50, 500), scale = FALSE, ...)
{
# argument checks accept only 'flowSet'
	if (class(fs) != "flowSet")
		stop("flowSet required")
# check for channel
	if (!chan %in% colnames(fs))
		stop('"', chan, '" is not found in ', deparse(substitute(fs)))
# extract G1 and G2 peaks if not provided in 'peaks'
	if (is.null(peaks)) {
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

# adjust any invisible G1 or G2 peaks before possible scaling
	missingG2 <- is.na(peaks[,2])
	peaks[missingG2, 2] <- peaks[missingG2, 1] * ratio
	missingG1 <- is.na(peaks[,1])
	peaks[missingG1, 1] <- peaks[missingG1, 2] / ratio

# assign limits to transformed values
	if (is.null(limits) || limits == TRUE)
		lim <- range(fsApply(fs, function(v) range(exprs(v)[,chan],
			na.rm = TRUE)), na.rm = TRUE)
	else if (is.numeric(limits) && length(limits) == 2)
		lim <- limits
	else
		lim <- c(-Inf, Inf)

# create copy of flowSet in order to change values 'in situ'
	x <- Subset(fs, TRUE)
# act only on non-NA G1 peaks
	sel <- which(!is.na(peaks[, 1]))
###
### NO - must scale by multiplying. Perhaps only G1 peak should
### be adjusted!! Scale option should be for only G1 or both G1/G2
###
# translate and (optionally) scale data
	for (i in sel) {
		v <- exprs(x[[i]])[, chan]
		if (scale == TRUE) # set scale factor
			a <- (g2 - g1)/(peaks[i, 2] - peaks[i, 1])
		else # no scaling
			a <- 1
		v <- a * v # scale
		v <- v + (g1 - a * peaks[i, 1]) # translate
		v <- pmax(v, lim[1]) # trim
		v <- pmin(v, lim[2])
		exprs(x[[i]])[, chan] <- v # replace with transformed values
	}
	return(x)
}
