#' Translate and Scale G1 and G2 Peaks
#'
#' Normalize "G1" and "G2" peaks in DNA profile data
#'
#' @md
#' @details
#' This function attempts to translate a linear DNA profile data to align
#' peaks corresponding to a "G1" population of cells. The "G2" peak will
#' be adjusted by scaling to be positioned at the value given by `ratio *
#' g1` where `g1` is the position of the "G1" peak. Although this ratio
#' should theoretically be 2.0, it is empirically closer to 1.9.
#'
#' The position of the actual peaks is specified in the optional argument
#' `peaks`. If provided, this can either be a vector of "G1" peak
#' positions or a matrix of "G1" and "G2" peak positions. A larger matrix
#' can be accepted, such as one produced with [peakFind()], but only the
#' first two columns will be used with "G1" and "G2" in the first and second
#' columns, respectively. Missing values in either "G1" or "G2"
#' will be replaced by a value derived from `ratio`. If the argument `
#' peaks` is `NULL`, the function [peakFind()] will be called with
#' arguments `range.search` and `...` to identify the peaks.
#'
#' If the `logical` argument `only.G1 = TRUE`, a complete set of "G1"
#' values must be provided and the data will only be translated to place
#' the "G1" peak at the selected value.
#'
#' The transformed values will be trimmed to the range specified by `
#' limits`, which has a default value of `c(0, 1023)`. The full range of
#' original values will be used if `limits = NULL` and no trimming will
#' occur if `limits = FALSE`.
#'
#' @param fs A `flowSet` with data to be transfomred in channel in `chan`
#' @param peaks An optional matrix of "G1" and "G2" peak values *or* a
#'   vector "G1" values. The default value of `NULL` will use
#'	 [peakFind()] to identify the peaks
#' @param g1 Value of the rescaled G1 peak (`200`)
#' @param g2 Value of the rescaled G2 peak (`ratio * g1`)
#' @param ratio G2 to G1 ratio, used to scale data and place missing peaks
#'   (`1.90`)
#' @param chan Character vector to specify the channel to normalize
#'   (`"FL2.A"`)
#' @param limits Numeric vector of length 2 specifying the limits for the
#'   transformed data (default value of `c(0, 1023)`). If `NULL`, the range
#'   will the the original limits of the data. If `FALSE` the range will be
#'   `-Inf` to `+Inf`.
#' @param range.search A numeric vector of length 2 specifying the range
#'   of values in the orignal data in which peaks will be accepted. The
#'   default value of `c(50, 500)` is chosen to ignore apoptotic values and
#'   polyploid peaks
#' @param only.G1 `Logical` value to prevent scaling by only using the
#'   "G1" peak to anchor the data
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
		range.search = c(50, 500), only.G1 = FALSE, ...)
{
# accept only 'flowSet'
	if (class(fs) != "flowSet") stop("flowSet required")

# extract G1 and G2 peaks if not provided in 'peaks'
	if (is.null(peaks)) {
		cat("Identifying peaks between ", range.search[1], " and ",
		range.search[2], "...")
		flush.console()
		peaks <- fsApply(fs, peakFind, chan = chan, range.search = range.search, ...)
	}
# assume vector of peaks specifies either a vector of G1 peaks or a "single" G1 peak
	if (is(peaks, "vector"))
		if (length(peaks) > 1)
			peaks <- cbind(peaks, NA)
		else
			peaks <- cbind(rep(peaks, length(fs)), NA)
	peaks <- peaks[, 1:2] # use first two peaks

# adjust any invisible G1 or G2 peaks
	missingG2 <- is.na(peaks[,2])
	peaks[missingG2, 2] <- peaks[missingG2, 1] * ratio
	missingG1 <- is.na(peaks[,1])
	peaks[missingG1, 1] <- peaks[missingG1, 2] / ratio

# assign limits to transformed values
	if (is.null(limits) || limits == TRUE)
		lim <- range(fsApply(fs, function(v) range(exprs(v)[,chan], na.rm=T)), na.rm=T)
	else if (is.numeric(limits) && length(limits) == 2)
		lim <- limits
	else
		lim <- c(-Inf, Inf)

# create copy of flowSet then change values 'in situ'
	x <- Subset(fs, TRUE)
	sel <- which(!is.na(peaks[,1])) # act only on non-NA G1 peaks
	for (i in sel) { # scale factor
		v <- exprs(x[[i]])[, chan]
		if (only.G1 == FALSE)
			a <- (g2 - g1)/(peaks[i,2] - peaks[i,1]) # scale factor
		else
			a <- 1
	# scale
		v <- a * v
	# translate
		v <- v + (g1 - a * peaks[i,1])
	# trim
		v[v < lim[1]] <- lim[1]
		v[v > lim[2]] <- lim[2]
		exprs(x[[i]])[, chan] <- v
	}
	return(x)
}
