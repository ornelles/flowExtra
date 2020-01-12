#' Translate and Scale Linear Data
#'
#' Scale linear DNA content to have identical "G1" peak values
#'
#' @md
#' @details
#' This function scales linear DNA content in `chan`
#' to align peaks corresponding to a "G1" population
#' of cells. Data will be adjusted by scaling the values
#' in each `flowFrame` to place the "G1" peak at the value in `g1`.
#'
#' The position of the actual "G1" population (and "G2/M" population)
#' is provided in the optional argument `peaks`. This can either be a
#' vector of "G1" peak positions or a matrix of "G1" and "G2/M" peak
#' positions. A larger matrix will be accepted, such as one produced
#' with [peakFind()], but only the first two columns will be used and
#' treated as "G1" and "G2/M" positions, respectively. Missing values in
#' "G1" that have a "G2/M" value will be replaced by a value derived
#' from `ratio`. If the argument `peaks` is missing, the function
#' [peakFind()] will be called with arguments `range.search` and `...`
#' to identify the peaks.
#'
#'
#' Transformed values will be trimmed to the range specified by
#' `limits` or the instrument range for `chan` if `limits = NULL`.
#' No trimming will occur if `limits = FALSE`.
#'
#' @param fs A `flowSet` with data to be transformed in channel in `chan`
#' @param chan Character vector to specify the channel to normalize (`"FL2.A"`)
#' @param peaks An optional matrix of "G1" and "G2/M" peak values *or* a
#'   vector of "G1" values. If missing, [peakFind()] to identify the peaks
#'   within the range specified by `range.search`
#' @param g1 Value of the rescaled G1 peak (`200`)
#' @param ratio G2 to G1 ratio to place missing peaks "G1" peaks given
#'   a corresponding "G2/M" peak, default value of `1.90`
#' @param limits Optional numeric vector of length 2 specifying the limits
#'   for the transformed data. If `NULL`, the instrument range for data in
#'   `chan` will be used. If `FALSE` the range will be `-Inf` to `+Inf`.
#' @param range.search A numeric vector of length 2 specifying the range
#'   of values in the original data in which peaks will be accepted. The
#'   default value of `c(50, 500)` is chosen to ignore apoptotic values and
#'   polyploid peaks in typical data
#' @param ... Optional arguments passed to [peakFind()] and then to
#'   [flowStats::curv1Filter()] including `range.search`, `probs`,
#'   `bwFac` and `gridsize`
#' 
#' @return
#' A transformed and scaled `flowSet`.
#' 
#' @examples
#'  fs <- readSet(system.file("extdata", "synch/", package = "flowExtra"))
#'  bf <- boundaryFilter("FL2.A")
#'  lg <- linearGate(fs)
#'  fs <- Subset(fs, bf & lg)
#'  o1 <- dnaplot(~FL2.A, fs, xlim = c(120, 380), main = "Raw", plot = FALSE)
#'  pks <- peakFind(fs, range.search = c(100, 400))
#'  fs.adj <- peakAdjust(fs, peaks = pks, g1 = 180)
#'  o2 <- dnaplot(~FL2.A, fs.adj, xlim = c(120, 380), main = "Adjusted", plot = FALSE)
#'  plot(o1, split = c(1, 1, 2, 1), more = TRUE)
#'  plot(o2, split = c(2, 1, 2, 1), more = FALSE)
#'
#' @importFrom flowCore fsApply Subset exprs
#' 
#' @export
#' 
peakAdjust <- function(fs, chan = "FL2.A", peaks, g1 = 200, ratio = 1.92, 
		limits = NULL, ...)
{
# argument checks accept only 'flowSet'
	if (class(fs) != "flowSet")
		stop("flowSet required")
# check for channel
	if (!chan %in% colnames(fs))
		stop('"', chan, '" is not found in ', deparse(substitute(fs)))
# extract G1 and G2 peaks if not provided in 'peaks'
	if (missing(peaks)) {
		dots <- list(...)
		if (length(dots) == 0)
			cat("Identifying peaks with peakFind...\n")
		else
			cat("Identifying peaks with additional arguments for peakFind...\n")
		flush.console()
		peaks <- fsApply(fs, peakFind, chan = chan, ...)
	}
# allow peaks to be a vector of G1 peaks or a "single" G1 peak
	if (is(peaks, "vector"))
		if (length(peaks) > 1)
			peaks <- cbind(peaks, NA)
		else
			peaks <- cbind(rep(peaks, length(fs)), NA)
	peaks <- peaks[, 1:2] # use first two peaks

# adjust any invisible G1 peaks that have a G2 counterpart
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
	fsCopy <- Subset(fs, TRUE)

# calculate scale factor, only on non-missing G1 peaks
	scale.factor <- peaks[,1]/g1
	sel <- which(!is.na(peaks[, 1]))

# scale data to align G1 peak with 'g1'
	for (i in sel) {
		v <- exprs(fsCopy[[i]])[, chan]
		v <- v/scale.factor[i] # scale
		v <- pmax(v, lim[1]) # trim
		v <- pmin(v, lim[2])
		exprs(fsCopy[[i]])[, chan] <- v # replace with transformed values
	}
	return(fsCopy)
}
