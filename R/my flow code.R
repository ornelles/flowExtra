#
# my flow functions - gradually being converted to package
#
# March 2013 - use of "add.text" as a function name clashed with S3/S4 convention
# changing to addText. Same thing with "add.gate" to be safe (addGate). Expanded
# mfi() to handle vector of ranges
#
## flowViz.par.set(theme = trellis.par.get(), reset = TRUE)
## xyplot(..., axis = axis.default)
##  path <- "C:/Users/Owner/Dropbox/docs/flow/2011 data/2011_0124 rpe1 synch roberta"
##
###
## Shift G1 peak without using peakNormalize
##
# library(flowExtra)
# flowViz.par.set(theme = trellis.par.get(), reset = TRUE)
# dropbox <- if (.Platform$OS.type == "unix") "~/Dropbox" else "~/../Dropbox"
# path <- file.path(dropbox, "/docs/flow/2011 data/2011_0124 rpe1 synch roberta")
# fs <- readSet(path)
# pd <- read.csv(file.path(path, "phenoData.csv"), stringsAsFactors = FALSE)
# pData(fs)$hpr <- factor(pd$hpr)
# pData(fs)$type <- factor(c("asynch", rep("synch", 14), "asynch"))
# pData(fs)$sample <- factor(pd$sample, levels = pd$sample[order(pd$hpr)])
# 
# fs <- Subset(fs, boundaryFilter("FL2.A"))
# lg <- linearGate(fs)
# fs2 <- Subset(fs, lg)
# fs3 <- peakAdjust(fs2, g1 = 170)
# dnaplot(sample ~ FL2.A, fs, main = "Raw data", xlim = c(0, 500))
# dnaplot(sample ~FL2.A, fs2, main = "Singlets", xlim = c(0, 500))
# dnaplot(sample ~FL2.A, fs3, main = "Correct G1 aligned", xlim = c(0, 500))
# 
#
# extract ellipsoid gate parameters from filter result for drawing or calculating
# example:
#	fres <- filter(flowSet, norm2Filter(c("SSC.H","FSC.H")))
#	fList <- getGate(fres)
#
# flowCore glitch: use type = "actual" for true values, type = "draw" to draw ellipse...
#
getGate <- function(f, type=c("actual", "draw")) {
	type <- match.arg(type)
	.fun <- function(v, type) {
		a <- v@filterDetails[[1]]
		if (type == "actual")
			mult <- a$radius
		else
			mult <- 1
		eg <- ellipsoidGate(.gate = a$cov * mult, mean = a$center, distance = a$radius)
		return(eg)
	}
	if (class(f) == "filterResultList")
		return(lapply(f, .fun, type))
	else if (class(f) == "logicalFilterResult")
		return(.fun(f, type))
	else
		stop("cannot use argument of class '", class(f), "'")
}
#
# rangeCut - similar to rangeGate but returns rectangleGate at desired cutoff.
# Pool data for all flowFrames in x. Values are trimmed at left and right edges
# by trim before determine the quantile at cut. Include the positive (true)
# or negative population. If absolute is TRUE, the trimming will be applied to the
# theoretical range of the instrument rather than the actual range of the data.
# Additional arguments in ... are passed to base plot function
#
rangeCut <- function(x, stain, cut = 0.99, plot = FALSE, positive = TRUE,
		trim = 0.001, absolute = TRUE, filterId = "defaultRectangleGate", ...)
{
    flowCore:::checkClass(x, c("flowFrame", "flowSet"))
    flowCore:::checkClass(stain, "character", 1)
    if (!stain %in% colnames(x))
        stop("'", stain, "' is not a valid parameter in this flowFrame")
    flowCore:::checkClass(plot, "logical", 1)
    if (is(x, "flowSet"))
        x <- as(x, "flowFrame")
	if (length(cut) != 1 || (cut < 0 | cut > 1))
		stop("cut must be a single value between 0 and 1")
    if (absolute) {
        vrange <- range(x, na.rm = TRUE)[, stain]
        vrange[is.nan(vrange)] <- 0
        vrange[1] <- min(vrange[1], min(exprs(x)[, stain], na.rm = TRUE))
    }
    else
		vrange <- range(exprs(x)[, stain], na.rm = TRUE)

    inc <- diff(vrange) * trim
    exprf <- char2ExpressionFilter(sprintf("`%s` > %s & `%s` < %s",
        stain, vrange[1] + inc, stain, vrange[2] - inc))
	tmp <- Subset(x, exprf)
	loc <- quantile(exprs(tmp)[, stain], cut)

    if (plot) {
		dens <- density(exprs(tmp)[, stain])
		plot(dens, main = paste(round(cut, 4), "breakpoint for parameter", stain),
			cex.main = 1, ...)
		lines(dens, ...)
		abline(v = loc, col = 2, lwd = 2)
		legend("topright", "breakpoint", fill = "red", bty = "n")
    }

	bounds <- if (positive)
		list(c(loc, Inf))
	else
		list(c(-Inf, loc))
	names(bounds) <- stain
	rectangleGate(bounds, filterId = filterId)
}
#
# Add text to an existing lattice plot where x can be an x, y point structure
#
addText <- function(x, y, labels, col = "red", cex = 2/3, adj = 0.5) {

	mat <- trellis.currentLayout("panel")
	N <- sum(mat != 0)
	if (N == 0)
		stop ("no lattice object found")

	if (missing(y)) {			# x should be a list
		if (all(c("x", "y") %in% names(x))) {
			xx <- x$x
			yy <- x$y
		}
	}
	else if (is.character(y)) {	# y is assumed to have labels
		labels <- y
		if (all(c("x", "y") %in% names(x))) {
			xx <- x$x
			yy <- x$y
		}
	}
	else {		# x and y better be numeric...
		xx <- as.numeric(x)
		yy <- as.numeric(y)
	}
	if (missing(labels))
		labels <- seq_along(xx)

	if (length(xx) == 1) {
		xx <- rep(xx, N)
		yy <- rep(yy, N)
	}
	if (length(xx) != N)
		stop("number of panels (", N, ") and number of points (", length(xx), ") differ")

	for (i in 1:N) {
		column <- which(apply(mat, 2, function(x) any(x == i)))
		row <- which(apply(mat, 1, function(x) any(x == i)))
		trellis.focus("panel", column = column, row = row, highlight = FALSE)
		do.call(panel.text, list(xx[i], yy[i], labels[i], cex = cex, col = col, adj = adj))
		trellis.unfocus()
	}
}

#
# Use base graphics to generate colorized density plot as for cell cycle
# Arguments
#	x		vector of values or a flowFrame or flowSet with a single flowFrame
# 	trim	drop highest (and/or) lowest values before density calculation
#	xrange	use in place of xlim to specify plotting limits
# 	scale	if TRUE, density (y) values will be scaled between 0 and 1
# 	breaks	vector of right endpoint values to divide up density plot
#			typical use is to specify the upper limit to subG1, G1, S, G2/M
#			default value is generated by quantile(): 0%,25%,50%,75%,100%
#	chan	character specifying channel if x is a flowFrame/flowSet
#
densityplot.color <- function(x, breaks = NULL, xrange = c(0, 2^10),
		col = NULL, darg = list(), main = NULL, scale = TRUE, chan = "FL2.A",
		trim = c("both", "low", "high", "none"), ...) {

# allow for flowSet (of one), flowFrame, or numeric vector

	if (class(x) == "flowSet" && length(x) != 1)
		stop("'x' as a flowSet must contain a single flowFrame")
	if (class(x) == "flowSet" && length(x) == 1)
		x <- exprs(x[[1]])[, chan]
	if (class(x) == "flowFrame")
		x <- exprs(x)[, chan]

# trim values as indicated

	if (is.logical(trim))
		trim <- ifelse(trim, "both", "none")
	else
		trim <- match.arg(trim)
	
	if (trim %in% c("low", "both"))
		x <- x[x != min(x)]
	if (trim %in% c("high", "both"))
		x <- x[x != max(x)]

# calculate density argument list in 'darg'

	h <- do.call("density", c(list(x = x), darg))

	if (is.null(col))
		col <- c("#E31CFF", "#E3AA00", "#118E00", "#1800C7", "#E30000", "#0080FF") # RGB
#		col <- c("#A771A6", "#F1B833", "#1E8E3B", "#7EBFE6", "#D92D2A", "#30539B") # CMYK

# assign xrange and (optionally) rescale density values

	if (is.null(xrange))
		xrange <- range(h$x, na.rm = T)
	if (scale == TRUE)
		h$y <- scale(h$y, center = FALSE, scale = max(h$y))

# plot density with xlim = xrange, using arguments in ...

	if (is.null(main))
		main <- "default density"
	plot(h, xlim = xrange, type = "n", main = main, ...)

# determine y-limits and adjust breaks to include endpoints

	ymin <- min(h$y, na.rm = T)
	if (is.null(breaks))
		breaks <- quantile(h$x, na.rm = TRUE)
	else
		breaks <- c(xrange[1], breaks, xrange[2])

# paint polygons for each zone

	for (n in 1:(length(breaks) - 1)) {
		id <- h$x >= breaks[n] & h$x < breaks[n + 1]
		xy <- cbind(x = c(h$x[id][1], h$x[id], rev(h$x[id])[1]),
					y = c(ymin, h$y[id], ymin))
		myCol <- col[(n - 1)%%length(col) + 1]
		polygon(xy, col = myCol, border = myCol, ...)
	}
}

#
# plot overlapping DNA profiles
#
# Arguments:
# 	fs - flowSet
#	sel - selection index for flowFrames
#	chan - channel to combine (FL2.A)
#	key.text - labels for auto.key
#	... - passed to density() rather than to densityplot() (use adj=0.2, etc.)
#
densityplot.overlap <- function(fs, sel, chan = "FL2.A",
	xlab = deparse(substitute(chan)),
	ylab = "Relative density",
	main = NULL,
	auto.key = NULL,
	key.text = NULL,
	columns = 2,
	n = 512,
	...)
{
	if (is.character(sel) & missing(chan)) {	# assume second argument IS chan
		chan <- sel
		sel <- TRUE
		cat("using entire flow set...\n")
	}

	temp <- fsApply(fs[sel], function(x) exprs(x)[, chan], simplify = FALSE)

	if (is.null(auto.key)) {
		if (is.null(key.text))
			auto.key <- list(columns = columns, size = 2.5)
		else
			auto.key <- list(text = as.character(key.text),
				columns = columns, size = 2.5)
	}
	if (is.list(auto.key) & !is.null(key.text))
		auto.key <- c(text = list(as.character(key.text)), auto.key)
	
	obj <- densityplot( ~ data, do.call(make.groups, temp), groups = which,
			plot.points = FALSE, auto.key = auto.key, n = n,
			xlab = xlab, ylab = ylab, main = main, ...)
	obj
}

#
# wrapper to select groups in flowset fs by virus, cell, and/or moi (in phenoData for fs)
#
overlap.wrapper.fun <- function(fs, virus, cell, moi, ...) {
	sel <- rep(TRUE, length(fs))
	if (!missing(virus))
		sel <- sel & pData(fs)$virus %in% virus
	if (!missing(cell))
		sel <- sel & pData(fs)$cell %in% cell
	if (!missing(moi))
		sel <- sel & pData(fs)$moi %in% moi

	if (!any(sel))
		warning("nothing to show for this grouping")
	else {
		key.text <- NULL
		if (missing(virus))
			key.text <- paste("virus:", as.character(pData(fs)$virus[sel]))
		if (missing(cell))
			key.text <- paste("cell:", as.character(pData(fs)$cell[sel]))
		if (missing(moi))
			key.text <- paste("moi:", as.character(pData(fs)$moi[sel]))
	
		overlap(fs, sel, key.text = key.text, ...)
	}
}

# NOTE - this needs to be revised to define the ranges in one channel and then
# perform the MFI in another - this is pretty standard stuff...need to find it in
# flowCore docs...
#
# extract the mfi for channel "chan" in the flowset "fs" defined by breaks
# return matrix of values
#
# Arguments
#	fs		flowSet to be analyzed
#	chan	character specifying channel to analyze
#	breaks	breakpoint(s) to which lower and upper limits will be added
#	FUN		actual function (for example, mean or median)
#	...		additional arguments to go with FUN (na.rm=T, for example)
#	lower.limit	lower limit of ranges
#	upper.limit	upper limit of ranges
#
mfi <- function(fs, chan, breaks, FUN=mean, ..., lower.limit=-Inf, upper.limit=Inf)
{
	if (class(fs) != "flowSet")
		stop("'fs' must be a flowSet")
	if (lower.limit >= upper.limit)
		stop("lower.limit cannot exceed upper.limit")
	if (!is(FUN, "function"))
		stop("FUN must be of class function")

# define and label break points

	breaks <- sort(c(lower.limit, breaks, upper.limit))
	labs <- character()
	for (n in 1:(length(breaks)-1))
		labs <- c(labs, paste("[",breaks[n],"-",breaks[n+1],")",sep=""))

# internal function to perform mean fluorescence measurements

	.mfi <- function(x, lim, mfi.fun=FUN, ...) {
		mfi.fun(x[x >= lim[1] & x < lim[2]], ...)
	}

# accumulate results

	res <- numeric()
	for (n in 1:(length(breaks)-1)) {
		lim <- c(breaks[n], breaks[n+1])
		res <- cbind(res, fsApply(fs, function(x) .mfi(x[, chan], lim), use.exprs = TRUE))
	}
	colnames(res) <- labs
	return(res)
}

#
# generate a 2D-matrix of "cell cycle" distribution values from flowset and breaks
# Each (named) value in 'breaks' identifies the upper limit of the cell cycle stage
#
dnap <- function(fs, breaks, chan = "FL2.A") {
	if (is.null(names(breaks)))
		names(breaks) <- LETTERS[1:length(breaks)]

	res <- numeric()
	for (i in 1:(1+length(breaks)))
		res <- cbind(res, ccp(i, breaks, fs, result = FALSE, chan = chan))
	colnames(res) <- c(names(breaks), paste(">", names(breaks)[length(breaks)], sep=""))
	res
}

#
# cell cycle percentages - extract percentages from flowset defined by endpoint values
# in breaks where 'n' chooses the interval (1, 2, ..., k +1) defined by k breaks
#
ccp <- function(n, breaks, fs, result = TRUE,
		chan = "FL2.A", lower.limit = -Inf, upper.limit = Inf)
{
	if (class(fs) != "flowSet")
		stop("'fs' must be a flowSet")
	if (n > length(breaks) + 1 | n < 1)
		stop("n must be between 1 and", length(breaks))
	if (lower.limit >= upper.limit)
		stop("lower.limit cannot exceed upper.limit")
	breaks <- sort(c(lower.limit, breaks, upper.limit))

	boundary <- list(c(breaks[n], breaks[n+1]))
	names(boundary) <- chan
	selectGate <- rectangleGate(filterId="rectangleGate", boundary)

	fres <- filter(fs, selectGate)
	if (result)
		return(fres)
	else
		return(sapply(fres, function(x) summary(x)$p))
}
