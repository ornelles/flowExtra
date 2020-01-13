# path to github
	github <- if (.Platform$OS.type == "unix") "~/Documents/github" else "~/github"

# working git
	path <- file.path(github, "flowExtra")
	setwd(path)

# libraries 
	library(flowCore)
	library(flowStats)
	library(flowViz)
	flowViz.par.set(theme = trellis.par.get(), reset = TRUE)

# R files
	ff <- list.files("R/", full = TRUE)
	ff <- setdiff(ff, c("R/my flow code.R", "R/zzz.R"))
	for(f in ff) source(f)

# load and process REP synch data
	fs <- readSet("inst/extdata/synch/")
	bf <- boundaryFilter(c("SSC.H","FSC.H","FL2.H","FL2.A"))
	lg <- linearGate(fs)
	fs <- Subset(fs, bf)
	fs <- Subset(fs, lg)

## 'mixtools' library for mixture models
## works pretty well but forces data to fit given model, need
## to add error/debris component with more explicit modeling

##
## multigaussian function applied to density distribution
##
mfun <- function(x, mean, sigma = NULL, lambda = NULL)
{
# argument checks
	CV <- 0.05
	N <- length(mean)
	if (identical(lambda, NULL))
		lambda <- rep(1, N)
	lambda <- rep(lambda, N)[1:N]
	if (sum(lambda > 1))
		lambda <- lambda/sum(lambda)
	if (identical(sigma, NULL))
		sigma <- abs(CV * mean)
	sigma <- rep(sigma, N)[1:N]

# single gaussian distribution
	gfun <- fun <- function(x, m, s)
		1/(s * sqrt(2 * pi))*exp(-0.5 * ((x - m)/s)^2)

# assemble multiple gaussian function
	xdat <- rep(list(x), N)
	ans <- Map(function(x, mean, sigma, lambda) {
		lambda * gfun(x, mean, sigma)}, xdat, mean, sigma, lambda)

# sum contribution at each x-value
	ans <- do.call(cbind, ans)
	ans <- apply(ans, 1, sum)
	return(ans)
}

#
# works with high quality starting parameters, use other algorithms if need be
# tweak bwFac for best peak discrimination
#
# The better bandwidth selection for these data-rich density plots is NOT
# the default, "nrd0". Rather the method of Sheather and Jones implemented by
# Venables and Ripley and then folded into R >= 3.4.0 is bw.SJ(). The default
# method yields an overly smooth estimate for areas of high density, the SJ
# is best for all data and large amounts of data. 
#
# May need to adjust logic in peakFind()
#
	peaks <- peakFind(fs) 
	par(ask = TRUE)
	for (i in 1:16) {
		d <- density(exprs(fs[[i]][,"FL2.A"]), n = 1024)
		myPeaks <- peaks[i,1:2] # use first two peaks
		sel <- !is.na(myPeaks)
		myPeaks <- myPeaks[sel]
		n <- sum(sel) # number of peaks (currently 2)

	# non-linear fit to multiple Gaussian peaks
		fmx <- nls(y ~ mfun(x, mean, sigma, lambda), data = data.frame(x = d$x, y = d$y),
			start = list(mean = myPeaks, sigma = rep(5, n), lambda = rep(0.9, n)/n))

	# extract and separate coefficients
		cf <- coef(fmx)
		mu <- cf[grep("mean", names(cf))]
		lambda <- cf[grep("lambda", names(cf))]
		sigma <- cf[grep("sigma", names(cf))]

	# visualize with smoothed curve
		plot(d)
		yp <- predict(fmx)
		lines(d$x, yp, col = 2)

	# extract difference, make non-negative
		yd <- ifelse(d$y > yp, d$y - yp, 0)

	# create factor defining areas between peaks
		sf <- 1 # sigma factor to extend from peak
		lo <- c(-Inf, mu - sf * sigma)
		hi <- c(mu + sf * sigma, Inf)
		cuts <- cbind(lo, hi)

	# split d$x and yd according to cuts
		sel <- apply(cuts, 1, function(v) d$x > v[1] & d$x < v[2])
		x.spl <- apply(sel, 2, function(s) d$x[s])
		x.spl <- if (is.matrix(x.spl)) split(x.spl, col(x.spl)) else x.spl
		y.spl <- apply(sel, 2, function(s) yd[s])
		y.spl <- if (is.matrix(y.spl)) split(y.spl, col(y.spl)) else y.spl

	# sum values outside of peaks
		dif <- sapply(y.spl, sum)

	# adjust 'dif' to make area under curve = 1
		adj <- sum(lambda) + sum(dif) - 1
		dif <- dif *(sum(dif) - adj)/sum(dif)

	# add smooth curves for inter-peak stuff for dif > 5%
		fml <- Map(function(x, y) loess(y ~ x, span = 1/3), x.spl, y.spl)

	# add curves for 'dif' values > 5% OR S phase (as here)
	#	sel <- which(dif > 0.05)
		sel <- 2
		for (i in sel) {
			yy <- predict(fml[[i]])
			yy <- ifelse(yy < 0, 0, yy)
			xx <- x.spl[[i]]
			lines(xx, yy, col = 4)
		}
		
	# names...for two peaks
		cc <- c(dif[1], coef(fmx)[5], dif[2], coef(fmx)[6], dif[3])
		names(cc) <- c("<G1", "G1", "S phase", "G2/M", ">G2")
		txt <- c(names(cc), sprintf("%3.1f%%", 100*cc))
		legend("topright", legend = txt, xjust = 1, ncol = 2)
		print(cc)
	}

## Works again with "improved" peakFind
#
# folded into function to convert FL2.A values to density, fit, return lambda
# values and "inter-peak" fraction

