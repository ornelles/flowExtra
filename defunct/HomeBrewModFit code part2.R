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
## single gaussian distribution
##
	gfun <- fun <- function(x, m, s)
		1/(s * sqrt(2 * pi))*exp(-0.5 * ((x - m)/s)^2)

##
## multigaussian function applied to density distribution
##
	mgfun <- function(x, mean, sigma = NULL, lambda = NULL)
	{
	# argument checks
		CV <- 0.05
		N <- length(mean)
		if (identical(lambda, NULL))
			lambda <- rep(1, N)
		lambda <- rep(lambda, N)[seq_len(N)]
		if (sum(lambda > 1))
			lambda <- lambda/sum(lambda)
		if (identical(sigma, NULL))
			sigma <- abs(CV * mean)
		sigma <- rep(sigma, N)[seq_len(N)]

	# assemble multiple gaussian function
		xdat <- rep(list(x), N)
		ans <- Map(function(x, mean, sigma, lambda) {
			lambda * gfun(x, mean, sigma)}, xdat, mean, sigma, lambda)

	# sum contribution at each x-value
		ans <- do.call(cbind, ans)
		ans <- apply(ans, 1, sum)
		return(ans)
	}
##
## left/right half peak fitting
##
# Possibly asymmetric Gaussian fit for left or right side of gaussian
# Take difference and use var = (actual - fit)^2 for variance for weighting
#
	gfit <- function(x, limit = 0.1) {
		d <- density(x)
		xmid <- d$x[which.max(d$y)]
		skew <- (mean(x) - xmid)/sd(x)
		if (skew > limit) { # right tail
			tail <- "right"
			xl <- x[x <= xmid]
			xx <- c(xl, 2*xmid  - xl)
		}
		else if (skew < -limit) { # left tail
			tail <- "left"
			xr <- x[x >= xmid]
			xx <- c(2*xmid - xr, xr)
		}
		else {
			tail <- "none"
			xx <- x
		}
		ans <- MASS::fitdistr(xx, "normal")
		cf <- ans$estimate
		yf <- dnorm(d$x, cf[["mean"]], cf[["sd"]])
		scale <- max(d$y)/max(yf)
		return(list(mean = cf[["mean"]], sd = cf[["sd"]], scale = scale,
			skew = skew, tail = tail))
	}

# generate skewed distribution with mgfun
	xd <- seq(1, 100, len = 1000)
	yd <- mgfun(xd, c(40, 50, 60, 70), c(10, 20, 20, 20), lambda = c(12, 6, 4, 2))

# generate 1e4 X values for distribution and Z for reversed distribution
	set.seed(321)
	n <- round(2e4*(yd + rnorm(yd, mean = 0, sd = mean(yd)/20)))
	n <- ifelse(n < 0, 0, n)
	X <- sample(rep(xd, n), 1e4)

	n <- round(2e4*(yd + rnorm(yd, mean = 0, sd = mean(yd)/20)))
	n <- ifelse(n < 0, 0, n)
	Y <- sample(rep(rev(xd), n), 1e4)

	yn <- dnorm(xd, 40, 10)
	n <- round(2e4*(yn + rnorm(yn, mean = 0, sd = mean(yd)/20)))
	n <- ifelse(n < 0, 0, n)
	Z <- sample(rep(xd, n), 1e4)
	
# Show that it works
	plot(density(X))
	ans <- gfit(X)
	lines(xd, dnorm(xd, ans$mean, ans$sd)*ans$scale, col = 2)

	plot(density(Y))
	ans <- gfit(Y)
	lines(xd, dnorm(xd, ans$mean, ans$sd)*ans$scale, col = 2)

# Take difference between actual and symmetric/scaled result from gfit as
# the variance to form weights as 1/variance^2
#
# Decide on how to form (possibly overlapping) weights to pass to nls with
# in mgfun (below) 
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

##
## Figure out how to weight with asymmetric fits
##
	# non-linear fit to multiple Gaussian peaks
		fmx <- nls(y ~ mgfun(x, mean, sigma, lambda), data = data.frame(x = d$x, y = d$y),
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
		names(cc) <- c("<G1", "G1", "S phase", "G2", ">G2")
		txt <- c(names(cc), sprintf("%3.1f%%", 100*cc))
		legend("topright", legend = txt, xjust = 1, ncol = 2)
		print(cc)
	}

## Works again with "improved" peakFind
#
# folded into function to convert FL2.A values to density, fit, return lambda
#
# values and "inter-peak" fraction
##
## parameters to capture for peakfit object?
##
# peak parameters (mean, sigma, lambda)
# other parameters (loess model with model frame?)
# channel


##
## parameters to report for peakfit object?
##
# File
# G1 % at peak
# G2 % at peak
# G2/G1 ratio
# S
# %CV (of G1 peak)
# total events
# events per channel

# multi-peak generating cell cycle faking!
	xp <- 1:1024
	mu <- seq(200, 400, len = 9)
	lam <- rep(0, length(mu))
	lam[c(1,9)] <- c(70, 15)
	lam[2:8] <- c(5,3,1,1,1,1,3)
	yp <- mgfun(xp, mu, lambda = lam)
	y <- round(2e4*yp)
	err <- sample(c(0,0,0,1,1,1,2,2,3), length(xp), TRUE)
	C <- rep(xp, y + err)
	C <- sample(C, 1e4)
	plot(density(C), ylim = c(0, 0.03))

	x1 <- subset(xp, xp > 100 & xp < 300)
	y1 <- f1$scale * gfun(xp, f1$mean, f1$sd)
	f1 <- gfit(C[C > 100 & C < 300])

	x1 <- subset(xp, xp > 300 & xp < 500)
	f2 <- gfit(C[C > 300 & C < 500])
	y2 <- f2$scale * gfun(xp, f2$mean, f2$sd)

# Not fitting, but could be good for extracting variance
	plot(density(C), ylim = range(y1, y2))
	lines(xp, y1 + y2, col = 2)

