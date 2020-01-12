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
		
# explore
	d <- density(exprs(fs[[16]][,"FL2.A"]), adj = 0.5, n = 1024)
	peaks <- peakFind(fs[[16]])[1:2]
	start <- list(a = 0.5, b = 0.5, m1 = peaks[1], s1 = 10, m2 = peaks[2], s2 = 10)
	fm <- nls(y ~ fun(x, a, b, m1, s1, m2, s2), data.frame(x = d$x, y = d$y), 
		start = start, algorithm = "port")
	coef(fm)
	# "s-phase" for clean profile
	1 - sum(coef(fm)[1:2])
	plot(d)
	lines(d$x, predict(fm), col = 2)
# change coefficients to vectors and you got multinomal fitting

# 'mixtools' library for mixture models
# works pretty well but forces data to fit given model, need
# to add error/debris component with more explicit modeling

##
## multigaussian function applied to density distribution
##
fun <- function(x, mean, sigma = NULL, lambda = NULL, cv = 0.05)
{
# argument checks
	N <- length(mean)
	if (N < 1 | N > length(x) / 10)
		stop("length 'mean' should be between 1 and ", round(length(x)/10))
	if (identical(lambda, NULL))
		lambda <- rep(1, N)/N
	if (sum(lambda > 1))
		lambda <- lambda/sum(lambda)
	if (length(lambda) != N)
		stop("length 'lambda' should be ", N)
	if (cv < 0 | cv > 2)
		warning("'cv' should be between 0 and 1")
	if (identical(sigma, NULL))
		sigma <- abs(cv * mean)
	if (any(sigma < 0))
		stop("'sigma' should be greater than 0")
	sigma <- rep(sigma, N)[1:N] # replicate along peaks

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
	peaks <- peakFind(fs)[,1:2]
	par(ask = TRUE)
	for (i in 1:16) {
		d <- density(exprs(fs[[i]][,"FL2.A"]), adj = 0.5, n = 1024)
		myPeaks <- peaks[i,]
		sel <- !is.na(myPeaks)
		myPeaks <- myPeaks[sel]
		n <- sum(sel)
		fmx <- nls(y ~ fun(x, mean, sigma, lambda), data = data.frame(x = d$x, y = d$y),
			start = list(mean = myPeaks, sigma = rep(5, n), lambda = rep(1, n)/n))
		plot(d)
		lines(d$x, predict(fmx), col = 2)
	}

## Once worked -- now fails with "improved" peakFind...
#
# fold into function to convert FL2.A values to density, fit, return labda
# values and "inter-peak" fraction

