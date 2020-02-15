#' Fit G1 and G2 peaks
#' 
#' Fit an arbitrary number of Gaussian peaks 
#'
#' @md
#' @description 
#'

#' @name gfun
#' @rdname peakFit
#'
#' Density distribution function for a single Gaussian.
#' This is a simpler version of [stats::dnorm] and runs a little faster.
#'
#' @return The Gaussian density for given `x` values
#'
	gfun <- function(x, mean, sd)
		1/(sd * sqrt(2 * pi))*exp(-0.5 * ((x - meam)/sd)^2)

#' Replacement Gaussian dnorm function
#'
#' @name mgfun
#' @rdname peakFit
#'
#' Density distribution function for multiple Gaussian function of 
#' given means (`mean`) and standard deviations (`sd`), each adjusted  by
#' the multiplier `lambda` where `sum(lambda) = 1`.
#'
#' @param x Values at which to evaluate the function
#' @param mean Numeric vector of mean for each Gaussian component
#' @param sd Numeric vector of standard deviation (sigma) for each Gaussian
#'   component with a default value of `0.05 * mean`
#' @param lambda Numeric vector serving as the multiplier for
#'   the relative area of each Gaussian component with a default value of 
#'   `1/length(mean)` for each. The value will be adjusted to
#'   `lambda/sum(lambda)` to permit a value such as `c(2, 5)` to express the
#'   relative **area** of each peak.
#'
#' @return The combined Gaussian density for given `x` values and parameters
#'
	mgfun <- function(x, mean, sd = 0.05 * mean, lambda = 1/length(mean))
	{
	# argument checks
		N <- length(mean)
		sd <- rep(sd, N)[seq_len(N)] # ensure number matches mean
		lambda <- rep(lambda, N)[seq_len(N)] # ensure number matches mean
		if (sum(lambda > 1))
			lambda <- lambda/sum(lambda)

	# apply single Gaussian for each value of mean/sd/lambda
		ans <- Map(function(x, m, s, k) {
			k * gfun(x, m, s)}, rep(list(x), N), mean, sd, lambda)

	# sum contribution from each Gaussian
		ans <- do.call(cbind, ans)
		ans <- apply(ans, 1, sum)
		return(ans)
	}
