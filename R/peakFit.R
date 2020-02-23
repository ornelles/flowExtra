#' Fit G1 and G2 peaks
#' 
#' Fit an arbitrary number of Gaussian peaks.
#'
#' @md
#' @details
#' Currently a stub for the `peakFit` function.
#'
#' @param x Argument check
#'
#' @return Echo argument
#'
#' @export
#'
peakFit <- function(x) {
	cat ("peakFit:", x, "\n")
	invisible(x)
}
#' @name mgfun
#' @rdname peakFit
#'
#' @title
#' Multiple Gaussian density distribution
#'
#' @description
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
#' @export
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
#' @name gfun
#' @rdname peakFit
#'
#' @title
#' Gaussian density distribution
#'
#' @description
#' This is a simpler version of [stats::dnorm] and runs a little faster.
#'
#' @param x Numeric values at which to evaluate the function
#' @param mean Mean value for the Gaussian distribution
#' @param sd Numeric value of standard deviation (sigma) for Gaussian
#'   function
#'
#' @return The Gaussian density for given `x` values
#'
#' @export
#'
gfun <- function(x, mean, sd) {
	1/(sd * sqrt(2 * pi))*exp(-0.5 * ((x - mean)/sd)^2)
}
