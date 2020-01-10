	f <- function(x, A1, mod1, A2, mod2) {
		m1 <- mod1$estimate[1]
		s1 <- mod1$estimate[2]
		m2 <- mod2$estimate[1]
		s2 <- mod2$estimate[2]
		K <- sqrt(2 * pi) 
		p1 <- A1/(s1*K)*exp(-0.5*((x - m1)/s1)^2)
		p2 <- A2/(s2*K)*exp(-0.5*((x - m2)/s2)^2)
		ans <- p1 + p2

f <- function(x, a, b, m1, s1, m2, s2) {
	K <- sqrt(2 * pi)
	p1 <- a/(s1*K)*exp(-0.5*((x-m1)/s1)^2)
	p2 <- b/(s2*K)*exp(-0.5*((x-m2)/s2)^2)
	ans <- p1 + p2
}

d <- density(exprs(fs[[16]][,"FL2.A"]), adj = 0.5, n = 1024)
pks <- peakFind(fs[[16]])[1:2]
start <- list(a = 0.5, b = 0.5, m1 = pks[1], s1 = 10, m2 = pks[2], s2 = 10)
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
#
fmx <- nls(y ~ fun(x, mean, sigma, lambda), data = data.frame(x = d$x, y = d$y),
	start = list(mean = pks[1,], sigma = c(5, 8), lambda = c(0.7,0.3)))

## Works!
#
# fold into function to convert FL2.A values to density, fit, return labda
# values and "inter-peak" fraction
