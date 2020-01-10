	f <- function(x, A1, mod1, A2, mod2) {
		m1 <- mod1$estimate[1]
		s1 <- mod1$estimate[2]
		m2 <- mod2$estimate[1]
		s2 <- mod2$estimate[2]
		K <- sqrt(2 * pi) 
		p1 <- A1/(s1*K)*exp(-0.5*((x - m1)/s1)^2)
		p2 <- A2/(s2*K)*exp(-0.5*((x - m2)/s2)^2)
		ans <- p1 + p2

fun <- function(x, a, b, m1, s1, m2, s2) {
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