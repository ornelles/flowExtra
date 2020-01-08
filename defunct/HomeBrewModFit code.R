# Homebrewed ModFit
	f <- function(x, A1, mod1, A2, mod2) {
		m1 <- mod1$estimate[1]
		s1 <- mod1$estimate[2]
		m2 <- mod2$estimate[1]
		s2 <- mod2$estimate[2]
		K <- sqrt(2 * pi) 
		p1 <- A1/(s1*K)*exp(-0.5*((x - m1)/s1)^2)
		p2 <- A2/(s2*K)*exp(-0.5*((x - m2)/s2)^2)
		ans <- p1 + p2
		return(ans)
	}

	D <- density(x)

	x.g1 <- x[x>50 & x < 220]
	d <- density(x.g1)
	xmid <- d$x[which.max(d$y)]
	x1 <- x.g1[x.g1 <= xmid]
	xx <- c(x1, 2*xmid - x1)
	mod1 <- MASS::fitdistr(xx, "normal")
	m1 <- mod1$estimate[1]
	s1 <- mod1$estimate[2]
	A <- sum(xx > m1 - 2.5*s1 & xx < m1 + 2.5*s1)

	x.g2 <- x[x>310 & x < 400]
	d <- density(x.g2)
	xmid <- d$x[which.max(d$y)]
	x2 <- x.g2[x.g2 >= xmid]
	xx <- c(x2, 2*xmid - x2)
	mod2 <- MASS::fitdistr(xx, "normal")
	m2 <- mod2$estimate[1]
	s2 <- mod2$estimate[2]
	B <- sum(xx > m2 - 2.5*s2 & xx < m2 + 2.5*s2)

	xp <- 100:500
	yp <- f(xp, A, fit.g1, B, fit.g2)

	plot(D$x, (A+B)*D$y, type = "l")
	lines(xp, yp, col = 2)

	total <- length(x)
	g1 <- A/total
	s <- (total - A - B)/total
	g2 <- B/total
	cbind(g1, s, g2)