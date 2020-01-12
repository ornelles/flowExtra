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

# load and process brdu data
	fs <- readSet("inst/extdata/brdu/")
	xyplot(SSC.H ~ FSC.H, fs, nbin = 256)
	n2f <- norm2Filter(c("SSC.H", "FSC.H"), scale = 2.5)
	xyplot(SSC.H ~ FSC.H, fs, filter = n2f, stats = TRUE, nbin = 256)
	fs <- Subset(fs, n2f)
	lg <- linearGate(fs, "FL3.A", "FL3.H", zero = FALSE, gRange = FALSE)
	xyplot(FL3.H ~ FL3.A, fs, filter = lg, nbin = 256, stats = TRUE)
	fs <- Subset(fs, lg)

# adjust FL3.H data DNA/7-AAD
	dnaplot(fs, chan = "FL3.H")
	fs <- peakAdjust(fs, "FL3.H")
	dnaplot(fs, chan = "FL3.H")

# adjust FL1.H data BrdU/FITC
	dnaplot(fs, chan = "FL1.H")
	fs <- warpSet(fs, "FL1.H")
	dnaplot(fs, chan = "FL1.H")

# assign gates from data in Data.002
	xyplot(FL1.H ~ FL3.H, fs, nbin = 256)
	s <- pGate("par")
	g2 <- pGate("par")
	g1 <- pGate("par")

# assemble list of gates and show results
	gates <- filters(list(G1 = g1, S = s, G2 = g2))
	flist <- setNames(rep(list(gates), 3), pData(fs)$name)	
	xyplot(FL1.H ~ FL3.H, fs, filter = flist, nbin = 256, stats = TRUE)

# extract the results with getFres and plot
	res <- sapply(gates, function(g) getFres(fs, g))
	res
	barchart(res, horizontal = FALSE, auto.key = T)