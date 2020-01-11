.onLoad <- function(libname, pkgname) {
	message("executing: flowViz.par.set(theme = trellis.par.get(), reset = TRUE)") 
	val <- dev.cur()
	flowViz.par.set(theme = trellis.par.get(), reset = TRUE)
	if (val == 1)
		dev.off()
}
