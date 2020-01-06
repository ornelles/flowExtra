.onLoad <- function(libname, pkgname) {
	message("executing: flowViz.par.set(theme = trellis.par.get(), reset = TRUE)") 
	flowViz.par.set(theme = trellis.par.get(), reset = TRUE)
	dev.off()
}
