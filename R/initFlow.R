#
# initFlow 
#
initFlow <- function() {
	library(MASS)
	library(flowCore)
	library(flowStats)
	library(flowViz)

	env <- attach(NULL, name = "flow")
	dropbox <- if (.Platform$OS.type == "unix") "~/Dropbox" else "~/../Dropbox"
	f <- file.path(dropbox, "/docs/flow/flow functions/my flow code.R")

	if (file.exists(f))
			sys.source(f, env)
	else
		stop("Unable to find 'my flow code.R' in Dropbox\n")

	cat("Attached 'flow'\n")
	invisible(env)
}
