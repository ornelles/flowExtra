#' Locator Function for Lattice
#' 
#' Reads the position of the graphics cursor from a lattice panel
#' 
#' @details This provides a functional replacement for
#' \code{\link{locator}} that works with \code{lattice} objects. A
#' \code{lattice} panel will be selected, redrawn as a single panel and
#' used to interact with the user. The position of the graphics cursor
#' will be recorded with each click of the (first) mouse button.
#' 
#' By default, a point will be added at that position unless otherwise
#' specified by the argument \code{type}.
#' 
#' Like the \code{\link{locator}} function, this function returns 
#' values in \code{"usr"} coordinates as a list.
#' 
#' @param n The maximum number of points to locate
#' @param type One of \code{"dot", "v", "h", "both"} or 
#'   \code{"none"} to add a point, vertical line, horizontal line,
#'   both vertical and horizontal lines, or no mark at each selected
#'   point
#' @param ... Additional graphic parameters passed to either
#'   \code{\link{panel.points}} or \code{\link{panel.abline}} if
#'   \code{type != "none"}
#' 
#' @return
#' A list containing \code{x} and \code{y} coordinates of the
#' identified points in the user coordinate system.
#' 
#' @import lattice
#' 
#' @export
#' 
locator2 <- function(n = 512, type = c("dot", "v", "h", "both", "none"),
		...)
{
	type <- match.arg(type)

	if (current.row() == 0) {
		if (length(trellis.currentLayout()) > 1) {
			cat("Select panel for focus...\n")
			flush.console()
		}
		do.call(trellis.focus, list())
	}
	if (length(trellis.currentLayout()) != 1) {
		cat("Redrawing lattice object...\n")
		flush.console()
		plot(trellis.last.object()[which.packet()])
		do.call(trellis.focus, list())
	}
# loop until interrupted or n == 0
	x <- y <- numeric(0) # collected coordinates
	while (n > 0 && !is.null(pp <- grid::grid.locator())) {
		xc <- as.numeric(pp$x) # coordinates come with unit character
		yc <- as.numeric(pp$y)
		if (type == "v" | type == "both")
			do.call(panel.abline, args = list(v = xc, ...))
		if (type == "h" | type == "both")
			do.call(panel.abline, args = list(h = yc, ...))
		if (type == "dot")
			do.call(panel.points, args = list(x = xc, y = yc, ...))
		x <- c(x, xc)
		y <- c(y, yc)
		n <- n - 1
	}
# transform to 'usr' coordinates
	logX <- trellis.last.object()$x.scales$log
	logY <- trellis.last.object()$y.scales$log
	x <- if (logX) logX^x else x
	y <- if (logY) logY^y else y

	if (length(x))
		return(list(x = x, y = y))
}
