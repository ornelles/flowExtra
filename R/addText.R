#' Add text to an existing lattice plot
#'
#' This serves the same role as the base function [graphics::text()] for
#' use in `lattice` plots.  (Also can be called as `text2()`.)
#'
#' @md
#' @param x,y Numeric vectors of coordinates where text `labels` should
#'  be written. If `y` is missing, `x` is treated as a list with 'x' and 'y'
#' @param labels Character vector or expression specifying the text to
#'  be written. If `x` is specified as a list and `y` is a character
#'  vectors, `labels` will use the values in `y`
#' @param col,cex,adj Default values passed to [lattice::panel.text()]
#' @param ... Other parameters are passed to [lattice::panel.text()] 
#' 
#' @details
#' 
#' Text labels can be added to an existing simple lattice plot with this
#' function. Values in `x`, `y` and `labels` are recycled as needed. 
#' Interaction with the existing plot occurs through the function
#' [lattice::trellis.currentLayout()]. See that function for more details. 
#'
#' @return
#' Nothing is returned. This function is called for the side effect of
#' adding labels. 
#' 
#' @import
#' lattice
#'
#' @name addText
#' @aliases text2
#'
#' @examples
#'
#' @export
#' 
addText <- function(x, y, labels, col = "red", cex = 2/3, adj = 0.5, ...) {

	mat <- trellis.currentLayout("panel")
	N <- sum(mat != 0)
	if (N == 0)
		stop ("no lattice object found")

	if (missing(y)) {			# x should be a list
		if (all(c("x", "y") %in% names(x))) {
			xx <- x$x
			yy <- x$y
		}
	}
	else if (is.character(y)) {	# y is assumed to have labels
		labels <- y
		if (all(c("x", "y") %in% names(x))) {
			xx <- x$x
			yy <- x$y
		}
	}
	else {		# x and y better be numeric...
		xx <- as.numeric(x)
		yy <- as.numeric(y)
	}
	if (missing(labels))
		labels <- seq_along(xx)

	if (length(xx) == 1) {
		xx <- rep(xx, N)
		yy <- rep(yy, N)
	}
	if (length(xx) != N)
		stop("number of panels (", N, ") and number of points (", length(xx), ") differ")

	for (i in 1:N) {
		column <- which(apply(mat, 2, function(x) any(x == i)))
		row <- which(apply(mat, 1, function(x) any(x == i)))
		trellis.focus("panel", column = column, row = row, highlight = FALSE)
		do.call(panel.text, list(xx[i], yy[i], labels[i], cex = cex, col = col, adj = adj, ...))
		trellis.unfocus()
	}
}

# alias
text2 <- addText
