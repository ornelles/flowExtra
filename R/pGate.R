#' Create Arbitrary Gate from Lattice Panel
#' 
#' Interact with lattice flow plot to define a gate
#' 
#' @param type Type of gate as a character string (see details)
#' @param vertical,horizontal Logical values to impose vertical or
#'   horizontal boundaries for \code{type = "parallel"} gate
#'   (default values of \code{vertical = TRUE} and
#'   \code{horizontal = FALSE})
#' @param filterId A character vector to identify this filter.
#' @param adjust.ellipse.distance Fudge factor for
#'   \code{\link{ellipsoidGate}} (default of \code{sqrt(2)})
#' 
#' @details
#' This function is called after an \code{xyplot} of flow data in order
#' to define a gate. 
#' The user will select a panel that will be redrawn and used to select
#' points that define the gate for the data in the remaining panels
#' (`flowSet`).
#' The type of gate is specific by the argument \code{type} and can be
#' one of the following character strings. Only the first letter is
#' required since partial matching is used to assign the value.
#' 
#' \tabular{ll}{
#' 	\code{"arbitrary"} \tab The selected points will define the vertices of an 
#' 		arbitrary polygon. \cr
#' 	\code{"parallel"} \tab A parallelogram specified by 3 points with optional 
#' 		vertical and horizontal constraints specified by \code{vertical} and
#' 		\code{horizontal}. \cr
#' 	\code{"convex"} \tab The convex hull of the selected points. \cr
#' 	\code{"ellipse"} \tab The ellipsoid hull (minimal ellipse) encompassing the
#' 		selected points. \cr
#' }
#' 
#' The implementation of \code{\link{ellipsoidGate}} does not seem to honor 
#' the covariance parameter, so the distance parameter has been arbitrarily  
#' adjusted by the square root of 2 to make up for this quirk.
#' 
#' @return
#' A \code{polygonGate} or \code{ellipsoidGate} object is returned of the  
#' same dimensions as the active lattice graph.
#' 
#' @import flowCore cluster
#'
#' @export
#'
pGate <- function(type = c("arbitrary", "parallel", "convex", "ellipse"),
			vertical = TRUE, horizontal = FALSE, filterId  = NULL,
			adjust.ellipse.distance = sqrt(2)) {

	if (current.row() == 0) {
		if (length(trellis.currentLayout()) > 1)
			cat("Select panel for focus...\n")
		fp <- do.call(trellis.focus, list())
	}
	if (length(trellis.currentLayout()) != 1) {
		cat("Redrawing lattice object...\n")
		plot(trellis.last.object()[which.packet()])
		do.call(trellis.focus, list())
	}

	type <- match.arg(type)

	if (type == "arbitrary" | type == "convex" | type == "ellipse") {
		p <- matrix(numeric(0), ncol = 2)
		while(!is.null(gp <- grid::grid.locator())) {
			do.call(panel.points, args = list(gp, cex = 0.5))
			p <- rbind(p, c(grid::convertX(gp$x, "native", TRUE),
					grid::convertY(gp$y, "native", TRUE)))
		}
		if (type == "convex") {
			idx <- chull(p)
			bd <- p[idx,]
		}
		else if (type == "ellipse") {
			exy <- cluster::ellipsoidhull(p)
			bd <- cluster::predict.ellipsoid(exy, n.out = 51)
		}
		else
			bd <- p
		do.call(panel.lines, args = list(rbind(bd, bd[1,])))

		if (is.null(filterId)) {
			if (type == "arbitrary")
				filterId <- "defaultPolygonGate"
			else if (type == "ellipse")
				filterId <- "defaultEllipsoidGate"
			else
				filterId <- "defaultConvexHull"
		}
	}
	else if (type == "parallel") {
		delta <- function(p1, p2) {
			x <- grid::convertX(p2$x, "npc", TRUE) - grid::convertX(p1$x, "npc", TRUE)
			y <- grid::convertY(p2$y, "npc", TRUE) - grid::convertY(p1$y, "npc", TRUE)
			list(x = x, y = y)
		}
	
		p1 <- grid::grid.locator()
		x1 <- as.numeric(p1$x); y1 <- as.numeric(p1$y)
			do.call(panel.points, args = list(p1, cex=0.5))
	
		p2 <- grid::grid.locator()
		x2 <- as.numeric(p2$x); y2 <- as.numeric(p2$y)
		d <- delta(p1, p2)
	
		if (abs(d$x) >= abs(d$y) & horizontal)
			y2 <- y1
		if (abs(d$y) >= abs(d$x) & vertical)
			x2 <- x1
		do.call(panel.points, args = list(x = x2, y = y2, cex=0.5))
		do.call(panel.lines, args = list(cbind(x=c(x1,x2), y=c(y1,y2))))
	
		p3 <- grid::grid.locator()
		x3 <- as.numeric(p3$x); y3 <- as.numeric(p3$y)
		d <- delta(p3, p2)
		if (abs(d$x) >= abs(d$y) & horizontal)
			y3 <- y2
		if (abs(d$y) >= abs(d$x) & vertical)
			x3 <- x2
		do.call(panel.points, args = list(x = x3, y = y3, cex=0.5))
		do.call(panel.lines, args = list(cbind(x=c(x2,x3), y=c(y2,y3))))
	
		x4 <- y4 <- numeric(0)
	
		if (vertical) {
			if (abs(x3 - x2) < abs(x2 - x1))
				x4 <- x1
			else
				x4 <- x3
		}
		if (horizontal) {
			if (abs(y3 - y2) < abs(y2 - y1))
				y4 <- y1
			else
				y4 <- y3
		}
		if (length(x4) == 0)
			x4 <- x1 + (x3 - x2)
		if (length(y4) == 0)
			y4 <- y1 + (y3 - y2)
	
		do.call(panel.points, args = list(x = x4, y = y4, cex=0.5))
		do.call(panel.lines, args = list(cbind(x=c(x3,x4), y=c(y3,y4))))
		do.call(panel.lines, args = list(cbind(x=c(x4,x1), y=c(y4,y1))))
	
		bd <- cbind(x=c(x1,x2,x3,x4), y=c(y1,y2,y3,y4))
		ord <- order(bd[,1], bd[,2])
		bd <- bd[ord[c(1,3,4,2)],]

		if (is.null(filterId))
			filterId <- "defaultParallelGate"
	}
	else
		stop("unexpected value for 'type' = ", type)

	val <- do.call(trellis.panelArgs, list())[c("channel.x.name","channel.y.name")]
	colnames(bd) <- unlist(val)

	if (type == "ellipse") {
		cov <- exy$cov
		rownames(cov) <- colnames(bd)
		colnames(cov) <- colnames(bd) 
		mean <- exy$loc
		names(mean) <- colnames(bd)
		return(ellipsoidGate(.gate = cov, mean = exy$loc, distance = adjust.ellipse.distance))
	}
	else
		return(polygonGate(bd, filterId = filterId))
}
