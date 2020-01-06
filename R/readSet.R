#' Read and Transform Flow Set
#' 
#' A wrapper to \code{\link{read.flowSet}} to read a flow set.
#' 
#' @param path A character vector identifying the path to search for FSC
#'   files, defaults to current directory
#' @param pattern Grep pattern to identify FSC files, defaults to
#'   \code{"\\d\\d\\d$"}
#' @param log.chan Character vector of channels to be log-transformed (after
#'  the \code{"alter.names"} option is applied)
#' @param alter.names Logical value passed to \code{read.flowSet} to 
#'   change "-" to "." in names, default value of \code{TRUE}
#' @param name.keyword Character string for the  FCS keyword used for names,
#'   default of \code{"SAMPLE ID"}
#' @param phenoData List to serve as the rudiments of the \code{phenoData} 
#'    built from the keywords \code{"$DATE"} and \code{"SAMPLE ID"}
#' @param tfun Transformation function to apply to log-transformed values,
#'   either \code{asinh} (default) or \code{log10}
#' @param depth Integer indicating bit depth (\code{2^depth}), 10 for
#'   FACSCaliber, 24 for Acurri
#' @param log.cutoff Values greater than this cutoff force transformation if 
#'   the argument \code{log.chan} is \code{NULL}, defaults to \code{2^depth}
#' @param ...	Additional arguments to \code{\link[flowCore]{read.flowSet}}
#'
#' @details
#' Files in the directory specified by \code{path} that match 
#' the \code{grep} pattern specified in \code{pattern} will be read with
#' \code{\link{read.flowSet}}. Parameter names will be adjusted if
#' \code{alter.names = TRUE} to replace \code{"-"} with \code{"."}. A 
#' minimal amount of \code{phenoData} will be assigned from the $DATE
#' and "SAMPLE ID" keywords.
#'
#' @return
#' 
#' Flow Set with transformed and name-adjusted parameters.
#' 
#' @import flowCore
#'
#' @export
#' 
readSet <- function(path = ".", pattern = "\\d\\d\\d$", log.chan = NULL,
	alter.names = TRUE, name.keyword = "SAMPLE ID",
	phenoData = list(date = "$DATE", sample = "SAMPLE ID"),
	tfun = c("asinh", "log10"), depth = 10, log.cutoff = 2^depth, ...) 
{
	tfun <- match.arg(tfun)
	if (missing(path))
		path <- dirname(file.choose())
	if (any(grepl("fcs$", list.files(path), ignore.case = TRUE))) {
		depth <- 24
		pattern <- "fcs$"
	}
	else if (any(grepl(pattern, list.files(path), ignore.case = TRUE))) {
		depth <- depth
		log.cutoff <- 2^depth
	}
	else
		stop("unable to find suitable FACS data file with pattern: ", pattern)

	fs <- flowCore::read.flowSet(path = path,
			pattern = pattern,
			alter.names = alter.names, 
			phenoData = phenoData, ...)
	key.names <- names(flowCore::keyword(fs[[1]]))

# adjusted selection key to match revision of flowCore ?
	channels <- flowCore::colnames(fs[[1]])

	if (is.null(log.chan)) { # transform based on value
		sel <- grep("P\\dRmax", key.names, value = TRUE)
		log.chan <- channels[which(keyword(fs[[1]], sel) > log.cutoff)]
	}
	else if (is.character(log.chan)) { # transform named channels
		log.chan <- channels[channels %in% log.chan]
		if (length(log.chan) == 0)
			stop("unable to identify any 'log.chan' values in data")
	}
	else if (log.chan == TRUE) # transform ALL channels
		log.chan <- channels

	if (length(log.chan) > 0) {
		tl <- flowCore::transformList(log.chan, tfun)
		fs <- flowCore::transform(fs, tl)
	}
	return(fs)
}
