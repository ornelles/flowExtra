#' Read and Transform Flow Set
#' 
#' Read a flow set and transform values with maximum values above given cutoff.
#' 
#' @param path A character vector identifying the path to search for the FSC files, 
#'  defaults to current directory
#' @param pattern Grep pattern to identify fsc files, defaults to \code{"\\d\\d\\d$"}
#'  log.chan Character vector of channels to be log-transformed (after 
#'  "alter.names")
#' @param alter.names Logical value passed to \code{read.flowSet} to change "-" to 
#'  "." in names
#' @param name.keyword Character string for the  FCS keyword used for names, default 
#'  "SAMPLE ID"
#' @param  phenoData List to serve as the rudiments of the phenoData built from 
#'  keywords
#' @param  tfun transformation function to apply to log-transformed values, either \
#'  code{asinh} or \code{log10}
#' @param depth Integer indicating bit depth (2^depth), 10 for FACVSCaliber, 24 for 
#'  Acurri
#' @param log.cutoff Values greater than this value force log-transfromation if 
#'  log.chan is \code{NULL}, defaults to \code{2^depth}
#' @param ...	Additional arguments to \code{\link{read.Flowset}}
#' 
#'
#' @return
#' 
#' Flow Set with transformed and name-adjusted parameters.
#'
#' @export
#'
readSet <- function(path = ".", pattern = "\\d\\d\\d$", log.chan = NULL,
	alter.names = TRUE, name.keyword = "SAMPLE ID", phenoData = list(date = "$DATE",
	sample = "SAMPLE ID"), tfun = c("asinh", "log10"), depth = 10,
	log.cutoff = 2^depth, ...) 
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

	fs <- read.flowSet(path = path,
			pattern = pattern,
			alter.names = alter.names, 
			phenoData = phenoData, ...)
	key.names <- names(keyword(fs[[1]]))

# adjusted selection key to match revision of flowCore ?
	channels <- colnames(fs[[1]])

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
		tl <- transformList(log.chan, tfun)
		fs <- transform(fs, tl)
	}
	return(fs)
}
