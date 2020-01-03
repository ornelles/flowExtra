#' Extract a Filter Result
#' 
#' Extract a single summary statistic from a \code{filterResult} or
#' \code{flowFrame/flowSet} plus \code{filter}
#' 
#' This function will extract a selected statistic from a filter
#' result. If needed, it will apply a \code{filter} to a
#' \code{flowFrame} or \code{flowSet} and then extract the desired
#' statistic. The positive or negative population can be specified.
#' 
#' @details
#' This is a convenience function to extract information normally produced
#' by sequentially calling \code{filter}, \code{summary} and an extraction
#' function to retrieve slot values. The extracted statistic is either the
#' positive fraction (\code{p}), positive number (\code{true}) or
#' total number (\code{count}).
#' 
#' If the first argument (\code{x}) is a \code{filterResult} or \code{
#' filterResultList}, the argument \code{filter} is ignored \emph{unless}
#' it is a single character string that will be interpreted as \code{what}.
#' This allows for lazy applications of the function to filter results
#' such as \code{getFres(res, "n")} to extract the total number of
#' events from the object \code{res}. If the first argument is 
#' a \code{flowFrame} or \code{flowSet}, \code{filter} must be
#' provided in order to apply it to the data before extracting
#' the statistic.
#' 
#' The specific statistic can be specified by the character value in
#' \code{what} as indicated below:
#' \tabular{ll}{
#' 	\code{"p"} \tab The positive fraction of events. \cr
#' 	\code{"true"} \tab Number of positive events. \cr
#' 	\code{"positive"} \tab Number of positive events. \cr
#' 	\code{"negative"} \tab Number of negative events. \cr
#' 	\code{"n"} \tab Total number of events. \cr
#' 	\code{"count"} \tab Total number of events. \cr
#' }
#' Note that with \code{what = "n"} or \code{what = "count"},
#' the \code{excluded} option is ignored.
#' 
#' @param x Either a \code{filterResult} or a \code{flowFrame}
#'   or a \code{flowSet}
#' @param filter An optional \code{filter}, required if \code{x} is not a
#'   \code{filterResult}
#' @param what A single character string specifying the statistic to
#'   return. The default value of \code{"p"} extracts the positive
#'   fraction. (See the Details section for other options)
#' @param as.percent If \code{TRUE}, return the positive fraction (or
#'   negative fraction) as a formatted character string
#' @param excluded If \code{TRUE}, return the negative fraction or
#'   negative population
#' @param fmt Character format string for \code{\link{sprintf}} if the
#'    argument \code{as.percent = TRUE}
#'
#' @return
#' For \code{getFres}: a numeric vector of the desired statistic \emph{or}
#' a character string formatted as the percent of positive or negative events.
#'
#' @import flowCore
#' 
#' @export
#' 
getFres <- function(x, filter,
		what = c("p", "true", "count", "n", "positive", "negative"),
		as.percent = FALSE, excluded = FALSE, fmt = "%4.1f%%")
{
# match what initially
	what <- match.arg(what)

# dispatch based on arguments 'x' and 'filter'
	if (class(x) %in% c("flowFrame", "flowSet")) {
		if (missing(filter))
			stop("'filter' must be provided if 'x' is a flowFrame or flowSet")
		if (class(filter) %in% c("filterResultList", "logicalFilterResult"))
			fres <- filter
		else if (is(filter, "filter"))
			fres <- filter(x, filter)
		else
			stop("'filter' must be a filter object")
	}
	if (class(x) %in% c("filterResultList", "logicalFilterResult"))
		fres <- x
# special case to accept 'what' in 2nd argument position
	if (!missing(filter) && is.character(filter) && length(filter) == 1) {
		choices <- c("p", "true", "count", "n", "positive", "negative")
		ans <- try(pmatch(filter, choices), silent = TRUE)
		if (class(ans) == "try-error")
			stop ("unable to handle second argument: ",	deparse(substitute(filter)))
		else
			what <- choices[ans]
	}
# allow for several forms of 'what'
	if (what == "n") what <- "count"
	if (what == "positive") what <- "true"
	if (what == "negative") {
		what <- "true"
		excluded <- TRUE
	}
# extract summary from filter result
	junk <- capture.output(ss <- summary(fres), file = NULL) # suppress messages
# process as list or as single object
	if (is(fres, "list")) {
		n <- sapply(ss, slot, "count")
		true <- sapply(ss, slot, "true")
	}
	else {
		n <- slot(ss, "count")
		true <- slot(ss, "true")
	}
# use negative population if 'excluded' is TRUE
	if (excluded == TRUE)
		true <- n - true
# extract return value
	if (what == "p" && as.percent) ret <- sprintf(fmt, 100*true/n)
	else if (what == "p" && !as.percent) ret <- true/n
	else if (what == "true") ret <- true
	else if (what == "count") ret <- n

	return(ret)
}
#' @name get.p
#' @rdname getFres
#' 
#' @return
#' For \code{get.p}: a numeric vector of the positive fraction, \emph{or} 
#' a character string formatted as the percent positive events.
#'
#' @export
get.p <- function(x, filter, as.percent = FALSE, excluded = FALSE)
	{
	if (missing(filter))
		getFres(x = x, what = "p", as.percent = as.percent,	
			excluded = excluded)
	else
		getFres(x = x, what = "p", filter = filter, as.percent = as.percent,
			excluded = excluded)
}
