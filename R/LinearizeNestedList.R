#' Flatten Nested Lists to Get Index  Paths to Items
#'
#' Recursively linearize a nested list. If the list is flat, do nothing.
#'
#' @author Akhil S. Bhel
#' @param NList A list of objects.
#' @param LinearizeDataFrames A boolean denoting whether NList is a data frame to be converted to a list.
#' @param NameSep A character that serves as separator between levels of nested list
#' @param ForceNames A boolean denoting whether LinearizeNestedList should return a named path to items in NList.
#'
#' @return LinearList: Flat named list
#'
#' @references {
#' Mahto A Turner R, Bhel AS. (2016) mrdwabmisc: Miscellaneous R functions, mostly for data processing.
#' }
#'
#' @source \url{https://github.com/mrdwab/mrdwabmisc/tree/master/R}
LinearizeNestedList <- function(NList,
                                LinearizeDataFrames = FALSE,
                                NameSep = "/",
                                ForceNames = FALSE) {
    # Sanity checks:
    #
    stopifnot(is.character(NameSep), length(NameSep) == 1)
    stopifnot(is.logical(LinearizeDataFrames),
              length(LinearizeDataFrames) == 1)
    stopifnot(is.logical(ForceNames), length(ForceNames) == 1)
    if (!is.list(NList))
        return(NList)
    #
    # If no names on the top-level list coerce names. Recursion shall handle
    # naming at all levels.
    #
    if (is.null(names(NList)) | ForceNames == TRUE)
        names(NList) <- as.character(1:length(NList))
    #
    # If simply a dataframe deal promptly.
    #
    if (is.data.frame(NList) & LinearizeDataFrames == FALSE)
        return(NList)
    if (is.data.frame(NList) & LinearizeDataFrames == TRUE)
        return(as.list(NList))
    #
    # Book-keeping code to employ a while loop.
    #
    A <- 1
    B <- length(NList)
    #
    # We use a while loop to deal with the fact that the length of the nested
    # list grows dynamically in the process of linearization.
    #
    while (A <= B) {
        Element <- NList[[A]]
        EName <- names(NList)[A]
        if (is.list(Element)) {
            #
            # Before and After to keep track of the status of the top-level trunk
            # below and above the current element.
            #
            if (A == 1) {
                Before <- NULL
            } else {
                Before <- NList[1:(A - 1)]
            }
            if (A == B) {
                After <- NULL
            } else {
                After <- NList[(A + 1):B]
            }
            #
            # Treat dataframes specially.
            #
            if (is.data.frame(Element)) {
                if (LinearizeDataFrames == TRUE) {
                    #
                    # `Jump` takes care of how much the list shall grow in this step.
                    #
                    Jump <- length(Element)
                    NList[[A]] <- NULL
                    #
                    # Generate or coerce names as need be.
                    #
                    if (is.null(names(Element)) |
                        ForceNames == TRUE)
                        names(Element) <-
                        as.character(1:length(Element))
                    #
                    # Just throw back as list since dataframes have no nesting.
                    #
                    Element <- as.list(Element)
                    #
                    # Update names
                    #
                    names(Element) <-
                        paste(EName, names(Element), sep = NameSep)
                    #
                    # Plug the branch back into the top-level trunk.
                    #
                    NList <- c(Before, Element, After)
                }
                Jump <- 1
            } else {
                NList[[A]] <- NULL
                #
                # Go recursive! :)
                #
                if (is.null(names(Element)) | ForceNames == TRUE)
                    names(Element) <-
                        as.character(1:length(Element))
                Element <-
                    LinearizeNestedList(Element,
                                        LinearizeDataFrames,
                                        NameSep,
                                        ForceNames)
                names(Element) <-
                    paste(EName, names(Element), sep = NameSep)
                Jump <- length(Element)
                NList <- c(Before, Element, After)
            }
        } else {
            Jump <- 1
        }
        #
        # Update book-keeping variables.
        #
        A <- A + Jump
        B <- length(NList)
    }
    return(NList)
}
