# TODO import individual functions and not entire library when not necessary
# in particular: beware of using data.table and dplyr together.

#' @useDynLib treestruct, .registration = TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarize mutate pull
#' @importFrom foreach "%dopar%"
#' @importFrom stats sd

.onUnload <- function (libpath) {
    library.dynam.unload("treestruct", libpath)
}
