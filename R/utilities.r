#' @export
list_cols <- function(df, return_names = TRUE) {
    classes = sapply(df, class)
    if (return_names) {
        return( colnames(df)[which(classes %in% "list")] )
    } else {
        return(which(classes %in% "list"))
    }
}
