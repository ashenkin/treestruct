#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param return_names PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname list_cols

list_cols <- function(df, return_names = TRUE) {
    classes = sapply(df, class)
    if (return_names) {
        return( colnames(df)[which(classes %in% "list")] )
    } else {
        return(which(classes %in% "list"))
    }
}


#' @title set_valid_col
#' @description set overall validity property of treestruct/branchstruct members
#' @param obj TreeStruct or BranchStruct object
#' @return TreeStruct or BranchStruct object
#' @details Removes 'valid' column if it already exists, so we don't include it in determining whether a treestruct is valid,
#' and then looks for all columns with "valid" in their name.  Performs a logical 'and' on those columns in a new column, 'valid'.
#' @export

set_valid_col <- function(obj) {
    #accepts treestruct or branchstruct object
    # remove 'valid' column if it already exists, so we don't include it in determining whether a treestruct is valid
    obj$treestructs = obj$treestructs %>%
        rowwise() %>%
        mutate(valid = all(c_across(contains("valid")), na.rm = T))
    return(obj)
}
