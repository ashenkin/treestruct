#' @title BranchStruct
#' @description Constructor for BranchStruct S3 object
#' @param treestructdf branch structure; see DESCRIPTION for column names, Default: data.frame()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname BranchStruct

BranchStruct <- function(treestructdf = data.frame()) {
    stopifnot(is.data.frame(treestructdf))

    branch <-
        list(
            id <- "",
            treestruct = treestructdf,
            surface_area_total = NA
        )

    structure(branch, class = "BranchStruct")
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param newVal PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname setID

setID <- function(obj, newVal) {
    UseMethod("setID", obj)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param newVal PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname setID.default

setID.default <- function(obj, newVal) {
    warning("Object BranchStruct required")
    return(obj)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param newVal PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname setID.BranchStruct

setID.BranchStruct <- function(obj, newVal) {
    obj$id <- newVal
    return(obj)
}

