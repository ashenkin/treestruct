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

BranchStruct <- function(treestructdf = data.frame(), id = NA) {
    stopifnot(is.data.frame(treestructdf))
    stopifnot(length(id) == 1)

    branch <-
        list(
            id = id,
            treestruct = NA,
            surface_area_total = NA
        )

    branch = structure(branch, class = "BranchStruct")
    branch = setTreestruct(branch, treestructdf)
    return(branch)
}



