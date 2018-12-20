#' @title BranchStructs
#' @description Constructor for BranchStructs S3 object
#' @param dataset name of dataset
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname BranchStructs

BranchStructs <- function(dataset = NA) {
    stopifnot(length(dataset) == 1)

    branches <-
        list(
            dataset = dataset,
            surface_area_total = NA,
            pathlen_total = NA,
            treestructs = NA
        )

    branches = structure(branches, class = "BranchStructs")
    return(branches)
}



