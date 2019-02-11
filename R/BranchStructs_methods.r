# Accessors ####

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

setDataset <- function(obj, newVal) {
    UseMethod("setDataset", obj)
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
#' @rdname setDataset.default

setDataset.default <- function(obj, newVal) {
    warning("Object BranchStructs required")
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
#' @rdname setDataset.BranchStructs

setDataset.BranchStructs <- function(obj, newVal) {
    obj$Dataset <- newVal
    return(obj)
}

#' @title setTreestruct.BranchStructs
#' @description This validates and creates the nested dataframe and populates the object with it.
#' @param obj BranchStructs object
#' @param treestructs treestruct data frame
#' @return BranchStructs object with validated, nested, treestructs dataframe.
#' @details We use the column names defined in the object properties.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname setTreestruct.BranchStructs
#'
#' @import dplyr
#' @import tidyr
setTreestruct.BranchStructs <- function(obj, treestructs) {
    newobj = obj

    # create ignore_error col if it doesn't exist
    if (! any(obj$ignore_error_col %in% names(treestructs)))
        treestructs[[obj$ignore_error_col]] = rep(F, length(nrow(treestructs)))
browser()
    # set NA ignore errors to F
    treestructs = treestructs %>% tidyr::replace_na(setNames(list(F),obj$ignore_error_col)) %>%
        dplyr::mutate(!! sym(obj$ignore_error_col) := as.logical(sym(obj$ignore_error_col)))  # make sure we end up with a logical column

    # create nested dataframe that is the central piece of the object
    newobj$treestructs = treestructs %>%
        dplyr::group_by_(newobj$idcol) %>%
        tidyr::nest(.key = "treestruct")

    # validate treestructs
    valid_treestruct = validate_treestruct(newobj)
    if (valid_treestruct) return(newobj)
    else {
        warning(paste("treestructs not valid,
                       returning unmodified BranchStructs object ID ", obj$id))
        return(obj)
    }

}

#' @export

getTreestruct <- function(obj, treestruct) {
    UseMethod("getTreestruct", obj)
}

#' @export

getTreestruct.BranchStructs <- function(obj) {
    return(obj$treestructs)
}

#' @export

getTreestruct.default <- function(obj) {
    warning("Doesn't apply to this class")
    return(obj)
}

# Validators ####

validate_treestruct.BranchStructs <- function(obj) {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE
    valid = T
    treestructs = getTreestruct(obj)
    for (this_row in 1:nrow(treestructs)) { # loop over all treestruct dfs in object
        thisbranch = treestructs[[this_row, obj$idcol]]
        this_treestruct = treestructs[[this_row, "treestruct"]]
        if (verbose) warning(paste("Validating branch",thisbranch))
        valid = valid & validate_parents(this_treestruct[[obj$internodeid_col]], this_treestruct[[obj$parentid_col]],
                                         ignore_errors = this_treestruct[[obj$ignore_error_col]])
        valid = valid & is.data.frame(this_treestruct)
            if (!is.data.frame(this_treestruct)) warning("Treestruct not a dataframe error")
        browser()
        valid = valid & validate_internodes(this_treestruct, obj$internodeid_col, obj$parentid_col,
                                            ignore_error_col = obj$ignore_error_col)
    }
    # assume columns are all there.  TODO validate column names.
    return(valid)
}

# Structure Analysis ####

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname calc_surfarea.BranchStructs

calc_surfarea.BranchStructs <- function(obj) {
    attach(obj$treestruct)
    obj$treestruct$surf_area = pi * (d_child/10 + d_parent/10) * sqrt((d_parent/10 - d_child/10)^2 + len^2)
    obj$surface_area_total = sum(obj$treestruct$surf_area, na.rm = T)
    detach(obj$treestruct)
    return(obj)
}

#' @title calc_pathlen.BranchStructs
#' @description calc path lengths
#' @param obj BranchStructs object
#' @return BranchStructs object
#' @details Details
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname calc_pathlen.BranchStructs

calc_pathlen.BranchStructs <- function(obj) {
    # overwrite treestruct::calc_pathlen
    attach(obj$treestruct)
    browser()
    print(obj$id)
    parent_row = match(parent_id, internode_id)
    pathlen_vec = calc_pathlen_cpp(len,
                                   parent_row)
    obj$treestruct$path_len = pathlen_vec
    detach(obj$treestruct)
    return(obj)
}

