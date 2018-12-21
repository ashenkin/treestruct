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

setTreestruct.BranchStructs <- function(obj, treestructs) {
    # this creates the nested dataframe.  We use the column names defined in the object properties.
    newobj = obj
    newobj$treestructs = treestructs

    # validate treestructs
    valid_treestruct = validate_treestruct(newobj)
    if (valid_treestruct) return(newobj)
    else {
        warning(paste( "treestructs not valid,
                       returning unmodified BranchStructs object ID ", obj$id))
        return(obj)
    }

    # create nested dataframe that is the central piece of the object

}

# Validators ####

validate_treestruct.BranchStructs <- function(obj) {
    valid = validate_parents(obj$treestructs$internode_id, obj$treestructs$parent_id)
    valid = valid & is.data.frame(obj$treestructs)
    valid = valid & validate_internodes(obj$treestructs, "internode_id", "parent_id")
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
#' @rdname calc_surfarea

calc_surfarea <- function(obj) {
    UseMethod("calc_surfarea", obj)
}

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
#' @rdname calc_surfarea.default

calc_surfarea.default <- function(obj) {
    warning("Object BranchStruct required")
    return(obj)
}

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
#' @rdname calc_surfarea.BranchStruct

calc_surfarea.BranchStruct <- function(obj) {
    attach(obj$treestruct)
    obj$treestruct$surf_area = pi * (d_child/10 + d_parent/10) * sqrt((d_parent/10 - d_child/10)^2 + len^2)
    obj$surface_area_total = sum(obj$treestruct$surf_area, na.rm = T)
    detach(obj$treestruct)
    return(obj)
}

#' @title calc_pathlen
#' @export
#' @rdname calc_pathlen
calc_pathlen <- function(obj) {
    UseMethod("calc_pathlen", obj)
}

#' @title calc_pathlen.BranchStruct
#' @description calc path lengths
#' @param obj BranchStruct object
#' @return BranchStruct object
#' @details overwrites treestruct::calc_pathlen
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname calc_pathlen.BranchStruct

calc_pathlen.BranchStruct <- function(obj) {
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

#' @title calc_pathlen.Default
#' @description Wrapper for C++ pathlength code
#' @param tree_structure a tree structure dataframe
#' @param length_col chr specifiying length column. default:"len"
#' @param parent_row_col chr specifiying parent row column. default:"parent_row"
#' @param path_len_col chr specifiying path length column. default:"path_len"
#' @return tree structure dataframe with a path_len column populated
#' @details DETAILS
#' @export
#' @rdname calc_pathlen.Default
calc_pathlen.Default <- function(tree_structure, length_col = "len",
                         parent_row_col = "parent_row",
                         path_len_col = "path_len") {
  # wrap cpp function
  pathlen_vec = calc_pathlen_cpp(tree_structure[[length_col]],
                                 tree_structure[[parent_row_col]])
  tree_structure[[path_len_col]] = pathlen_vec
  return(tree_structure)
}
