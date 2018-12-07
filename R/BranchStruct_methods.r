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
    if (grepl("\\s", newVal)) warning("whitespace removed from ID")
    obj$id <- gsub("\\s", "", newVal) # remove whitespace
    return(obj)
}

setTreestruct <- function(obj, treestruct) {
    UseMethod("setTreestruct", obj)
}

setTreestruct.BranchStruct <- function(obj, treestruct) {
    newobj = obj
    newobj$treestruct = treestruct
    valid_treestruct = validate_treestruct(newobj)
    if (valid_treestruct) return(newobj)
    else {
        warning(paste( "treestruct not valid,
                returning unmodified BranchStruct object ID ", obj$id))
        return(obj)
    }

}

# Validators ####

validate_treestruct <- function(obj) {
    UseMethod("validate_treestruct", obj)
}

#' @importFrom settings clone_and_merge
validate_treestruct.BranchStruct <- function(obj) {
    valid = validate_parents(obj$treestruct$internode_id, obj$treestruct$parent_id)
    valid = valid & is.data.frame(obj$treestruct)
    valid = valid & validate_internodes(obj$treestruct, "internode_id", "parent_id")
    return(valid)
}

validate_treestruct.Default <- function(...) {
    valid = treestruct::validate_parents(...)
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
