# Accessors ####

setDataset.TreeStructs <- function(obj, newVal) {
    obj$Dataset <- newVal
    return(obj)
}

#' @title setTreestruct.TreeStructs
#' @description This validates and creates the nested dataframe and populates the object with it.
#' @param obj TreeStructs object
#' @param treestructs treestruct data frame
#' @return TreeStructs object with validated, nested, treestructs dataframe.
#' @details We use the column names defined in the object properties.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname setTreestruct.TreeStructs
#'
#' @import dplyr
#' @import tidyr
setTreestruct.TreeStructs <- function(obj, treestructs) {
    newobj = obj

    # create nested dataframe that is the central piece of the object
    newobj$treestructs = treestructs %>%
        dplyr::group_by_(newobj$idcol) %>%
        tidyr::nest(.key = "treestruct")

    # validate treestructs
    if (!getOption("skip_validation", default = FALSE)) {
        valid_treestruct = validate_treestruct(newobj)
        if (valid_treestruct) return(newobj)
        else {
            warning(paste("treestructs not valid,
                           returning unmodified TreeStructs object ID ", obj$id))
            return(obj)
        }
    } else {
        warning("validation turned off, returning unvalidated TreeStructs")
        return(newobj)
    }

}

#' @export

getTreestruct.TreeStructs <- function(obj) {
    return(obj$treestructs)
}

# Validators ####

validate_treestruct.TreeStructs <- function(obj) {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE
    valid = T
    treestructs = getTreestruct(obj)
    for (this_row in 1:nrow(treestructs)) { # loop over all treestruct dfs in object
        thisTree = treestructs[[this_row, obj$idcol]]
        this_treestruct = treestructs[[this_row, "treestruct"]]
        if (verbose) warning(paste("Validating Tree",thisTree))
        valid = valid & validate_parents(1:nrow(this_treestruct), this_treestruct[[obj$parent_row_col]],
                                         parents_are_rows = T)
        valid = valid & is.data.frame(this_treestruct)
            if (!is.data.frame(this_treestruct)) warning("Treestruct not a dataframe error")

        valid = valid & validate_internodes(this_treestruct %>% mutate(internode_id = 1:nrow(this_treestruct)), ignore_error_col = NA)
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
#' @rdname calc_surfarea.TreeStructs

calc_surfarea.TreeStructs <- function(obj) {
    sa_fun <- function(x) {
        # surface area of a truncated cone
        x$surf_area = with(x, pi * (d_child/20 + d_parent/20) * sqrt((d_parent/20 - d_child/20)^2 + len^2))
        return(x)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, sa_fun)
    obj$treestructs$surface_area_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$surf_area, na.rm = T))
    return(obj)
}

#' @export

calc_vol.TreeStructs <- function(obj) {
    vol_fun <- function(x) {
        x$vol = with(x, pi * ((d_child/10 + d_parent/10)/4)^2 * len)
        return(x)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, vol_fun)
    obj$treestructs$vol_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$vol, na.rm = T))
    return(obj)
}

#' @title calc_pathlen.TreeStructs
#' @description calc path lengths
#' @param obj TreeStructs object
#' @return TreeStructs object
#' @details Details
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname calc_pathlen.TreeStructs

calc_pathlen.TreeStructs <- function(obj) {
    # overwrite treestruct::calc_pathlen
    pathlen_vec <- function(treestruct) {
        parent_row = match(treestruct$parent_id, treestruct$internode_id)
        treestruct$pathlen = calc_pathlen_cpp(treestruct$len,
                                       parent_row)
        return(treestruct)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, pathlen_vec)
    obj$treestructs$pathlen_max = map_dbl(obj$treestructs$treestruct, ~max(.$pathlen, na.rm = T))
    obj$treestructs$pathlen_mean = map_dbl(obj$treestructs$treestruct, ~mean(.$pathlen, na.rm = T))
    obj$treestructs$pathlen_frac = map_dbl(obj$treestructs, ~.$pathlen_mean/.$pathlen_max)
    return(obj)
}

#' @export

calc_len.TreeStructs <- function(obj) {
    obj$treestructs$len_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$len, na.rm = T))
    return(obj)
}
