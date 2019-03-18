# Accessors ####

setDataset.TreeStructs <- function(obj, newVal) {
    obj$Dataset <- newVal
    return(obj)
}

#' @title setTreestruct.TreeStructs
#' @description This validates and creates the nested dataframe and populates the object with it.
#' @param obj TreeStructs object
#' @param treestructs treestruct data frame
#' @param convert_to_meters numeric.  TreeQSM uses meters as a standard unit.  If units are not in meters, pass the conversion factor here such that meters = QSMunits * convert_to_meters.  Only radius and length values are converted, not cylinder positions or axes.  A value of NA results in no conversion.  Default: NA
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
setTreestruct.TreeStructs <- function(obj, treestructs, convert_to_meters = NA) {
    newobj = obj

    # convert units to meters
    if (!is.na(convert_to_meters))
        newobj$treestructs = newobj$treestructs %>%
            dplyr::mutate(
                !!rlang::sym(obj$radius_col) := !!rlang::sym(obj$radius_col) * convert_to_meters,
                !!rlang::sym(obj$length_col) := !!rlang::sym(obj$length_col) * convert_to_meters
            )

    # create nested dataframe that is the central piece of the object
    newobj$treestructs = treestructs %>%
        dplyr::group_by_(newobj$idcol) %>%
        # add id columns to match hand branch measurements
        dplyr::mutate(!!rlang::sym(obj$internodeid_col) := 1:n(),
                      #!!rlang::sym(obj$parentid_col) :=
                          #!!rlang::sym(obj$internodeid_col)[if_else(rlang::UQ(rlang::sym(obj$parent_row_col)) %in% 0, NA, !!rlang::sym(obj$parent_row_col))]) %>%
                      # TODO fix the direct reference to column names below.  can't figure out how to make it dynamic in ifelse...
                      !!rlang::sym(obj$parentid_col) := internode_id[ifelse(parent_row %in% 0, NA, parent_row)], # vector index ignores 0
                      d_parent = !!rlang::sym(obj$radius_col)*2, # to make compliant with branchstructs.  TODO remove once treestructs no longer inherits from branchstructs
                      d_child = d_parent) %>%
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
setTips.TreeStructs <- function(obj) {
    istip <- function(ts) {
        # a branch is a tip if no other branch claims it as a parent
        ts$tip = ! ts[[obj$internodeid_col]] %in% ts[[obj$parentid_col]]
        return(ts)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, istip)
    obj$tips_set = T
    return(obj)
}

# Validators ####

validate_treestruct.TreeStructs <- function(obj) {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE
    valid = T
    treestructs = obj$treestructs
    for (this_row in 1:nrow(treestructs)) { # loop over all treestruct dfs in object
        thisTree = treestructs[[this_row, obj$idcol]]
        this_treestruct = treestructs[[this_row, "treestruct"]]
        if (verbose) message(paste("Validating Tree",thisTree))
        valid = valid & validate_parents(1:nrow(this_treestruct), this_treestruct[[obj$parent_row_col]],
                                         parents_are_rows = T)
        valid = valid & is.data.frame(this_treestruct)
            if (!is.data.frame(this_treestruct)) warning("Treestruct not a dataframe error")

        valid = valid & validate_internodes(this_treestruct %>% mutate(internode_id = 1:nrow(this_treestruct)), ignore_error_col = NA)
    }
    # assume columns are all there.  TODO validate column names.
    return(valid)
}

#' @export
setTips.TreeStructs <- function(obj) {
    istip <- function(ts) {
        # a branch is a tip if it has no daughters
        ts$tip = ts[[obj$daughter_row_col]] == 0
        return(ts)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, istip)
    obj$tips_set = T
    return(obj)
}

# Housekeeping ####

#' @export
make_compatible.TreeStructs <- function(obj) {
    getTreestruct(obj) %>% select(file:branch, len, parent_row, daughter_row, internode_id:pathlen)
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
        # x is a treestruct dataframe
        # surface area of a cylinder
        # QSM units are meters.
        x$surf_area = with(x, 2 * pi * get(obj$radius_col) * get(obj$length_col))
        return(x)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, sa_fun)
    obj$treestructs$surface_area_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$surf_area, na.rm = T))
    return(obj)
}

#' @export

calc_vol.TreeStructs <- function(obj) {
    vol_fun <- function(x) {
        x$vol = with(x, pi * get(obj$radius_col)^2 * get(obj$length_col))
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

    valid_internode_order = map_lgl(obj$treestructs$treestruct, ~validate_internode_order(.[[obj$parent_row_col]], parents_are_rows = T))
    if (any(!valid_internode_order)) stop(paste("Invalid internode order:", paste(obj$treestructs[[obj$idcol]], collapse = ", ")))

    pathlen_vec <- function(treestruct) {
        treestruct$pathlen = calc_pathlen_cpp(treestruct[[obj$length_col]],
                                              treestruct[[obj$parent_row_col]])
        return(treestruct)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, pathlen_vec)
    obj$treestructs$pathlen_max = map_dbl(obj$treestructs$treestruct, ~max(.$pathlen, na.rm = T))
    # mean pathlength to tips of QSM
    obj$treestructs$pathlen_mean = map_dbl(obj$treestructs$treestruct,
                                           .f = function(x) x %>%
                                               dplyr::filter(!!sym(obj$daughter_row_col) == 0) %>%
                                               dplyr::summarize(pathlen_mean = mean(pathlen, na.rm = T)) %>%
                                               dplyr::pull(pathlen_mean))

    obj$treestructs$pathlen_frac =  with(obj$treestructs, pathlen_mean/pathlen_max)
    obj$treestructs$len_tot =  map_dbl(obj$treestructs$treestruct, ~sum(.[[obj$length_col]], na.rm = T))
    return(obj)
}

