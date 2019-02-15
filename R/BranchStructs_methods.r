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
#' @param convert_to_meters Hand-measured branches typically use mm for diameter and cm for length.  treestruct uses meters as a standard unit.
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
setTreestruct.BranchStructs <- function(obj, treestructs, convert_to_meters = T) {
    newobj = obj

    # create ignore_error col if it doesn't exist
    if (! any(obj$ignore_error_col %in% names(treestructs)))
        treestructs[[obj$ignore_error_col]] = rep(F, length(nrow(treestructs)))

    # set NA ignore errors to F
    # TODO make ignore_error dynamic based on obj$ignore_error_col value
    treestructs = treestructs %>% tidyr::replace_na(setNames(list(F),obj$ignore_error_col)) %>%
                  dplyr::mutate(ignore_error := as.logical(ignore_error))  # make sure we end up with a logical column

    # convert units to meters
    if (convert_to_meters) treestructs = treestructs %>%
                                         dplyr::mutate(
                                             !!rlang::sym(obj$length_col) := !!rlang::sym(obj$length_col) / 10,
                                             !!rlang::sym(obj$d_child_col) := !!rlang::sym(obj$d_child_col) / 100,
                                             !!rlang::sym(obj$d_parent_col) := !!rlang::sym(obj$d_parent_col) / 100
                                         )

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
                           returning unmodified BranchStructs object ID ", obj$id))
            return(obj)
        }
    } else {
        warning("validation turned off, returning unvalidated BranchStructs")
        return(newobj)
    }

}

#' @export

getTreestruct <- function(obj, treestruct) {
    UseMethod("getTreestruct", obj)
}

#' @export

getTreestruct.BranchStructs <- function(obj) {
    return(obj$treestructs$treestruct)
}

#' @export

getTreestruct.default <- function(obj) {
    warning("Doesn't apply to this class")
    return(obj)
}

#' @export
setTips <- function(obj) {
    UseMethod("setTips", obj)
}

#' @export
setTips.BranchStructs <- function(obj) {
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

validate_treestruct.BranchStructs <- function(obj) {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE
    valid = T
    treestructs = obj$treestructs
    for (this_row in 1:nrow(treestructs)) { # loop over all treestruct dfs in object
        thisbranch = treestructs[[this_row, obj$idcol]]
        this_treestruct = treestructs[[this_row, "treestruct"]]
        if (verbose) message(paste("Validating branch",thisbranch))
        valid = valid & validate_parents(this_treestruct[[obj$internodeid_col]], this_treestruct[[obj$parentid_col]],
                                         ignore_errors = this_treestruct[[obj$ignore_error_col]])
        valid = valid & is.data.frame(this_treestruct)
            if (!is.data.frame(this_treestruct)) warning("Treestruct not a dataframe error")

        valid = valid & validate_internodes(this_treestruct, obj$internodeid_col,
                                            ignore_error_col = obj$ignore_error_col)
    }
    # assume columns are all there.  TODO validate column names.
    return(valid)
}

#' @export

reorder_internodes <- function(obj) {
    UseMethod("reorder_internodes", obj)
}

#' @export

reorder_internodes.BranchStructs <- function(obj) {
    if (! obj$tips_set) obj = setTips(obj)

    reorder_treestruct <- function(ts) {
        tscopy = ts
        ntips = sum(tscopy$tip)
        tscopy$orig_row = 1:nrow(ts)
        tscopy$parent_row = match(tscopy[[obj$parentid_col]], tscopy[[obj$internodeid_col]])
        # which branch does this cylinder belong to
            tscopy$branchnum = NA
        # the sequence order of the cylinders within the branch
            tscopy$branchorder = NA
        # initiate branches from tip downwards
            tscopy[tscopy$tip,]$branchnum = 1:ntips
            tscopy[tscopy$tip,]$branchorder = 1
        # tscopy$daughter_row = match(tscopy[[obj$internodeid_col]], tscopy[[obj$parentid_col]])
        # follow each tip downwards, stopping when you hit the base or a cylinder already claimed by another branch
        for (thisbranch in 1:ntips) {
            curr_row = which(tscopy$branchnum == thisbranch)
            thisbranchorder = 2
            repeat {
                next_row = tscopy[curr_row,]$parent_row
                if (is.na(next_row) |                           # we've hit the base
                    ! is.na(tscopy[next_row,]$branchnum)) break # cylinder already claimed by another branch
                tscopy[next_row,]$branchnum = thisbranch
                tscopy[next_row,]$branchorder = thisbranchorder
                thisbranchorder = thisbranchorder + 1
                curr_row = next_row
            }
        }
        # now reorder ts from last branch to first, and tip downwards within branch
            tscopy = tscopy %>% dplyr::arrange(branchnum, desc(branchorder))# %>%
                #dplyr::arrange(dplyr::desc(dplyr::row_number()))
            ts = ts[tscopy$orig_row,]    # use copied dataframe to reorder original (we've added columns to tscopy, etc)
            return(tscopy)
            #return(ts)
    }

    obj$treestructs$treestruct = purrr::map(getTreestruct(obj), reorder_treestruct)
    return(obj)
}



# reorder_internodes.BranchStructs <- function(obj) {
#     reorder_treestruct <- function(ts) {
#
#         parent_row = match(ts[[obj$parentid_col]], ts[[obj$internodeid_col]])
#         if (all(parent_row - 1:length(parent_row) < 0, na.rm = T)) return(ts)
#
#         # need to reorder rows
#         i = 1
#         repeat{
#             ts = ts[order(parent_row, na.last = F),]
#             parent_row = match(ts[[obj$parentid_col]], ts[[obj$internodeid_col]])
#             correct_order = all(parent_row - 1:length(parent_row) < 0, na.rm = T)
#             i = i + 1
#             if (i > 10) { browser() }
#             if (correct_order) break
#         }
#         return(ts)
#     }
#     obj$treestructs$treestruct = map(getTreestruct(obj), reorder_treestruct)
#     return(obj)
# }

#' @export

reorder_internodes.default <- function(obj) {
    warning("Doesn't apply to this class")
    return(obj)
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
    sa_fun <- function(x) {
        # TODO make colnames dynamic
        # surface area of a truncated cone
        x$surf_area = with(x, pi * (d_child/20 + d_parent/20) * sqrt((d_parent/20 - d_child/20)^2 + len^2))
        return(x)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, sa_fun)
    obj$treestructs$surface_area_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$surf_area, na.rm = T))
    return(obj)
}

#' @export
calc_vol <- function(obj) {
    UseMethod("calc_vol", obj)
}

#' @export
calc_vol.BranchStructs <- function(obj) {
    vol_fun <- function(x) {
        x$vol = with(x, pi * ((d_child/10 + d_parent/10)/4)^2 * len)
        return(x)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, vol_fun)
    obj$treestructs$vol_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$vol, na.rm = T))
    return(obj)
}

#' @export
calc_vol.default <- function(obj) {
    warning("Doesn't apply to this class")
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
    pathlen_vec <- function(treestruct) {
        parent_row = match(treestruct$parent_id, treestruct$internode_id)
        print(treestruct$len)
        print(parent_row)
        treestruct$pathlen = calc_pathlen_cpp(treestruct$len,
                                       parent_row)
        return(treestruct)
    }

    valid_internode_order = map_lgl(obj$treestructs$treestruct, ~validate_internode_order(.[[obj$parentid_col]], .[[obj$internodeid_col]], parents_are_rows = F))
    if (any(!valid_internode_order)) stop(paste("Invalid internode order:", paste(obj$treestructs[[obj$idcol]][!valid_internode_order], collapse = ", ")))

    obj$treestructs$treestruct = map(obj$treestructs$treestruct, pathlen_vec)
    obj$treestructs$pathlen_max = map_dbl(obj$treestructs$treestruct, ~max(.$pathlen, na.rm = T))
    obj$treestructs$pathlen_mean = map_dbl(obj$treestructs$treestruct, ~mean(.$pathlen, na.rm = T))
    obj$treestructs$pathlen_frac = map_dbl(obj$treestructs, ~.$pathlen_mean/.$pathlen_max)
    obj$treestructs$len_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$len, na.rm = T))
    return(obj)
}

#' @export
calc_len <- function(obj) {
    UseMethod("calc_len", obj)
}

#' @export
calc_len.BranchStructs <- function(obj) {
    obj$treestructs$len_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$len, na.rm = T))
    return(obj)
}
