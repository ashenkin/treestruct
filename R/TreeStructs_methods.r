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
        dplyr::mutate(!!rlang::sym(obj$internodeid_col) := 1:dplyr::n(),
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
getCylSummary <- function(obj, ...) {
    UseMethod("getCylSummary", obj)
}

#' @export
getCylSummary.TreeStructs <- function(obj, idx = NA, concat = T) {
    if (! is.na(idx))
        return(getTreestructs(obj)$cyl_summ[[idx]])
    else if (concat) {
        return(obj$treestructs %>% tidyr::unnest(cyl_summ))
    } else {
        return(getTreestructs(obj)$cyl_summ)
    }
}

# Validators ####

#' @import crayon
validate_treestruct.TreeStructs <- function(obj) {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE
    valid = T
    treestructs = obj$treestructs
    for (this_row in 1:nrow(treestructs)) { # loop over all treestruct dfs in object
        thisTree = treestructs[[this_row, obj$idcol]]
        this_treestruct = treestructs[[this_row, "treestruct"]]
        if (verbose) message(paste("Validating Tree",thisTree))
        this_valid = T
        this_valid = this_valid & validate_parents(1:nrow(this_treestruct), this_treestruct[[obj$parent_row_col]],
                                         parents_are_rows = T)
        this_valid = this_valid & is.data.frame(this_treestruct)
            if (!is.data.frame(this_treestruct)) warning("Treestruct not a dataframe error")

        this_valid = this_valid & validate_internodes(this_treestruct %>% dplyr::mutate(internode_id = 1:nrow(this_treestruct)), ignore_error_col = NA)
        if (verbose) message(paste(thisTree, ifelse(this_valid, "passed", crayon::red("failed")), "validation"))
        valid = valid & this_valid

    }
    # assume columns are all there.  TODO validate column names.
    return(valid)
}

# Utilities ####

#' @export
make_compatible.TreeStructs <- function(obj) {
    getTreestruct(obj) %>% dplyr::select(file:branch, len, parent_row, daughter_row, internode_id:pathlen)
}

#' @export
calc_summary_cyls <- function(obj, ...) {
    UseMethod("calc_summary_cyls", obj)
}

#' @export
calc_summary_cyls.TreeStructs <- function(obj) {
    if (! check_property(obj, "branchnums_assigned")) obj = assign_branch_num(obj)

    obj$treestructs$cyl_summ = purrr::map(getTreestruct(obj, concat = FALSE), calc_summary_cyls, check_property(obj, "furcations_corrected"))
    return(obj)
}

#' @export
calc_summary_cyls.default <- function(ts, furcations_corrected = F) {
    # assign_branch_num must be run first

    if (! furcations_corrected) {
        ts = correct_furcations(ts)
    }

    # collapse by branchnum
    cyl_summ = ts %>% dplyr::group_by(branchnum) %>%
        dplyr::summarize(
            parent_branchnum = dplyr::first(parent_branchnum),
            rad = mean(rad, na.rm = T),
            len = sum(len, na.rm = T),
            pathlen = mean(pathlen, na.rm = T),
            n_furcation = max(n_furcation, na.rm = T),
            num_cyls_in_branch = dplyr::n(),
            surf_area = {if ( "surf_area" %in% names(ts) ) sum(surf_area) else 2*pi*rad*len},
            vol = {if ( "vol" %in% names(ts) ) sum(vol) else pi*rad^2*len}
        )

    # join parent branch metrics to each row to facilitate branch scaling calculations
    cyl_summ = cyl_summ %>%
        dplyr::left_join(cyl_summ %>%
                             dplyr::select(branchnum,
                                           rad_parent = rad,
                                           len_parent = len,
                                           n_furcation_parent = n_furcation),
                         # rename(parent_rad = rad_mean,
                         #        parent_len = len,
                         #        parent_n_furcation = n_furcation),
                         by = c("parent_branchnum" = "branchnum")) %>%
        # rename cols for compatibility with treestruct dataframe
        dplyr::rename(parent_id = parent_branchnum, internode_id = branchnum) %>%
        # add parent_row for compatibility with some routines
        dplyr::mutate(parent_row = parent_row(parent_id, internode_id))

    cyl_summ = setTips(cyl_summ)

    return(cyl_summ)

}

#' @export
find_first_branch <- function(obj, ...) {
    UseMethod("find_first_branch", obj)
}

#' @title find_first_branch.TreeStructs
#' @description Finds the first branches of treestructs
#' @param obj Treestructs object
#' @param daughter_threshold Proportion of the radius of the parent branch that daughter
#'   branches must attain for a furction to be counted as the first branch, Default: 0.25
#' @return internode_id of parent branch where first furcation occurs
#' @details Returns NA if no significant furcations found
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname find_first_branch.default
find_first_branch.TreeStructs <- function(obj, daughter_threshold = 0.25) {
    if (! check_property(obj, "internodes_reordered")) obj = reorder_internodes(obj)
    if (! check_property(obj, "furcations_corrected")) obj = correct_furcations(obj)

    # TODO: should this work on cyl_summ instead of treestruct?  Not sure what constraints treeQSM puts on daughter radii...
    # Potential issue is mapping back from cyl_summ to treestruct..
    obj$treestructs$first_branch_id = purrr::map_int(getTreestruct(obj, concat = FALSE), find_first_branch, daughter_threshold = daughter_threshold)
    obj$first_branch_assigned = TRUE
    return(obj)
}

#' @title find_first_branch.default
#' @description Finds the first branch in a treestruct dataframe
#' @param ts treestruct dataframe
#' @param daughter_threshold Proportion of the radius of the parent branch that daughter
#'   branches must attain for a furction to be counted as the first branch, Default: 0.25
#' @return internode_id of parent branch where first furcation occurs
#' @details Returns NA if no significant furcations found
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname find_first_branch.default

find_first_branch.default <- function(ts, daughter_threshold = 0.25) {
    # assume ts is ordered properly
    # climb from bottom upwards, looking for first major furcation
    ts = ts[nrow(ts):1,]
    furc_rows = which(ts$n_furcation > 1)
    for (this_furc_row in furc_rows) {
        daughter_radii = ts %>% dplyr::filter(parent_id == ts[this_furc_row,]$internode_id) %>% dplyr::pull(rad)
        # major furcation if two daughters at least 1/4 as wide as parent
        num_large_daughters = sum(daughter_radii > ts[this_furc_row,]$rad * daughter_threshold)
        if (num_large_daughters >= 2) {
            return(ts[this_furc_row,]$internode_id)
        }
    }
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE
    if (verbose) message("No significant furcation found that could be called a first branch")
    return(NA)
}

#' @export
assign_cyls_to_crown <- function(obj, ...) {
    UseMethod("assign_cyls_to_crown", obj)
}

#' @title assign_cyls_to_crown.TreeStructs
#' @description Assigns cylinders to crown and stem based on where the first branch occurs
#' @param obj Treestructs object
#' @param daughter_threshold Proportion of the radius of the parent branch that daughter
#'   branches must attain for a furction to be counted as the first branch, Default: 0.25
#' @return treestruct object with crown column assigned
#' @details Returns NA if no significant furcations found
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname assign_cyls_to_crown.default
assign_cyls_to_crown.TreeStructs <- function(obj, daughter_threshold = 0.25) {
    if (! check_property(obj, "internodes_reordered")) obj = reorder_internodes(obj)
    if (! check_property(obj, "furcations_corrected")) obj = correct_furcations(obj)
    if (! check_property(obj, "first_branch_assigned")) obj = find_first_branch(obj, daughter_threshold)

    # TODO figure out how to pass arguments by name (tried with pmap, but wasn't working)
    obj$treestructs$treestruct = purrr::map2(getTreestruct(obj, concat = FALSE),
                                             getTreestructs(obj)$first_branch_id,
                                             assign_cyls_to_crown.default)
    obj$cyls_assigned_to_crown = T
    return(obj)
}

#' @title assign_cyls_to_crown.default
#' @description Assigns logical crown column in a treestruct dataframe
#' @param ts treestruct dataframe
#' @return treestruct dataframe
#' @details Returns NA in crown column if no first branch was found
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname assign_cyls_to_crown.default

assign_cyls_to_crown.default <- function(ts, first_branch_id) {
    # assume ts is ordered properly, from tip downwards
    # assign all cylinders to crown until reaching first branch
    if (is.na(first_branch_id)) {
        ts$crown = NA
        warning("No first branch found, returning treestruct with crown = NA")
        return(ts)
    }
    ts$crown = T
    first_branch_row = which(ts$internode_id == first_branch_id)
    ts[first_branch_row:nrow(ts),]$crown = F
    return(ts)
}

# Structure Analysis ####

#' @export
calc_dbh <- function(obj) {
    UseMethod("calc_dbh", obj)
}

#' @export
calc_dbh.TreeStructs <- function(obj) {
    dbhlist = map(obj$treestructs$treestruct, calc_dbh_ts)
    obj$treestructs$dbh = unlist(lapply(dbhlist, `[[`, "dbh"))
    obj$treestructs$pom = unlist(lapply(dbhlist, `[[`, "pom"))
    return(obj)
}

#' @export
calc_max_height <- function(obj) {
    UseMethod("calc_max_height", obj)
}

#' @export
calc_max_height.TreeStructs <- function(obj) {
    obj$treestructs$treeheight = unlist(map(obj$treestructs$treestruct, calc_max_height_ts))
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


#' @export
calc_sa_above.TreeStructs <- function(obj) {
    obj$treestructs$cyl_summ = purrr::map(getCylSummary(obj, concat = FALSE), calc_sa_above)
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

    valid_internode_order = map_lgl(obj$treestructs$treestruct, ~validate_internode_order(.[[obj$parentid_col]], .[[obj$internodeid_col]], parents_are_rows = F))
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
                                               dplyr::ungroup() %>% # grouping is carrying through from somewhere
                                               dplyr::filter(tip) %>%
                                               dplyr::summarize(pathlen_mean = mean(pathlen, na.rm = T)) %>%
                                               dplyr::pull(pathlen_mean))

    obj$treestructs$pathlen_frac =  with(obj$treestructs, pathlen_mean/pathlen_max)
    obj$treestructs$len_tot =  map_dbl(obj$treestructs$treestruct, ~sum(.[[obj$length_col]], na.rm = T))
    return(obj)
}

# crown metrics

#' @export
make_convhull <- function(obj) {
    UseMethod("make_convhull", obj)
}

#' @export
make_convhull.TreeStructs <- function(obj) {

    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE

    if ((! check_property(obj, "cyls_assigned_to_crown")) & obj$trees_not_branches) obj = assign_cyls_to_crown(obj)

    convhulls = purrr::map(getTreestruct(obj, concat = FALSE), make_convhull.default, trees_not_branches = obj$trees_not_branches)

    if (verbose) message("convex hulls created")

    # make sure to overwrite columns
    suppressWarnings(
        obj$treestructs <- obj$treestructs %>% dplyr::select(-one_of("^convhull$", "^convhull2d$", "^convhull2d_vert$", "^crown_vol_convhull$",
                                                        "^crown_surfarea_convhull$", "^crown_proj_area_convhull$",
                                                        "^crown_proj_area_vert_convhull$"))
    )

    #TODO the assignment below is crazy slow.  do it with data.table maybe
    # from https://jennybc.github.io/purrr-tutorial/ls01_map-name-position-shortcuts.html
    obj$treestructs = obj$treestructs %>%
        dplyr::bind_cols(
            convhulls %>% {
            tibble::tibble(
                convhull = map(., "convhull"),
                convhull2d = map(., "convhull2d"),
                convhull2d_vert = map(., "convhull2d_vert"),
                crown_vol_convhull = map_dbl(., "crown_vol_convhull"),
                crown_surfarea_convhull = map_dbl(., "crown_surfarea_convhull"),
                crown_proj_area_convhull = map_dbl(., "crown_proj_area_convhull"),
                crown_proj_area_vert_convhull = map_dbl(., "crown_proj_area_vert_convhull")
            )}
        )

    return(obj)
}

#' @export
make_convhull.default <- function(ts, trees_not_branches) {
    #TODO include both start and end points of cylinders

    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE

    crown_defined = trees_not_branches & !any(is.na(ts$crown))

    tryCatch({
        this_convhull = NA
        this_vert2d_convhull = NA
        #  make 2d convhull even if couldn't find first branch
        this_2dconvhull = geometry::convhulln(ts[,c("x_start","y_start")], options = "FA")

        # convhulls for trees (same as branches, but don't make some if we can't ID first branch)
        if (crown_defined) {
            ts = subset(ts, ts$crown) # make crown convex hulls for full trees
            this_convhull = geometry::convhulln(ts[,c("x_start","y_start","z_start")], options = "FA")
            # TODO make vert2d convhull correct somehow.  it only takes one aspect.  good enough for depth, but not sail area.
            this_vert2d_convhull = geometry::convhulln(ts[,c("y_start","z_start")], options = "FA")

        } else if (! trees_not_branches) {
        # convhulls for branches
            this_convhull = geometry::convhulln(ts[,c("x_start","y_start","z_start")], options = "FA")
            this_vert2d_convhull = geometry::convhulln(ts[,c("y_start","z_start")], options = "FA")
        }

    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    return(list(convhull = this_convhull,
                convhull2d = this_2dconvhull,
                convhull2d_vert = this_vert2d_convhull,
                crown_vol_convhull = ifelse(crown_defined, this_convhull$vol, NA),
                crown_surfarea_convhull = ifelse(crown_defined, this_convhull$area, NA),
                crown_proj_area_convhull = ifelse(crown_defined, this_2dconvhull$vol, NA),
                crown_proj_area_vert_convhull = ifelse(crown_defined, this_vert2d_convhull$vol, NA)))
}

#' @export
radius_scaling.TreeStructs <- function(obj) {
    if (! "cyl_summ" %in% names(obj$treestructs)) obj = calc_summary_cyls(obj)
    obj$treestructs$cyl_summ = purrr::map(getCylSummary(obj, concat = FALSE), radius_scaling.default)
    obj$treestructs$a_median = purrr::map_dbl(getCylSummary(obj, concat = FALSE), function(x) median(x$a, na.rm = T))
    return(obj)
}

#' @export
length_scaling.TreeStructs <- function(obj) {
    if (! "cyl_summ" %in% names(obj$treestructs)) obj = calc_summary_cyls(obj)
    obj$treestructs$cyl_summ = purrr::map(getCylSummary(obj, concat = FALSE), length_scaling.default)
    obj$treestructs$b_median = purrr::map_dbl(getCylSummary(obj, concat = FALSE), function(x) median(x$b, na.rm = T))
    return(obj)
}

#' @export
run_all.TreeStructs <- function(obj, calc_dbh = T) {
    if (obj$trees_not_branches) {
        obj = assign_cyls_to_crown(obj) # only assign crown columns for whole trees, not hand-measured branches
    }
    obj = make_convhull(obj) # convhull won't work on hand measured branches
    obj = run_all.default(obj, calc_dbh)
    obj = calc_sa_above(obj)
    return(obj)
}
