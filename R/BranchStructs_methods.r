# Accessors ####

#' @title getFile
#' @description file list accessor
#' @param obj TreeStructs or BranchStructs object
#' @return vector of base file names
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname getFile

getFile <- function(obj) {
    UseMethod("getFile", obj)
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
#' @rdname getFile.default

getFile.default <- function(obj) {
    return(getTreestructs(obj)$file)
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

#' @export
setTreestruct <- function(obj, treestruct, ...) {
    UseMethod("setTreestruct", obj)
}

#' @title setTreestruct.BranchStructs
#' @description This validates and creates the nested dataframe and populates the object with it.
#' @param obj BranchStructs object
#' @param treestructs treestruct data frame
#' @param convert_to_meters Hand-measured branches typically use mm for diameter and cm for length.  treestruct uses meters as a standard unit.
#' @return BranchStructs object with validated, nested, treestructs dataframe.
#' @details We use the column names defined in the object properties.  convert_to_meters is the amount to divide by to get meters.  cm = 100, mm = 1000, etc.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname setTreestruct.BranchStructs
#'
setTreestruct.BranchStructs <- function(obj, treestructs, convert_to_meters = T, abort_when_invalid = F) {
    # TODO this should really be renamed importTreestruct.  setTreestruct should be a simple assignment method I think...
    newobj = obj

    # create ignore_error col if it doesn't exist
    if (! any(obj$ignore_error_col %in% names(treestructs)))
        treestructs[[obj$ignore_error_col]] = rep(F, length(nrow(treestructs)))

    # set NA ignore errors to F
    # TODO make ignore_error dynamic based on obj$ignore_error_col value
    treestructs = treestructs %>%
        tidyr::replace_na(setNames(list(F),obj$ignore_error_col)) %>%
        dplyr::ungroup() %>% # grouped variables from previous operations were giving errors here for some reason...
        dplyr::mutate(ignore_error := as.logical(ignore_error))  # make sure we end up with a logical column

    # convert units to meters
    if (convert_to_meters) {
        if (is.logical(convert_to_meters) & convert_to_meters) { #default conversion
            treestructs = treestructs %>%
                                         dplyr::mutate(
                                             !!rlang::sym(obj$length_col) := !!rlang::sym(obj$length_col) / 100, # hand-measured lengths are in cm.  Convert to meters.
                                             !!rlang::sym(obj$d_child_col) := !!rlang::sym(obj$d_child_col) / 1000, # hand-measured diams are in mm.  Convert to meters.
                                             !!rlang::sym(obj$d_parent_col) := !!rlang::sym(obj$d_parent_col) / 1000,
                                             # add radius column for compatibility with TreeStructs
                                             !!rlang::sym(obj$radius_col) := (!!rlang::sym(obj$d_parent_col) + !!rlang::sym(obj$d_child_col)) / 2 / 2 # mean of parent+child, div by 2 for radius instead of diameter
                                         )

        } else {
            # convert_to_meters is the amount to divide by to get meters.  cm = 100, mm = 1000, etc.
            treestructs = treestructs %>%
                dplyr::mutate(
                    !!rlang::sym(obj$length_col) := !!rlang::sym(obj$length_col) / convert_to_meters,
                    !!rlang::sym(obj$d_child_col) := !!rlang::sym(obj$d_child_col) / convert_to_meters,
                    !!rlang::sym(obj$d_parent_col) := !!rlang::sym(obj$d_parent_col) / convert_to_meters,
                    # add radius column for compatibility with TreeStructs
                    !!rlang::sym(obj$radius_col) := (!!rlang::sym(obj$d_parent_col) + !!rlang::sym(obj$d_child_col)) / 2 / 2 # mean of parent+child, div by 2 for radius instead of diameter
                )
        }
    }

    else { # convert_to_meters = F
        treestructs = treestructs %>%
            dplyr::mutate(
                # add radius column for compatibility with TreeStructs
                !!rlang::sym(obj$radius_col) := (!!rlang::sym(obj$d_parent_col) + !!rlang::sym(obj$d_child_col)) / 2 / 2 # mean of parent+child, div by 2 for radius instead of diameter
            )
    }

    # create nested dataframe that is the central piece of the object
    newobj$treestructs = treestructs %>%
        dplyr::group_by_(newobj$idcol) %>%
        tidyr::nest(.key = "treestruct")

    # reorder internodes - necessary for many operations
    if (check_property(newobj, "has_topology")) newobj = reorder_internodes(newobj)

    # validate treestructs
    if (!getOption("skip_validation", default = FALSE) & check_property(newobj, "has_topology")) {
        validation = validate_treestruct(newobj)
        valid_treestruct = validation[["all_valid"]]
        newobj = validation[["branchstructs"]]

        if (valid_treestruct | !abort_when_invalid) {

            if (!valid_treestruct)
                warning(paste("treestructs not valid, continuing with validation:", newobj$id))

            newobj = setGraph(newobj)
            valid_graph = validate_connectivity(newobj)
            if (!valid_graph) warning(crayon::red("Multiple components in one or more treestructs"))
            return(newobj)

        } else {
            warning(paste("treestructs not valid,
                           returning unmodified BranchStructs object ID ", obj$id))
            return(obj)
        }
    } else {
        warning("validation turned off, returning unvalidated BranchStructs")
        if (check_property(newobj, "has_topology")) newobj = setGraph(newobj)
        return(newobj)
    }

}

#' @export

getTreestruct <- function(obj, treestruct, ...) {
    UseMethod("getTreestruct", obj)
}

#' @export

getTreestruct.BranchStructs <- function(obj, concat = TRUE) {
    if (concat) {
        # # have to remove other nested columns before unnesting
        # colclasses = sapply(obj$treestructs, class)
        # listcols = colclasses %in% "list"
        # exludecols = names(colclasses)[listcols] != "treestruct"
        # ret = obj$treestructs
        # if (length(excludecols) > 1) {
        #     # more nested cols than just treestruct
        #     excludecolnames = names(colclasses)[excludecols]
        #     ret = ret %>% select(-excludecols) %>% unnest()
        # } else {
        #     # just treestruct is nested
        #     ret = ret %>% unnest()
        # }
        return(obj$treestructs %>% tidyr::unnest(treestruct, names_repair = "universal"))
    } else {
        return(obj$treestructs$treestruct)
    }
}

#' @export

getTreestruct.default <- function(obj) {
    warning("Doesn't apply to this class")
    return(obj)
}

#' @export
getTreestructs <- function(obj, ...) {
    UseMethod("getTreestructs", obj)
}

#' @export
getTreestructs.BranchStructs <- function(obj, list_cols = T) {
    if (list_cols) {
        return(obj$treestructs)
    } else {
        return(obj$treestructs[-list_cols(obj$treestructs, return_names = F)])
    }
}

#' @export
setTips <- function(obj) {
    UseMethod("setTips", obj)
}

#' @export
setTips.BranchStructs <- function(obj) {
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, setTips.default)
    obj$tips_set = T
    return(obj)
}

#' @export
setTips.default <- function(ts) {
    # a branch is a tip if no other branch claims it as a parent
    ts$tip = ! ts[["internode_id"]] %in% ts[["parent_id"]]
    return(ts)
}

#' @export
getTips <- function(obj) {
    UseMethod("getTips", obj)
}

#' @export
getTips.BranchStructs <- function(obj, make_compatible = T) {
    if (make_compatible) ret = obj %>% make_compatible()
    else ret = getTreestruct(obj, concat = T)
    ret = ret %>% dplyr::filter(tip)
    return(ret)
}

#' @export
setGraph <- function(obj, ...) {
    UseMethod("setGraph", obj)
}


#' @title setGraph.BranchStructs
#' @description FUNCTION_DESCRIPTION
#' @param obj BranchStructs object
#' @return OUTPUT_DESCRIPTION
#' @details Invalid treestructs will have graph set to NA.  This means that all functions making use of graphs will have to
#' be able to handle NA as an input in addition to a graph object.  This is a design decision.  An alternative is to use an empty graph object, but then
#' it becomes less obvious which treestructs are invalid and hence don't have graphs, and could lead to including erroneous
#' data in analyses.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname setGraph.BranchStructs
#'
setGraph.BranchStructs <- function(obj, ts_accessor = getTreestruct) {

    make_tidygraph <- function(ts, valid = T) {
        if (! valid) return(NA)
        thisedgelist = ts %>% dplyr::mutate( from = !!rlang::sym(obj$parentid_col),
                                             to = !!rlang::sym(obj$internodeid_col) ) %>%
            dplyr::select(from, to, everything())
        thisnodelist = ts %>% dplyr::mutate(id = !!rlang::sym(obj$internodeid_col), label = id) %>%
            dplyr::select(id, label, everything())
        allnodes = unique(c(thisedgelist$from, thisedgelist$to))
        missing_nodes = allnodes[ is.na(match(allnodes, thisnodelist$id)) ]
        thisnodelist = thisnodelist %>% dplyr::bind_rows(data.frame(id = missing_nodes))
        this_igraph = igraph::graph_from_data_frame(d = thisedgelist, vertices = thisnodelist, directed = TRUE)
        return(tidygraph::as_tbl_graph(this_igraph))
    }

    # do not construct graphs for treestructs that are invalid
    obj$treestructs = obj$treestructs %>%
        dplyr::rowwise() %>%
        dplyr::mutate(valid = all(dplyr::c_across(contains("valid")), na.rm = T))
            # mutate(graph = ifelse(valid, make_tidygraph(treestruct), NA)) #make_tidygraph(treestruct), NA))

    # set invalid treestruct graphs to NA
    obj$treestructs$graph = as.list(rep(NA, nrow(obj$treestructs)))
    obj$treestructs[obj$treestructs$valid,]$graph = map(ts_accessor(obj, concat = F)[obj$treestructs$valid], make_tidygraph)
    #obj$treestructs$graph = map(ts_accessor(obj, concat = F), make_tidygraph,  obj$treestructs$valid)
    obj = setNumComponents(obj)

    return(obj)
}

#' @export
setNumComponents <- function(obj) {
    UseMethod("setNumComponents", obj)
}

#' @export
setNumComponents.BranchStructs <- function(bss) {
    # num_components = NA if graph = NA.  Otherwise, count components.
    bss$treestructs$num_components = rep(NA_real_, nrow(bss$treestructs))

    bss$treestructs[!is.na(bss$treestructs$graph),]$num_components =
        purrr::map_dbl(bss$treestructs$graph[!is.na(bss$treestructs$graph)], igraph::count_components)

    return(bss)
}

#' @export
setNumComponents.default <- function(obj) {
    warning("method doesn't apply to this object, returning unmodified object")
    return(obj)
}

#' @export
getPerRadMetrics <- function(obj, ...) {
    UseMethod("getPerRadMetrics", obj)
}

#' @export
getPerRadMetrics.BranchStructs <- function(obj, idx = NA, concat = T) {
    if (! is.na(idx))
        return(getTreestructs(obj)$per_rad_metrics[[idx]])
    else if (concat) {
        # don't return other list columns
        return(obj$treestructs %>% tidyr::unnest(per_rad_metrics) %>% dplyr::select(-list_cols(., return_names = F)))
    } else {
        return(getTreestructs(obj)$per_rad_metrics)
    }
}

#' @export
getSummary <- function(obj) {
    UseMethod("getSummary", obj)
}

#' @export
getSummary.BranchStructs <- function(obj) {
    # remove nested (list) columns
    # listcols = sapply(hand_branches, function(x) {any(class(x) %in% "data.frame")})
    colclasses = sapply(obj$treestructs, typeof)
    is_listcol = colclasses %in% "list"
    listcols = names(colclasses)[!is_listcol]
    if (sum(is_listcol) > 0) {
        ret = obj$treestructs %>% dplyr::select(listcols)
    } else {
        ret = obj$treestructs
    }
    return(ret)
}

#' @export
check_property <- function(obj, ...) {
    UseMethod("check_property", obj)
}

#' @export
check_property.default <- function(obj, prop) {
    if (! prop %in% names(obj)) return(F)
    return(obj[[prop]])
}

#' @export
set_property <- function(obj, ...) {
    UseMethod("set_property", obj)
}

#' @export
set_property.default <- function(obj, prop, value) {
    obj[[prop]] = value
    return(obj)
}

# Validators ####

validate_treestruct <- function(obj) {
    UseMethod("validate_treestruct", obj)
}

#' @import crayon
validate_treestruct.BranchStructs <- function(obj) {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE

    # set validation flag columns
    obj$treestructs$valid_parents = NA
    obj$treestructs$valid_df = NA
    obj$treestructs$valid_internodes = NA

    all_valid = T
    treestructs = obj$treestructs
    for (this_row in 1:nrow(treestructs)) { # loop over all treestruct dfs in object

        thisbranch = treestructs[[this_row, obj$idcol]]
        this_treestruct = treestructs[[this_row, "treestruct"]][[1]]  # this is a list for hand-measured branches, so must add [[1]].
        if (verbose) message(paste("Validating branch",thisbranch))

        parents_valid = validate_parents(this_treestruct[[obj$internodeid_col]], this_treestruct[[obj$parentid_col]],
                                         ignore_errors = this_treestruct[[obj$ignore_error_col]])
        obj$treestructs[this_row,]$valid_parents = parents_valid
        if (!parents_valid) warning(crayon::red(paste("Parents NOT valid:",thisbranch)))

        df_valid = is.data.frame(this_treestruct)
        obj$treestructs[this_row,]$valid_df = df_valid
        if (!df_valid) warning(crayon::red(paste("treestruct not a dataframe:",thisbranch)))

        internodes_valid = validate_internodes(this_treestruct, obj$internodeid_col,
                                            ignore_error_col = obj$ignore_error_col)
        obj$treestructs[this_row,]$valid_internodes = internodes_valid
        if (!internodes_valid) warning(crayon::red(paste("Internodes NOT valid:",thisbranch)))

        this_valid = parents_valid & df_valid & internodes_valid
        all_valid = all_valid & this_valid
    }
    # assume columns are all there.  TODO validate column names.
    return(list(all_valid = all_valid, branchstructs = obj))
}

validate_treestruct.default <- function(...) {
    valid = treestruct::validate_parents(...)
}

#' @export
validate_connectivity <- function(obj) {
    UseMethod("validate_connectivity", obj)
}

#' @export
validate_connectivity.BranchStructs <- function(bss) {
    if (! "graph" %in% names(bss$treestructs)) {
        warning("setting Graph column for validation.  Run setGraph to keep the graph column in your object.")
        bss = SetGraph(bss)
    }
    valid = ! ( any(bss$treestructs$num_components != 1, na.rm = T) | any(is.na(bss$treestructs$num_components)) )
    if (!valid) warning("Multiple components in one or more treestructs")
    return(valid)
}

#' @export
validate_connectivity.default <- function(obj) {
    warning("method not applicable to this type")
    return(obj)
}

#' @export

reorder_internodes <- function(obj) {
    UseMethod("reorder_internodes", obj)
}

#' @title reorder_internodes.BranchStructs
#' @description reorders internodes in treestruct objects so pathlength and other accumulation routines work properly
#' @param obj TreeStructs object
#' @return TreeStructs object
#' @details This reorders all rows in the obj$treestructs$treestruct dataframes. It runs setTips if that hasn't already
#' been run.  This reorders rows such that, as you go down the dataframe, you can accumulate metrics (such as pathlength and
#' surface area).  That is, you are guaranteed that internodes above you (children, grandchildren, etc) will have been
#' assigned values.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname reorder_internodes.BranchStructs
reorder_internodes.BranchStructs <- function(obj) {
    if (! check_property(obj, "tips_set")) obj = setTips(obj)
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), reorder_internodes.default)
    obj$internodes_reordered = T
    return(obj)
}

#' @export
reorder_internodes.default <- function(ts) {
    tscopy = ts %>% dplyr::select(internode_id, parent_id, tip)
    ntips = sum(tscopy$tip)
    tscopy$orig_row = 1:nrow(ts)
    tscopy$parent_row = match(tscopy$parent_id, tscopy$internode_id)
    # which branch does this cylinder belong to
    tscopy$branchnum = NA_integer_
    # the sequence order of the cylinders within the branch
    tscopy$cyl_order_in_branch = NA_integer_
    # initiate branches from tip downwards
    tscopy[tscopy$tip,]$branchnum = 1:ntips
    tscopy[tscopy$tip,]$cyl_order_in_branch = 1
    # tscopy$daughter_row = match(tscopy[[obj$internodeid_col]], tscopy[[obj$parentid_col]])
    # follow each tip downwards, stopping when you hit the base or a cylinder already claimed by another branch
    for (thisbranch in 1:ntips) {
        curr_row = match(thisbranch, tscopy$branchnum)
        thiscylorder = 2
        repeat {
            next_row = tscopy[curr_row,]$parent_row
            if (is.na(next_row) |                           # we've hit the base
                ! is.na(tscopy[next_row,]$branchnum)) break # cylinder already claimed by another branch
            tscopy[next_row,]$branchnum = thisbranch
            tscopy[next_row,]$cyl_order_in_branch = thiscylorder
            thiscylorder = thiscylorder + 1
            curr_row = next_row
        }
    }
    # now reorder ts from last branch to first, and tip downwards within branch
    tscopy = tscopy %>% dplyr::arrange(desc(branchnum), cyl_order_in_branch)
    ts = ts[tscopy$orig_row,]    # use copied dataframe to reorder original (we've added columns to tscopy, etc)

    if ("parent_row" %in% names(ts)) {
        ts$parent_row = match(ts$parent_id, ts$internode_id) # reset parent_row now that we've moved things around.
    }
    return(ts)
}

#' @export
correct_furcations <- function(obj) {
    UseMethod("correct_furcations", obj)
}

#' @export
correct_furcations.BranchStructs <- function(obj) {
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), correct_furcations.default)
    obj$furcations_corrected = T
    return(obj)
}

#' @export
correct_furcations.default <- function(ts) {
    return(
        ts %>%
        dplyr::left_join(ts %>% dplyr::select(parent_id), by = c("internode_id" = "parent_id")) %>% # how many rows call you "parent"
            dplyr::add_count(internode_id, name = "n_furcation") %>%
            dplyr::group_by(internode_id) %>%
            dplyr::filter(dplyr::row_number() == 1)
    )
}

#' @export
assign_branch_num <- function(obj) {
    UseMethod("assign_branch_num", obj)
}

#' @title assign_branch_num.BranchStructs
#' @description assign branch numbers according to internodes occurring between furcations
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
#' @rdname assign_branch_num.BranchStructs
#' @seealso
#'
assign_branch_num.BranchStructs <- function(obj) {
    if (! check_property(obj, "tips_set")) obj = setTips(obj)
    if (! check_property(obj, "internodes_reordered")) obj = reorder_internodes(obj)
    if (! check_property(obj, "furcations_corrected")) obj = correct_furcations(obj)

    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), assign_branch_num.default)
    obj$branchnums_assigned = T
    return(obj)
}

#' @export
assign_branch_num.default <- function(ts) {
    branchnum_ret = assign_branchnum_cpp(ts$n_furcation, ts$tip)
    ts$branchnum = branchnum_ret[["branchnum"]]
    ts$order_in_branch = branchnum_ret[["order_in_branch"]]
    # build branch topology when assigning branch numbers
    ts$parent_branchnum = ts$branchnum[match(ts$parent_id, ts$internode_id)]
    branch_lookup = ts %>% dplyr::select(branchnum, parent_branchnum) %>% dplyr::filter(branchnum != parent_branchnum)
    ts$parent_branchnum = branch_lookup$parent_branchnum[match(ts$branchnum, branch_lookup$branchnum)]

    return(ts)
}

# Housekeeping ####

#' @export
parse_id <- function(obj, ...) {
    UseMethod("parse_id", obj)
}

#' @title parse_id.BranchStructs
#' @description parse the plot, tree tag, and branch identifier out from the concatenated identifier
#' @param obj treestruct object
#' @param treetag_regex regex to match tree tag, Default: `[^.]+`
#' @param branch_regex regex to match branchcode, Default: `B\\d+(SH|S|H)`
#' @param treetag_sec_no which section of the string (separated by `[_-]`) contains the treetag? Default: 2
#' @param branch_sec_no which section of the string (separated by `[_-]`) contains the branch code? Default: 3
#' @param nobranchcode set to TRUE if there is no branch code in the identifier, Default: FALSE
#' @return treestruct object
#' @details This function assumes the id column (idcol) is in the form of "plotcode-treecode-branchcode"
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname parse_id.BranchStructs
parse_id.BranchStructs <- function(obj, treetag_regex = "[^.]+", branch_regex = "B\\d+(SH|S|H)", treetag_sec_no = 2, branch_sec_no = 3, nobranchcode = F, ...) {
    split_treecode <- function(x) {
        codes = stringr::str_split(x, "[-_]")
        codes = purrr::map(codes, function(x) {
            x[treetag_sec_no] = stringr::str_extract(x[treetag_sec_no], treetag_regex)
            return(x)
        } )
        if (! nobranchcode) {
            codes = purrr::map(codes, function(x) {
                x[branch_sec_no] = stringr::str_extract(x[branch_sec_no], branch_regex)
                return(x)
            } )
        }
        return(codes)
    }

    assign_treecode_cols <- function(df, colname) {
        newcols = split_treecode(df[[colname]])
        df$plot = sapply(newcols, `[[`, 1)
        df$tag = sapply(newcols, `[[`, treetag_sec_no)
        if (! nobranchcode) {
            df$branch = sapply(newcols, `[[`, branch_sec_no)
        }
        return(df)
    }

    obj$treestructs = assign_treecode_cols(obj$treestructs, obj$idcol)

    return(obj)
}

#' @export
make_compatible <- function(obj) {
    UseMethod("make_compatible", obj)
}

#' @export
make_compatible.BranchStructs <- function(obj) {
    # note - returns dataframe, not object
    #TODO change column types to be compatible
    getTreestruct(obj, concat = T) %>% dplyr::select(branch_code:branch, internode_id:is_broken, d_child:tip, branchnum:pathlen)
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
#' @rdname calc_surfarea.BranchStructs

calc_surfarea.BranchStructs <- function(obj) {
    sa_fun <- function(x) {
        # TODO make colnames dynamic
        # surface area of a truncated cone
        x$surf_area = with(x, pi * (d_child/2 + d_parent/2) * sqrt((d_parent/2 - d_child/2)^2 + len^2))
        return(x)
    }
    obj$treestructs$treestruct = map(obj$treestructs$treestruct, sa_fun)
    obj$treestructs$surface_area_tot = map_dbl(obj$treestructs$treestruct, ~sum(.$surf_area, na.rm = T))
    return(obj)
}

#' @export
calc_sa_above <- function(obj) {
    UseMethod("calc_sa_above", obj)
}

#' @title Calculate surface area above every node
#'
#' @param obj Branchstruct object
#'
#' @return Branchstruct object
#' @details an 'sa_above' column is added to the treestruct nested dataframe
#' @export
#'
#' @examples
calc_sa_above.BranchStructs <- function(obj) {
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), calc_sa_above)
    return(obj)
}

#' @export
calc_sa_above.default <- function(ts) {
    if (! "parent_row" %in% names(ts)) {
        parent_row = parent_row(ts$parent_id, ts$internode_id)
    } else {
        parent_row = ts$parent_row
    }
    # wrap cpp function above
    if (! validate_internode_order(parent_row)) {
        warning(paste("Bad cyl file order in tree ", ts$tree[1]))
        ts$sa_above = NA
        return(ts)
    }
    ts$sa_above = calc_sa_above_cpp(ts$surf_area, parent_row)
    return(ts)
}

#' @export
calc_vol <- function(obj) {
    UseMethod("calc_vol", obj)
}

#' @export
calc_vol.BranchStructs <- function(obj) {
    vol_fun <- function(x) {
        x$vol = with(x, pi * ((d_child/2 + d_parent/2)/2)^2 * len)
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

#' @export
calc_vol_above <- function(obj) {
    UseMethod("calc_vol_above", obj)
}

#' @title Calculate wood volume above every node
#'
#' @param obj Branchstruct object
#'
#' @return Branchstruct object
#' @details a 'vol_above' column is added to the treestruct nested dataframe
#' @export
#'
#' @examples
calc_vol_above.BranchStructs <- function(obj) {
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), calc_vol_above)
    return(obj)
}

#' @export
calc_vol_above.default <- function(ts) {
    if (! "parent_row" %in% names(ts)) {
        parent_row = parent_row(ts$parent_id, ts$internode_id)
    } else {
        parent_row = ts$parent_row
    }
    # wrap cpp function above
    if (! validate_internode_order(parent_row)) {
        warning(paste("Bad cyl file order in tree ", ts$tree[1]))
        ts$vol_above = NA
        return(ts)
    }
    ts$vol_above = calc_total_x_above_internode_cpp(ts$vol, parent_row)
    return(ts)
}

#' @title calc_pathlen
#' @export
#' @rdname calc_pathlen
calc_pathlen <- function(obj) {
    UseMethod("calc_pathlen", obj)
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

calc_pathlen.BranchStructs <- function(bss) {
    # overwrite treestruct::calc_pathlen
    if (! "graph" %in% names(bss$treestructs)) {
        bss = setGraph(bss)
    }

    pathlen_vec <- function(treestructs) {
        # takes a 1-row treestructs dataframe
        # TODO maybe split out treetruct and treestructs operations into different functions....
        treestruct = treestructs$treestruct[[1]]
        parent_row = match(treestruct$parent_id, treestruct$internode_id)
        treestruct = treestruct %>% dplyr::mutate(len_tot = sum(len, na.rm = T))
        if (validate_internode_order(treestruct[[bss$parentid_col]], treestruct[[bss$internodeid_col]], parents_are_rows = F) &
            treestructs$num_components %in% 1) { # %in% instead of = since num_components is sometimes NA - in that case, don't calc pathlen
            # calc pathlength along entire network for each node
            if (treestructs$branch_code == "SAF05-T8-B1S") { browser() }
            treestruct = treestruct %>%
                dplyr::mutate(pathlen = calc_pathlen_cpp(treestruct$len, parent_row))
            # summary pathlength stats
            summary_pathlen_vec =
                treestruct %>% dplyr::filter(tip) %>%
                dplyr::summarize(pathlen_max = max(pathlen, na.rm = T),
                                 pathlen_min = min(pathlen, na.rm = T),
                                 pathlen_mean = mean(pathlen, na.rm = T),
                                 path_frac = pathlen_mean/pathlen_max )

            treestructs$treestruct = list(treestruct)
            treestructs = treestructs %>% dplyr::bind_cols(summary_pathlen_vec)
        } else {
            warning(paste("Invalid internode order or multiple components in", treestructs[[bss$idcol]], ". Just calculating total length."))
            treestruct$pathlen = NA
            treestructs$treestruct = list(treestruct)
            treestructs = treestructs %>%
                dplyr::mutate(pathlen_max = NA,
                              pathlen_min = NA,
                              pathlen_mean = NA,
                              path_frac = NA )
        }
        return(treestructs)
    }

    # bss$treestructs = bss$treestructs %>% rowwise() %>% do(pathlen_vec(.)) %>% bind_rows()
    # bss$treestructs = bss$treestructs %>% purrrlyr::by_row(pathlen_vec, .labels = F, .collate = c("cols")) # this messes with the column names
    bss$treestructs = bss$treestructs %>%
        dplyr::ungroup() %>%
        dplyr::mutate(tempcol = 1:dplyr::n()) %>%
        dplyr::group_by(tempcol) %>%
        #dplyr::rowwise() %>%
        dplyr::do(pathlen_vec(.)) %>%
        dplyr::bind_rows() %>%
        dplyr::ungroup() %>%
        dplyr::select(-tempcol) # rowwise turns "." into a list.  group_by(1:n()) doesn't...

    # calc total pathlength
    bss = calc_len(bss)

    return(bss)
}

#' @title calc_pathlen.default
#' @description Wrapper for C++ pathlength code
#' @param tree_structure a tree structure dataframe
#' @param length_col chr specifiying length column. default:"len"
#' @param parent_row_col chr specifiying parent row column. default:"parent_row"
#' @param path_len_col chr specifiying path length column. default:"path_len"
#' @return tree structure dataframe with a path_len column populated
#' @details DETAILS
#' @export
#' @rdname calc_pathlen.default
calc_pathlen.default <- function(tree_structure, length_col = "len",
                                 parent_row_col = "parent_row",
                                 path_len_col = "path_len") {
    # wrap cpp function
    pathlen_vec = calc_pathlen_cpp(tree_structure[[length_col]],
                                   tree_structure[[parent_row_col]])
    tree_structure[[path_len_col]] = pathlen_vec
    return(tree_structure)
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

#' @export
join_parents_to_treestruct <- function(obj) {
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = F), join_parents_to_treestruct.default)
    return(obj)
}

#' @export
join_parents_to_treestruct.default <- function(ts) {
    if (! "d_parent_internode" %in% names(ts)) {
        ts = ts %>%
            dplyr::left_join(
                ts %>% dplyr::select(c(internode_id, d_child, n_furcation)) %>%
                    dplyr::rename(d_parent_internode = d_child, n_furcation_parent = n_furcation),
                by = c("parent_id" = "internode_id"))
    } else {
        verbose <- getOption("treestruct_verbose")
        if(is.null(verbose)) verbose <- FALSE
        if (verbose) message("parents already joined to treestruct, skipping")
    }

    return(ts)
}

#' @export
radius_scaling <- function(obj) {
    UseMethod("radius_scaling", obj)
}

#' @export
radius_scaling.BranchStructs <- function(obj) {
    obj = join_parents_to_treestruct(obj)
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), radius_scaling.default)
    obj$treestructs$a_median = purrr::map_dbl(getTreestruct(obj, concat = FALSE), function(x) median(x$a, na.rm = T))
    return(obj)
}

#' @export
radius_scaling.default <- function(ts) {
    # add a radius scaling exponent to each row
    # n = n_child/n_parent (furcation)
    # beta = r_child/r_parent
    # a = - log(beta) / log(n)

    if ("d_child" %in% names(ts)) {
        ts = ts %>%
            dplyr::mutate(beta = ifelse(n_furcation_parent > 1, d_parent/d_parent_internode, NA),
                          a = ifelse(n_furcation_parent > 1, -log(beta)/log(n_furcation_parent), NA))
    } else {
        # this is working on cyl_summ from TreeStructs
        ts = ts %>%
            dplyr::mutate(beta = ifelse(n_furcation_parent > 1, rad/rad_parent, NA),
                          a = ifelse(n_furcation_parent > 1, -log(beta)/log(n_furcation_parent), NA))
    }
    return(ts)
}

#' @export
length_scaling <- function(obj) {
    UseMethod("length_scaling", obj)
}

#' @export
length_scaling.BranchStructs <- function(obj) {
    obj = join_parents_to_treestruct(obj)
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), length_scaling.default)
    obj$treestructs$b_median = purrr::map_dbl(getTreestruct(obj, concat = FALSE), function(x) median(x$b, na.rm = T))
    return(obj)
}

#' @export
length_scaling.default <- function(ts) {
    # add a length scaling exponent to each row
    # n = n_child/n_parent (furcation)
    # gamma = l_child/l_parent
    # b = - log(gamma) / log(n)

    if ("d_child" %in% names(ts)) {
        ts = ts %>%
            dplyr::left_join(
                ts %>% dplyr::select(c(internode_id, len)) %>% dplyr::rename(len_parent_internode = len),
                by = c("parent_id" = "internode_id")) %>%
            dplyr::mutate(gamma = ifelse(n_furcation_parent > 1, len/len_parent_internode, NA),
                          b = ifelse(n_furcation_parent > 1, -log(gamma)/log(n_furcation_parent), NA))
    } else {
        # this is working on cyl_summ from TreeStructs
        ts = ts %>%
            dplyr::mutate(gamma = ifelse(n_furcation_parent > 1, len/len_parent, NA),
                          b = ifelse(n_furcation_parent > 1, -log(gamma)/log(n_furcation_parent), NA))
    }
    return(ts)
}

#' @export
calc_per_rad_class_metrics <- function(obj, ...) {
    UseMethod("calc_per_rad_class_metrics", obj)
}

#' @export
calc_per_rad_class_metrics.BranchStructs <- function(obj, metrics = c("surf_area", "vol", "len"), bin_size = 0.005) {
    obj$treestructs$per_rad_metrics = purrr::map(getTreestruct(obj, concat = F), calc_per_rad_class_metrics.default, metrics = metrics, bin_size = bin_size)
    return(obj)
}

#' @export
calc_per_rad_class_metrics.default <- function(ts, metrics = c("surf_area", "vol", "len"), bin_size = 0.005) {

    perclass =
        ts %>%
        dplyr::ungroup() %>%
        dplyr::group_by(rad_class = as.numeric(as.character(cut(rad, breaks = seq(0, ceiling(max(rad, na.rm = T)*2*100)/200, by = bin_size),
                                                   labels = head(seq(0, ceiling(max(rad, na.rm = T)*2*100)/200, by = bin_size), -1))))) %>%
        dplyr::summarize_at(.vars = metrics, ~ sum(., na.rm = T))

    # fill classes with 0 if there aren't any branches in the smallest classes
    # hence guaranteed to have classes populated down to zero
    minrad = min(perclass$rad_class, na.rm = T)
    if (minrad > 0) {
        nrow = as.integer(minrad / bin_size)
        fillrows = data.frame(matrix(data = 0, nrow = nrow, ncol = length(metrics) + 1))
        names(fillrows) = c("rad_class", metrics)
        fillrows$rad_class = seq(minrad - bin_size, 0, -bin_size)
        perclass = rbind(perclass, fillrows)
    }

    return(perclass)
}

#' @export
run_all <- function(obj, ...) {
    UseMethod("run_all", obj)
}

#' @export
run_all.BranchStructs <- function(obj, calc_dbh = F, calc_summ_cyl = F, calc_max_height = F, make_graph_obj = T) {
    obj = run_all.default(obj, calc_dbh, calc_summ_cyl, calc_max_height, make_graph_obj)
    if (check_property(obj, "has_topology")) obj = calc_sa_above(obj)
    return(obj)
}

#' @export
run_all.default <- function(obj, calc_dbh = T, calc_summ_cyl = T, calc_max_height = T, make_graph_obj = T) {
    # if (! check_property(obj, "tips_set")) obj = setTips(obj)
    if (check_property(obj, "has_topology")) {
        obj = setTips(obj) # always set tips for now...  not necessary if reading from source...
        obj = calc_surfarea(obj)
        obj = calc_vol(obj)
        obj = calc_pathlen(obj)
        if (calc_max_height) obj = calc_max_height(obj)
        if (calc_dbh) obj = calc_dbh(obj)
        obj = correct_furcations(obj)
        obj = assign_branch_num(obj)
        if (calc_summ_cyl) obj = calc_summary_cyls(obj)
        obj = radius_scaling(obj)
        obj = length_scaling(obj)
        obj = calc_per_rad_class_metrics(obj)
        if (make_graph_obj) obj = setGraph(obj)
    } else {
        obj = calc_surfarea(obj)
        obj = calc_vol(obj)
        if (calc_max_height) obj = calc_max_height(obj)
        if (calc_dbh) obj = calc_dbh(obj)
        obj = calc_per_rad_class_metrics(obj)
    }
    return(obj)
}

# Utilities ####

#' @title save_cylfiles
#' @description Extract cylinder files from treestruct objects and save them to csv files.
#' @param obj PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname save_cylfiles

save_cylfiles <- function(obj, ...) {
    UseMethod("save_cylfiles", obj)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param path PARAM_DESCRIPTION, Default: '.'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname save_cylfiles.BranchStructs

save_cylfiles.BranchStructs <- function(obj, path = ".") {
    filenames = getFile(obj)
    ts_list = getTreestruct(obj, concat = F)
    for (i in 1:length(filenames)) {
        write.csv(ts_list[[i]], file.path(path, paste0(filenames[i],".csv")))
    }
}

#' @title truncate_branches
#' @description removes branches smaller than a given radius
#' @param obj PARAM_DESCRIPTION
#' @param rad_min minimum radius of branch to keep in meters.  0.05 would mean keeping all branches >= 5cm.
#' @return
#' @details make sure to re-run routines that calculate tree metrics after altering the tree
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname truncate_branches

truncate_branches <- function(obj, ...) {
    UseMethod("truncate_branches", obj)
}

#' @export
truncate_branches.BranchStructs <- function(obj, rad_min) {
    treestructs = getTreestructs(obj)
    obj$treestructs$treestruct = purrr::map(getTreestruct(obj, concat = FALSE), truncate_branches.default, rad_min = rad_min)
    obj = set_property(obj, "first_branch_assigned", FALSE) # need to re-find first branch when truncating, in case we've truncated first branch
    return(obj)
}

#' @export
truncate_branches.default <- function(ts, rad_min) {
    return(ts[ts$rad >= rad_min,])
}

#' @export
subset.Branchstructs <- function(obj, idx) {
    obj = obj %>%
        setTreestructs(getTreestructs(obj) %>%
                          dplyr::slice(idx))
    return(obj)
}

# Visualization ####

#' @export
visNetwork <- function(obj, ...) {
    UseMethod("visNetwork", obj)
}

#' @title visNetwork.BranchStructs
#' @description visualize branch network.  This overwrites visNetwork for BranchStructs objects.
#' @param bss BranchStructs object (requred)
#' @param index \code{character} or \code{integer} char matching BranchStructs id column, or index of branch in BranchStructs object (required)
#' @param hierarchical \code{logical}, plot in hierarchical layout. Default: T
#' @param width_factor \code{numeric} edge width factor, Default: 100
#' @param length_factor \code{numeric} edge length factor, Default: 10
#' @return visNetwork object
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  visNetwork(bss_obj) %>% visEdges(arrows = "middle")
#'  }
#' }
#' @export
#' @rdname visNetwork.BranchStructs
#' @seealso
#'  \code{\link[visNetwork]{visNetwork-igraph}},\code{\link[visNetwork]{visNetwork}}
visNetwork.BranchStructs <- function(bss, index, hierarchical = T, width_factor = 1000, length_factor = 10) {

    if (is.character(index)) {
        index = enquo(index)
        this_tidygraph = bss$treestructs %>%
            dplyr::filter(!!sym(bss$idcol) := !!index) %>%
            dplyr::pull(graph)
        this_id = bss$treestructs %>%
            dplyr::filter(!!sym(bss$idcol) := !!index) %>%
            dplyr::pull(!!sym(bss$idcol))
    } else if (is.numeric(index)) {
        # index is lookup number
        this_tidygraph = bss$treestructs[index,]$graph[[1]]
        this_id = bss$treestructs[index,bss$idcol][[1]]
    }
    netdata = visNetwork::toVisNetworkData(this_tidygraph)
    netdata$edges = netdata$edges %>% dplyr::mutate(width = d_parent * width_factor, length = len * length_factor)
    ret =
        visNetwork::visNetwork(netdata$nodes, netdata$edges, layout = "layout_with_fr", main = this_id) %>%
        visNetwork::visNodes(font = list(size = 25) ,
                             scaling = list(label = list(
                                 enabled = TRUE,
                                 min = 25, max = 100,
                                 maxVisible = 100,
                                 drawThreshold = 1
                             )))
    if (hierarchical) ret = ret %>% visNetwork::visHierarchicalLayout()
    return(ret)
}

#' @export
visNetwork.default <- function(obj) {
    return(visNetwork::visNetwork(obj))
}
