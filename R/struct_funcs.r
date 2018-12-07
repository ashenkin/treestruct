
# TODO import individual functions and not entire library when not necessary
# in particular: beware of using data.table and dplyr together.

#' @useDynLib treestruct, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarize mutate pull
#' @importFrom foreach "%dopar%"
#' @importFrom stats sd

# from https://github.com/STAT545-UBC/Discussion/issues/451
## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("."))

# Validation checks ####

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tree_structure PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname check_cylfile_internode_order

check_cylfile_internode_order <- function(tree_structure) {
    # Checks to make sure rows in cylfile are in the correct order for propagating
    #   cumulative calculations such as surface area above an internode
    good_order = all(tree_structure %>% mutate(rowdiff = parent_row - row_number(parent_row)) %>% pull(rowdiff) < 0)
    return(good_order)
}

#' @title validate_parents
#' @description check whether the parent/child structure tree structure data is valid (i.e., does every parent refer to an existing internode)
#' @param internode_ids vector
#' @param parent_ids vector
#' @param parents_are_rows If parent_ids are actually row numbers, not id's (common for QSM cyl_files), Default: F
#' @return TRUE if validation passes, FALSE if it fails
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname validate_parents

validate_parents <- function(internode_ids, parent_ids, parents_are_rows = F) {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE

    if(parents_are_rows) {
        parent_ids = internode_ids[parent_ids]
    }
    parents = match(parent_ids, internode_ids)
    valid = ifelse(sum(is.na(parents)) > 1, FALSE, TRUE) # only 1 NA allowed (should be the base)
    if (verbose & !valid) {
        warning(paste("parents don't exist:", paste(parent_ids[is.na(parents)], collapse = " ")))
    }

    return(valid)
}

validate_internodes <- function(treestruct_df, internode_col = "internode_id", parent_col = "parent_id") {
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE

    # check duplicated internodes...
    dups = duplicated(treestruct_df[[internode_col]])

    valid = !any(dups)
    if (verbose & !valid) {
        warning(paste("duplicated internodes",paste(treestruct_df[[internode_col]][dups], collapse = " ")))
    }
    return(valid)
}

# Structure Analysis ####

#' @title Astem_chambers_2004
#' @description Return surface area as estimate by Chambers' 2004 allometry
#' @param DBH vector; tree diameter at breast height (cm)
#' @return vector
#' @details See Chambers 2004 and GEM manual
#' @examples
#' Astem_chambers_2004(40)
#' @rdname Astem_chambers_2004
#' @export
Astem_chambers_2004 <- function(DBH) 10^(-0.105-0.686*log10(DBH)+2.208*(log10(DBH))^2-0.627*(log10(DBH))^3)

#' @title FUNCTION_TITLE
#' @description Calc surface areas of and above each internode
#' @param tree_structure PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname calc_sa_above

calc_sa_above <- function(tree_structure) {
    # wrap cpp function above
    if (! check_cylfile_internode_order(tree_structure)) {
        warning(paste("Bad cyl file order in tree ", tree_structure$tree[1]))
    }
    tree_structure$sa_above = calc_sa_above_cpp(tree_structure$surf_area, tree_structure$parent_row)
    return(tree_structure)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tree_structure PARAM_DESCRIPTION
#' @param start_order PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname pathlengths
#' @export
pathlengths <- function(tree_structure, start_order) {  # TODO implement start_order parameter
  # Calculate path length, return a tree_structure dataframe with a new path_len column
  # start_order defines the branch order to consider the "tip".  Defaults to the maximum branch order for each tip.
  # Start at each tip and go downwards to trunk, adding lengths along the way
  # Accepts a tree structre as created in the analyze_cyl_file function
  # First, find all tips
  tree_structure$path_len = NA
  tree_structure$tip = (tree_structure$daughter_row == 0)
  tree_structure = calc_pathlen(tree_structure)
  return(tree_structure)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param lidar_tree PARAM_DESCRIPTION
#' @param sapwood_depth PARAM_DESCRIPTION, Default: 2
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname sapwood_volume
#' @export
sapwood_volume <- function(lidar_tree, sapwood_depth = 2) {
  # Calcuate sapwood volume given a lidar tree with pathlengths already calculated
  # lidar_tree must have $tree_structure_full, $DBH, $POM, $mean_path_len fields
  # sapwood_depth is depth at DBH (1.3m).  sapwood_depth in cm, dbh in m, pom in m, path lengths in m.
  sapwood_depth = sapwood_depth/100 # work in meters
  mean_path_len = lidar_tree$mean_path_len
  r = lidar_tree$DBH/2
  pom = lidar_tree$POM
  sapwood_area = pi*r^2 - pi*(r - sapwood_depth)^2
  sapwood_vol_area_pres = sapwood_area * mean_path_len
  lidar_tree$sapwood_vol = list("area_pres" = sapwood_vol_area_pres)

  tree_structure = lidar_tree$tree_structure_full
  tree_structure$vol = pi*tree_structure$rad^2*tree_structure$len
  # const: constant sapwood depth until stem radii < sapwood depth, at which point all wood is considered sapwood.  No
  #   decisions about maintaining area per order, distributing sapwood, etc.
  tree_structure$heartwood_radius = with(tree_structure, rad - sapwood_depth)
  tree_structure$heartwood_radius[tree_structure$heartwood_radius < 0] = 0
  tree_structure$sw_vol_const = with(tree_structure, len * pi * (rad^2 - heartwood_radius^2))
  sw_vol_const_total = sum(tree_structure$sw_vol_const)

  sw_vol_const_per_order = tree_structure %>% group_by(branch_order) %>% summarize(sw_vol_const = sum(sw_vol_const)) %>%
    mutate(cum_trunk_upward = cumsum(sw_vol_const),
           cum_tips_downward = rev(cumsum(rev(sw_vol_const))))
  sw_vol_const_per_diam_class = tree_structure %>% group_by(size_class) %>% summarize(sw_vol_const = sum(sw_vol_const)) %>%
    mutate(cum_tips_downward = cumsum(sw_vol_const),
           cum_trunk_upward = rev(cumsum(rev(sw_vol_const))))
  # sw_vol_const_per_order = ddply(tree_structure, .(branch_order), summarize, sw_vol_const = sum(sw_vol_const))
  # sw_vol_const_per_diam_class = ddply(tree_structure, .(size_class), summarize, sw_vol_const = sum(sw_vol_const))

  lidar_tree$tree_structure_full = tree_structure
  lidar_tree$sapwood_vol[["sw_vol_const_total"]] = sw_vol_const_total
  lidar_tree$sapwood_vol[["sw_vol_const_per_order"]] = sw_vol_const_per_order
  lidar_tree$sapwood_vol[["sw_vol_const_per_diam_class"]] = sw_vol_const_per_diam_class

  return(lidar_tree)
}


# TODO: implement area-preserving (top-up and bottom-down) and dissipation-minimizing sapwood distribution algorithrms


#' @title analyze_cyl_file
#' @description Reads in a cyl file (which is usually created by treeqsm), and calculates path length, surface area, volume, and sapwood metrics.
#' @param cyl_file path to file with the cylinder data
#' @param calc_sapwood Whether or not sapwood should be calculated, Default: T
#' @param sapwood_depth Set assumed sapwood depth (cm) at the base of the tree, Default: 2
#' @param treeid identifier to use in outputs for this tree, Default: basename(cyl_file)
#' @param size_classes branch diameter size classes (cm) across which aggregated computations will be made, Default: c(0,1,3,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140)
#' @return a list of metrics and an updated tree structure dataframe with metrics appended to each row
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname analyze_cyl_file
#' @export
analyze_cyl_file <- function(cyl_file, calc_sapwood = T, sapwood_depth = 2, treeid, size_classes = c(0,1,3,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140)) {

    if (missing(treeid)) treeid = basename(cyl_file)

    # sapwood depth in cm
    sapwood_depth = sapwood_depth/100 # everything is in meters

    tree_structure <- utils::read.table(cyl_file,header=FALSE,sep="\t")

    if (ncol(tree_structure) == 14) {
        # older treeqsm versions
        colnames(tree_structure) <- c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                                      "daughter_row","branch_data_row","branch_order","index_num","added_after")
    } else {
        # newer treeqsm versions
        colnames(tree_structure) <- c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                                      "daughter_row","branch_data_row","branch_order","index_num","added_after","rad0")
    }

    tree_structure$tree = treeid

    tree_structure$z_corr = tree_structure$z_start - min(tree_structure$z_start) # start z at 0 if it doesn't already
    tree_structure$surf_area = 2*pi*tree_structure$rad * tree_structure$len

    surf_area_total = sum(tree_structure$surf_area)

    #surf_area_per_order = ddply(tree_structure, .(branch_order), summarize, surf_area = sum(surf_area))
    surf_area_per_order = tree_structure %>% group_by(branch_order) %>% summarize(surf_area = sum(surf_area))
    #surf_area_per_order = mutate(surf_area_per_order, cum_sa_trunk_upwards = cumsum(surf_area),
    #                             cum_sa_tips_downward = rev(cumsum(rev(surf_area))))
    surf_area_per_order = surf_area_per_order %>% mutate(cum_sa_trunk_upwards = cumsum(surf_area),
                                                         cum_sa_tips_downward = rev(cumsum(rev(surf_area))))

    surf_area_per_order$tree = treeid

    tree_structure$size_class = cut(tree_structure$rad*2*100, include.lowest = T, ordered_result = T,
                                    labels = size_classes[-length(size_classes)],
                                    breaks = size_classes)
    tree_structure$size_class = as.integer(as.character(tree_structure$size_class))
    surf_area_per_diam_class = tree_structure %>% group_by(size_class) %>% summarize(surf_area = sum(surf_area))
    #ddply(tree_structure, .(size_class), summarize, surf_area = sum(surf_area))

    # Calculate trunk-upwards and tips-downwards cumulative surface areas
    surf_area_per_diam_class = surf_area_per_diam_class %>% mutate(cum_sa_tips_downward = cumsum(surf_area),
                                                                   cum_sa_trunk_upwards = rev(cumsum(rev(surf_area))))

    surf_area_per_diam_class$tree = treeid

    # calculate path lengths
    tree_structure_full = pathlengths(tree_structure)

    # calculate surf_area supported by each internode  TODO: rename resulting df to something better
    tree_structure_full = calc_sa_above(tree_structure_full)

    # just get the row with z_start closest to 1.4 for now; you could get more sophisticated in the future if necessary
    # also make sure that you're getting a cylinder from the main stem (drooping branches or otherwise can cross the 1.4m line)
    main_stem = subset(tree_structure, branch_order == 0)
    dbh_row = which(abs(main_stem$z_corr - 1.4) == min(abs(main_stem$z_corr - 1.4)))
    DBH = main_stem[dbh_row,]$rad*2
    POM = main_stem[dbh_row,]$z_corr

    height = max(tree_structure$z_corr)

    # From GEM manual, DBH in cm, resulting surface area in m^2
    DBH_cm = DBH*100
    Astem_lidar_chambers_2004 = Astem_chambers_2004(DBH_cm)

    print(paste("about to return analyzed file", cyl_file))

    ret = list(file = basename(cyl_file),
               lidar_tree = treeid,
               surf_area_total = surf_area_total,
               surf_area_per_order = surf_area_per_order,
               surf_area_per_diam_class = surf_area_per_diam_class,
               Astem_lidar_chambers_2004 = c(treeid, Astem_lidar_chambers_2004),
               tree_structure_full = tree_structure_full,
               mean_path_len = mean(tree_structure_full$path_len, na.rm=T),
               DBH = DBH,
               POM = POM,
               height = height)

    if (calc_sapwood) {
        ret = sapwood_volume(ret, sapwood_depth)
    }

    return(ret)
}

#' @title Rename lidar tree
#' @description convenience function to name a QSM tree based on it's file name
#' @param name PARAM_DESCRIPTION
#' @param regex PARAM_DESCRIPTION
#' @param name_is_path PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @note This is only exported becuase it's used internally in foreach
#' @export
#' @rdname rename_tree

rename_tree <- function(name, regex, name_is_path = T) {
    # convenience function for renaming lidar_trees
    if (is.na(regex) | regex == "") {
        if (name_is_path) {
            return(basename(name))
        } else {
            return(name)
        }
    } else {
        return (sub(regex, "\\1", basename(name))) # always get rid of path
    }
}

#' @title process_qsm_dir
#' @description Runs analyze_cyl_file on all files in a directory
#' @param qsm_path PARAM_DESCRIPTION
#' @param parallel_process PARAM_DESCRIPTION, Default: T
#' @param cyl_file_pat PARAM_DESCRIPTION, Default: 'cyl.*.txt'
#' @param rename_pat Regex applied to filenames when naming ouput list elements (sub(rename_pat, '\\1', filename))
#' @param recursive PARAM_DESCRIPTION, Default: F
#' @param file_batching Number of files to process per batch, or 0 to process all at once, Default: 0
#' @return list, one element per cyl file, of outputs from analyze_cyl_file
#' @details This routine provides a facility to name each element of the list based on the name of the corresponding cyl file.  Use rename_pat parameter to customize the renaming.  File batching can be used to avoid memory errors, etc
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname process_qsm_dir
#' @export
process_qsm_dir <- function(qsm_path = ".", parallel_process = T,
                            cyl_file_pat = "cyl.*.txt", rename_pat = "",
                            recursive = F, file_batching = 0) {
    # rip through a directory, run analyze_cyl_file on all the qsm's
    # file_batching will process a number of files at a time, largely for debugging purposes.  Set to 0 for no batching.

    #TODO add progress bar
    lidar_trees = list()
    cyl_files = list.files(qsm_path, pattern = cyl_file_pat, full.names = T, recursive = recursive)
    if (parallel_process) {
        clust <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
        doParallel::registerDoParallel(clust)
        if (file_batching != 0) {
            cyl_batches = split(cyl_files, ceiling(seq_along(cyl_files)/file_batching))
            for (this_batch in cyl_batches) {
                print(paste("processing", this_batch))
                this_lidar_trees = foreach::foreach (i = 1:length(this_batch), .packages = c("treestruct", "dplyr")) %dopar% {
                    analyze_cyl_file(this_batch[i], treeid = rename_tree(i, rename_pat))
                }
                names(this_lidar_trees) = rename_tree(this_batch, rename_pat)
                lidar_trees = c(this_lidar_trees, lidar_trees)
            }
            rm(this_lidar_trees)
        } else {
            lidar_trees = foreach::foreach (i = 1:length(cyl_files), .packages = c("treestruct", "dplyr")) %dopar% {
                analyze_cyl_file(cyl_files[i], treeid = rename_tree(cyl_files[i], rename_pat))
            }
            names(lidar_trees) = rename_tree(cyl_files, rename_pat)
        }
        parallel::stopCluster(clust)
    } else {
        # non-parallel processing
        for (i in cyl_files) {
            print(i)
            lidar_trees[[rename_tree(i, rename_pat)]] = analyze_cyl_file(i, treeid = rename_tree(i, rename_pat))
        }
    }
    return(lidar_trees)
}


#' @title concat_lidar_trees
#' @description Takes a list of analyzed cyl files (such as would be produced by process_qsm_dir) and assembles dataframes useful for analysis and plotting.
#' @param lidar_trees PARAM_DESCRIPTION
#' @param pick_qsm PARAM_DESCRIPTION, Default: F
#' @param best_qsm PARAM_DESCRIPTION, Default: ''
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname concat_lidar_trees
#' @seealso \code{\link{concat_lidar_trees_parallel}}
#' @export
concat_lidar_trees <- function(lidar_trees, pick_qsm = F, best_qsm = "") {

    # Concatenate and further process analyzed QSMs
    # lidar_trees is a list, and must have the following fields populated in each element:
    # tree, qsm_idx, replicate, and qsm_setting

    # TODO redo with data.table. likely much quicker.  use data.table::rbindlist...

    first_run = T
    for (lidar_tree in names(lidar_trees)) {
        this_tree = lidar_trees[[lidar_tree]]
        if (first_run) {
            first_run = F
            surf_area_per_order_tot = this_tree[["surf_area_per_order"]]
            surf_area_per_order_tot$lidar_tree = lidar_tree
            surf_area_per_order_tot$tree = this_tree$tree
            surf_area_per_order_tot$qsm_idx = this_tree$qsm_idx
            surf_area_per_order_tot$replicate = this_tree$replicate
            surf_area_per_order_tot$qsm_setting = this_tree$qsm_setting

            surf_area_per_diam_class_tot = this_tree[["surf_area_per_diam_class"]]
            surf_area_per_diam_class_tot$lidar_tree = lidar_tree
            surf_area_per_diam_class_tot$tree = this_tree$tree
            surf_area_per_diam_class_tot$qsm_idx = this_tree$qsm_idx
            surf_area_per_diam_class_tot$replicate = this_tree$replicate
            surf_area_per_diam_class_tot$qsm_setting = this_tree$qsm_setting

            Astem_lidar_chambers_2004_tot = data.frame("lidar_tree" = lidar_tree,
                                                       "tree" = this_tree$tree,
                                                       "Astem" = as.numeric(this_tree[["Astem_lidar_chambers_2004"]][2]),   # TODO fix this (??)
                                                       "qsm_idx" = this_tree$qsm_idx,
                                                       "replicate" = this_tree$replicate,
                                                       "qsm_setting" = this_tree$qsm_setting)

            mean_path_len_tot = data.frame("lidar_tree" = lidar_tree,
                                           "tree" = this_tree$tree,
                                           "mean_path_len" = this_tree$mean_path_len,
                                           "qsm_idx" = this_tree$qsm_idx,
                                           "replicate" = this_tree$replicate,
                                           "qsm_setting" = this_tree$qsm_setting)

            sw_vol = data.frame("lidar_tree" = lidar_tree,
                                "tree" = this_tree$tree,
                                "sw_vol_area_pres" = this_tree[["sapwood_vol"]]$area_pres,
                                "sw_vol_const_total" = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                                "qsm_idx" = this_tree$qsm_idx,
                                "replicate" = this_tree$replicate,
                                "qsm_setting" = this_tree$qsm_setting)

            sw_vol_const_per_order_tot = this_tree[["sapwood_vol"]][["sw_vol_const_per_order"]]
            sw_vol_const_per_order_tot$lidar_tree = lidar_tree
            sw_vol_const_per_order_tot$tree = this_tree$tree
            sw_vol_const_per_order_tot$qsm_idx = this_tree$qsm_idx
            sw_vol_const_per_order_tot$replicate = this_tree$replicate
            sw_vol_const_per_order_tot$qsm_setting = this_tree$qsm_setting

            sw_vol_const_per_diam_class_tot = this_tree[["sapwood_vol"]][["sw_vol_const_per_diam_class"]]
            sw_vol_const_per_diam_class_tot$lidar_tree = lidar_tree
            sw_vol_const_per_diam_class_tot$tree = this_tree$tree
            sw_vol_const_per_diam_class_tot$qsm_idx = this_tree$qsm_idx
            sw_vol_const_per_diam_class_tot$replicate = this_tree$replicate
            sw_vol_const_per_diam_class_tot$qsm_setting = this_tree$qsm_setting

            tree_struct_full = this_tree[["tree_structure_full"]]

            tree_data = data.frame(tag = this_tree$tree,
                                   lidar_tree = lidar_tree,
                                   tree = this_tree$tree,
                                   replicate = this_tree$replicate,
                                   qsm_idx = this_tree$qsm_idx,
                                   qsm_setting = this_tree$qsm_setting,
                                   dbh = this_tree$DBH,
                                   pom = this_tree$POM,
                                   height = this_tree$height,
                                   mean_path_len = this_tree$mean_path_len,
                                   sw_vol_area_pres = this_tree[["sapwood_vol"]]$area_pres,
                                   sw_vol_const_total = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                                   Astem_chambers_2004 = this_tree$Astem_lidar_chambers_2004[2])

        } else {

            this_surf_area = this_tree[["surf_area_per_order"]]
            this_surf_area$lidar_tree = lidar_tree
            this_surf_area$tree = this_tree$tree
            this_surf_area$qsm_idx = this_tree$qsm_idx
            this_surf_area$replicate = this_tree$replicate
            this_surf_area$qsm_setting = this_tree$qsm_setting
            surf_area_per_order_tot = rbind(surf_area_per_order_tot, this_surf_area)

            this_surf_area = this_tree[["surf_area_per_diam_class"]]
            this_surf_area$lidar_tree = lidar_tree
            this_surf_area$tree = this_tree$tree
            this_surf_area$qsm_idx = this_tree$qsm_idx
            this_surf_area$replicate = this_tree$replicate
            this_surf_area$qsm_setting = this_tree$qsm_setting
            surf_area_per_diam_class_tot = rbind(surf_area_per_diam_class_tot, this_surf_area)

            this_Astem = data.frame("lidar_tree" = lidar_tree,
                                    "tree" = this_tree$tree,
                                    "Astem" = as.numeric(this_tree[["Astem_lidar_chambers_2004"]][2]),
                                    "qsm_idx" = this_tree$qsm_idx,
                                    "replicate" = this_tree$replicate,
                                    "qsm_setting" = this_tree$qsm_setting)
            Astem_lidar_chambers_2004_tot = rbind(Astem_lidar_chambers_2004_tot, this_Astem)

            this_mean_path_len = data.frame("lidar_tree" = lidar_tree,
                                            "tree" = this_tree$tree,
                                            "mean_path_len" = this_tree$mean_path_len,
                                            "qsm_idx" = this_tree$qsm_idx,
                                            "replicate" = this_tree$replicate,
                                            "qsm_setting" = this_tree$qsm_setting)
            mean_path_len_tot = rbind(mean_path_len_tot, this_mean_path_len)

            this_sw_vol_area_pres = data.frame("lidar_tree" = lidar_tree,
                                               "tree" = this_tree$tree,
                                               "sw_vol_area_pres" = this_tree[["sapwood_vol"]]$area_pres,
                                               "sw_vol_const_total" = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                                               "qsm_idx" = this_tree$qsm_idx,
                                               "replicate" = this_tree$replicate,
                                               "qsm_setting" = this_tree$qsm_setting)
            sw_vol = rbind(sw_vol, this_sw_vol_area_pres)

            this_sw_vol_const = this_tree[["sapwood_vol"]][["sw_vol_const_per_order"]]
            this_sw_vol_const$lidar_tree = lidar_tree
            this_sw_vol_const$tree = this_tree$tree
            this_sw_vol_const$qsm_idx = this_tree$qsm_idx
            this_sw_vol_const$replicate = this_tree$replicate
            this_sw_vol_const$qsm_setting = this_tree$qsm_setting
            sw_vol_const_per_order_tot = rbind(sw_vol_const_per_order_tot, this_sw_vol_const)

            this_sw_vol_const = this_tree[["sapwood_vol"]][["sw_vol_const_per_diam_class"]]
            this_sw_vol_const$lidar_tree = lidar_tree
            this_sw_vol_const$tree = this_tree$tree
            this_sw_vol_const$qsm_idx = this_tree$qsm_idx
            this_sw_vol_const$replicate = this_tree$replicate
            this_sw_vol_const$qsm_setting = this_tree$qsm_setting
            sw_vol_const_per_diam_class_tot = rbind(sw_vol_const_per_diam_class_tot, this_sw_vol_const)

            tree_struct_full = rbind(tree_struct_full, this_tree[["tree_structure_full"]])

            this_tree_data = data.frame(tag = this_tree$tree,
                                   lidar_tree = lidar_tree,
                                   tree = this_tree$tree,
                                   replicate = this_tree$replicate,
                                   qsm_idx = this_tree$qsm_idx,
                                   qsm_setting = this_tree$qsm_setting,
                                   dbh = this_tree$DBH,
                                   pom = this_tree$POM,
                                   height = this_tree$height,
                                   mean_path_len = this_tree$mean_path_len,
                                   sw_vol_area_pres = this_tree[["sapwood_vol"]]$area_pres,
                                   sw_vol_const_total = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                                   Astem_chambers_2004 = this_tree$Astem_lidar_chambers_2004[2])

            tree_data = rbind(tree_data, this_tree_data)

        }
    }

    surf_area_per_order_tot_avg = surf_area_per_order_tot %>% group_by(tree, branch_order) %>%
        summarize(cum_sa_trunk_upwards_mean = mean(cum_sa_trunk_upwards),
                  se = sd(cum_sa_trunk_upwards, na.rm=T)/sqrt(length(cum_sa_trunk_upwards)))


    surf_area_per_diam_class_tot_avg = surf_area_per_diam_class_tot %>% group_by(tree, size_class) %>%
        summarize(cum_sa_trunk_upwards_mean = mean(cum_sa_trunk_upwards),
                  se = sd(cum_sa_trunk_upwards, na.rm=T)/sqrt(length(cum_sa_trunk_upwards)))


    Astem_lidar_chambers_2004_tot_avg = Astem_lidar_chambers_2004_tot %>% group_by(tree) %>%
        summarize(Astem_mean = mean(Astem), se = sd(Astem, na.rm=T)/sqrt(length(Astem)))

    # this is where hand measurements are joined to field data.  will have to add this somewhere at some point ####
    # path_len_join = join(mean_path_len_tot, planchon_sa_tot, by = "tree")
    # sw_vol_join = join(sw_vol, planchon_sa_tot, by = "tree", type = "left")
    # sw_vol_join_avg = ddply(ifelse(pick_qsm, subset(sw_vol_join, qsm_setting == best_qsm), sw_vol_join), .(tree), summarize, sw_vol_area_pres_mean = mean(sw_vol_area_pres),
    #                         se = sd(sw_vol_area_pres, na.rm = T)/sqrt(length(sw_vol_area_pres)))

    return(list("surf_area_per_order_tot" = surf_area_per_order_tot,
                "surf_area_per_diam_class_tot" = surf_area_per_diam_class_tot,
                "Astem_lidar_chambers_2004_tot" = Astem_lidar_chambers_2004_tot,
                "mean_path_len_tot" = mean_path_len_tot,
                "sw_vol" = sw_vol,
                "sw_vol_const_per_order_tot" = sw_vol_const_per_order_tot,
                "sw_vol_const_per_diam_class_tot" = sw_vol_const_per_diam_class_tot,
                "surf_area_per_order_tot_avg" = surf_area_per_order_tot_avg,
                "surf_area_per_diam_class_tot_avg" = surf_area_per_diam_class_tot_avg,
                "Astem_lidar_chambers_2004_tot_avg" = Astem_lidar_chambers_2004_tot_avg,
                "tree_struct_full" = tree_struct_full,
                "tree_data" = tree_data))

}

#' @title concat_lidar_trees_parallel
#' @description As concat_lidar_trees, but parallelized.
#' @param lidar_trees PARAM_DESCRIPTION
#' @param pick_qsm PARAM_DESCRIPTION, Default: F
#' @param best_qsm PARAM_DESCRIPTION, Default: ''
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname concat_lidar_trees_parallel
#' @seealso \code{\link{concat_lidar_trees}}
#'
concat_lidar_trees_parallel <- function(lidar_trees, pick_qsm = F, best_qsm = "") {

    clust <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
    doParallel::registerDoParallel(clust)

    # Concatenate and further process analyzed QSMs
    # lidar_trees is a list, and must have the following fields populated in each element:
    # tree, qsm_idx, replicate, and qsm_setting

    list_of_dfs = foreach::foreach (lidar_tree = names(lidar_trees)) %dopar% {

        this_tree = lidar_trees[[lidar_tree]]

        surf_area_per_order_tot = this_tree[["surf_area_per_order"]]
        surf_area_per_order_tot$lidar_tree = lidar_tree
        surf_area_per_order_tot$tree = this_tree$tree
        surf_area_per_order_tot$qsm_idx = this_tree$qsm_idx
        surf_area_per_order_tot$replicate = this_tree$replicate
        surf_area_per_order_tot$qsm_setting = this_tree$qsm_setting

        surf_area_per_diam_class_tot = this_tree[["surf_area_per_diam_class"]]
        surf_area_per_diam_class_tot$lidar_tree = lidar_tree
        surf_area_per_diam_class_tot$tree = this_tree$tree
        surf_area_per_diam_class_tot$qsm_idx = this_tree$qsm_idx
        surf_area_per_diam_class_tot$replicate = this_tree$replicate
        surf_area_per_diam_class_tot$qsm_setting = this_tree$qsm_setting

        Astem_lidar_chambers_2004_tot = data.frame("lidar_tree" = lidar_tree,
                                                   "tree" = this_tree$tree,
                                                   "Astem" = as.numeric(this_tree[["Astem_lidar_chambers_2004"]][2]),   # TODO fix this (??)
                                                   "qsm_idx" = this_tree$qsm_idx,
                                                   "replicate" = this_tree$replicate,
                                                   "qsm_setting" = this_tree$qsm_setting)

        mean_path_len_tot = data.frame("lidar_tree" = lidar_tree,
                                       "tree" = this_tree$tree,
                                       "mean_path_len" = this_tree$mean_path_len,
                                       "qsm_idx" = this_tree$qsm_idx,
                                       "replicate" = this_tree$replicate,
                                       "qsm_setting" = this_tree$qsm_setting)

        sw_vol = data.frame("lidar_tree" = lidar_tree,
                            "tree" = this_tree$tree,
                            "sw_vol_area_pres" = this_tree[["sapwood_vol"]]$area_pres,
                            "sw_vol_const_total" = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                            "qsm_idx" = this_tree$qsm_idx,
                            "replicate" = this_tree$replicate,
                            "qsm_setting" = this_tree$qsm_setting)

        sw_vol_const_per_order_tot = this_tree[["sapwood_vol"]][["sw_vol_const_per_order"]]
        sw_vol_const_per_order_tot$lidar_tree = lidar_tree
        sw_vol_const_per_order_tot$tree = this_tree$tree
        sw_vol_const_per_order_tot$qsm_idx = this_tree$qsm_idx
        sw_vol_const_per_order_tot$replicate = this_tree$replicate
        sw_vol_const_per_order_tot$qsm_setting = this_tree$qsm_setting

        sw_vol_const_per_diam_class_tot = this_tree[["sapwood_vol"]][["sw_vol_const_per_diam_class"]]
        sw_vol_const_per_diam_class_tot$lidar_tree = lidar_tree
        sw_vol_const_per_diam_class_tot$tree = this_tree$tree
        sw_vol_const_per_diam_class_tot$qsm_idx = this_tree$qsm_idx
        sw_vol_const_per_diam_class_tot$replicate = this_tree$replicate
        sw_vol_const_per_diam_class_tot$qsm_setting = this_tree$qsm_setting

        tree_struct_full = this_tree[["tree_structure_full"]]

        tree_data = data.frame(tag = this_tree$tree,
                               lidar_tree = lidar_tree,
                               tree = this_tree$tree,
                               replicate = this_tree$replicate,
                               qsm_idx = this_tree$qsm_idx,
                               qsm_setting = this_tree$qsm_setting,
                               dbh = this_tree$DBH,
                               pom = this_tree$POM,
                               height = this_tree$height,
                               mean_path_len = this_tree$mean_path_len,
                               sw_vol_area_pres = this_tree[["sapwood_vol"]]$area_pres,
                               sw_vol_const_total = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                               Astem_chambers_2004 = this_tree$Astem_lidar_chambers_2004[2])

        return(list("surf_area_per_order_tot" = surf_area_per_order_tot,
                    "surf_area_per_diam_class_tot" = surf_area_per_diam_class_tot,
                    "Astem_lidar_chambers_2004_tot" = Astem_lidar_chambers_2004_tot,
                    "mean_path_len_tot" = mean_path_len_tot,
                    "sw_vol" = sw_vol,
                    "sw_vol_const_per_order_tot" = sw_vol_const_per_order_tot,
                    "sw_vol_const_per_diam_class_tot" = sw_vol_const_per_diam_class_tot,
                    "tree_struct_full" = tree_struct_full,
                    "tree_data" = tree_data))
    }

    parallel::stopCluster(clust)

    # rbind all dataframes in list together
    ret = list()
    for (this_df in names(list_of_dfs[[1]])) {
        ret[[this_df]] = data.table::rbindlist(lapply(list_of_dfs, `[[`, this_df))
    }

    ret[["surf_area_per_order_tot_avg"]] = ret[["surf_area_per_order_tot"]] %>% group_by(tree, branch_order) %>%
        summarize(cum_sa_trunk_upwards_mean = mean(cum_sa_trunk_upwards),
                  se = sd(cum_sa_trunk_upwards, na.rm=T)/sqrt(length(cum_sa_trunk_upwards)))


    ret[["surf_area_per_diam_class_tot_avg"]] = ret[["surf_area_per_diam_class_tot"]] %>% group_by(tree, size_class) %>%
        summarize(cum_sa_trunk_upwards_mean = mean(cum_sa_trunk_upwards),
                  se = sd(cum_sa_trunk_upwards, na.rm=T)/sqrt(length(cum_sa_trunk_upwards)))


    ret[["Astem_lidar_chambers_2004_tot_avg"]] = ret[["Astem_lidar_chambers_2004_tot"]] %>% group_by(tree) %>%
        summarize(Astem_mean = mean(Astem), se = sd(Astem, na.rm=T)/sqrt(length(Astem)))

    return(ret)

}

#' @title filter_opt_mods_df
#' @description keep only optimal models in a data.frame
#' @param df PARAM_DESCRIPTION
#' @param qsm_idx_col PARAM_DESCRIPTION, Default: 'qsm_idx'
#' @param opt_qsm_idxs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details If opt_qsm_idxs is empty, then the full dataframe is return unaltered
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname filter_opt_mods_df

filter_opt_mods_df <- function(df, qsm_idx_col = "qsm_idx", opt_qsm_idxs) {
    qsm_idx_col = quo(qsm_idx_col)
    if (is.na(opt_qsm_idxs) | length(opt_qsm_idxs) == 0) {
        return(df)
    } else {
        return(df %>% filter(!!qsm_idx_col %in% opt_qsm_idxs))
    }
}


#' @title add_opt_col
#' @description add "opt_set" and "opt_mod" columns to a dataframe, if the qsm_idx_col matches the opt_mod
#' @param lidar_trees concat_df[["dataset]]
#' @param dfs dataframes to add columns to, Default: 'all'
#' @param qsm_idx_col name of qsm index column within dataframes, Default: 'qsm_idx'
#' @param overwrite_existing_cols Should existing columns in the dataframe be overwritten those from opt_models?, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname add_opt_col

add_opt_col <- function(lidar_trees, dfs = "all", qsm_idx_col = "qsm_idx", overwrite_existing_cols = F) {

    if (! is.null(lidar_trees$all_optimal) && lidar_trees$all_optimal == T) {
        # handle special case of all trees being optimal
        return(add_all_opt_col(lidar_trees, dfs))
    }

    if (! "opt_models" %in% names(lidar_trees) | length(lidar_trees[["opt_models"]][["opt_mod"]]) == 0 ) {

        warning("optimal models data empty")
        return(lidar_trees)

    } else {
        if (dfs == "all") {
            # get all the elements of the list that are dataframes
            dfs = lidar_trees %>% map_lgl(is.data.frame)
            dfs = names(dfs)[which(dfs)]
            dfs = dfs[! dfs %in% "opt_models"]
        }
        opt_set_qsm_idxs = unlist(lidar_trees[["opt_models"]]$opt_set)
        opt_mods_qsm_idxs = unlist(lidar_trees[["opt_models"]]$opt_mod)
        opt_cols = c("opt_set", "opt_mod")

        for (this_df in dfs) {
            if ("data.frame" %in% class(lidar_trees[[this_df]]) & qsm_idx_col %in% colnames(lidar_trees[[this_df]])) {

                dup_cols = opt_cols[ opt_cols %in% names(lidar_trees[[this_df]]) ]

                if (length(dup_cols) > 0) {
                    if (overwrite_existing_cols) {
                        warning("Some columns to be added to", this_df, "already exist. Overwriting...")
                        lidar_trees[[this_df]] = lidar_trees[[this_df]] %>% select(-opt_cols)
                    } else {
                        warning("Some columns to be added to", this_df, "already exist. Skipping this element...")
                          next
                        }
                    }

                # see https://stackoverflow.com/questions/29678435/how-to-pass-dynamic-column-names-in-dplyr-into-custom-function
                # for using dynamic column names in dplyr
                lidar_trees[[this_df]] = lidar_trees[[this_df]] %>%
                                         mutate(opt_set = !!as.name(qsm_idx_col) %in% opt_set_qsm_idxs,
                                                opt_mod = !!as.name(qsm_idx_col) %in% opt_mods_qsm_idxs)
            } else {
                warning(paste(this_df, "either isn't a data.frame, or doesn't have column", qsm_idx_col))
            }
        }
    }
    return(lidar_trees)
}

# for internal use only, used when all the trees are optimal as signified by the "all_optimal" variable being true in the object
add_all_opt_col <- function(lidar_trees, dfs) {
    if (dfs == "all") {
        # get all the elements of the list that are dataframes
        dfs = lidar_trees %>% map_lgl(is.data.frame)
        dfs = names(dfs)[which(dfs)]
        dfs = dfs[! dfs %in% "opt_models"]
    }
    for (this_df in dfs) {
        lidar_trees[[this_df]]$opt_mod = T
        lidar_trees[[this_df]]$opt_set = T
    }
    return(lidar_trees)
}


#' @title add_tree_data
#' @description Adds data from the tree_data dataframe in a lidar_tree list/object to other dataframes in that list/object
#' @param lidar_trees PARAM_DESCRIPTION
#' @param dfs PARAM_DESCRIPTION, Default: 'all'
#' @param tree_data_df_name PARAM_DESCRIPTION, Default: 'tree_data'
#' @param data_cols PARAM_DESCRIPTION, Default: c("dbh", "pom", "height", "Astem_chambers_2004")
#' @param match_col PARAM_DESCRIPTION, Default: 'tree'
#' @param overwrite_existing_cols Should existing columns in the dataframe be overwritten those from tree_data?, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname add_tree_data

add_tree_data <- function(lidar_trees, dfs = "all", tree_data_df_name = "tree_data", data_cols = c("dbh", "pom", "height", "Astem_chambers_2004"), overwrite_existing_cols = F, match_col = "tree") {
    tree_data_df = lidar_trees[[tree_data_df_name]]

    # check that tree_data_df has what it needs
    if (! all(c(match_col, data_cols) %in% names(tree_data_df)) ) {
        warning(paste("Some columns missing from", tree_data_df_name, ".  Skipping adding tree data."))
        return(lidar_trees)
    }

    if (dfs == "all") {
        # get all the elements of the list that are dataframes
        dfs = lidar_trees %>% map_lgl(is.data.frame)
        dfs = names(dfs)[which(dfs)]
        dfs = dfs[! dfs %in% "opt_models"] # opt_models has a class dataframe for some reason, but shouldn't be processed here
    }

    for (this_df in dfs) {
        if ("data.frame" %in% class(lidar_trees[[this_df]]) & match_col %in% colnames(lidar_trees[[this_df]])) {
            dup_cols = data_cols[ data_cols %in% names(lidar_trees[[this_df]]) ]
            this_data_cols = data_cols

            if (length(dup_cols) > 0) {
                if (overwrite_existing_cols) {
                    warning("Some columns to be added to", this_df, "already exist. Overwriting...")
                    lidar_trees[[this_df]] = lidar_trees[[this_df]] %>% select(-dup_cols)
                } else {
                    warning("Some columns to be added to", this_df, "already exist. Skipping those columns...")
                    this_data_cols = data_cols[! data_cols %in% dup_cols]
                    if (length(this_data_cols == 0)) {
                        # all columns to be added aready exist, and we're not overwriting, so don't alter the dataframe
                        next
                    }
                }
            }

            lidar_trees[[this_df]] = lidar_trees[[this_df]] %>%
                left_join(select(tree_data_df, c(match_col, this_data_cols)), by = match_col)
        } else {
            warning(paste(this_df, "either isn't a data.frame, or doesn't have column", match_col))
        }
    }

    return(lidar_trees)

}
