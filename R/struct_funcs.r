library(dplyr)
library(Rcpp)

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

# functions ####
# from https://stackoverflow.com/questions/46368188/how-to-efficiently-calculate-path-lengths-in-tree-topology/46369457#46369457
Rcpp::cppFunction('
                  NumericVector calc_pathlen_cpp(NumericVector len, NumericVector idx){
                  int n = len.size();
                  NumericVector res(n);
                  double cumsum=0;
                  int j;
                  res[0] = len[0];

                  for(int i = 1; i < n; i++){
                  j = idx[ i ] - 1;
                  cumsum = res[ j ] + len[ i ];
                  res[ i ] = cumsum;
                  }
                  return res;
                  }' )

#' @title calc_pathlen
#' @description Wrapper for C++ pathlength code
#' @param tree_structure a tree structure dataframe
#' @return tree structure dataframe with a path_len column populated
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname calc_pathlen
#' @export
calc_pathlen <- function(tree_structure) {
  # wrap cpp function above
  pathlen_vec = calc_pathlen_cpp(tree_structure$len, tree_structure$parent_row)
  tree_structure$path_len = pathlen_vec
  return(tree_structure)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tree_structure PARAM_DESCRIPTION
#' @param start_order PARAM_DESCRIPTION
#' @param use_recurse_function PARAM_DESCRIPTION, Default: F
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
pathlengths <- function(tree_structure, start_order, use_recurse_function = F) {  # TODO implement start_order parameter
  # Calculate path length, return a tree_structure dataframe with a new path_len column
  # start_order defines the branch order to consider the "tip".  Defaults to the maximum branch order for each tip.
  # Start at each tip and go downwards to trunk, adding lengths along the way
  # Accepts a tree structre as created in the analyze_cyl_file function
  # First, find all tips
  tree_structure$path_len = NA
  tree_structure$tip = (tree_structure$daughter_row == 0)
  if (use_recurse_function) {
    for (tiprow in which(tree_structure$tip)) {
      # navigate from tip to trunk
      tree_structure[tiprow,]$path_len = recurse_pathlen(tree_structure, tiprow)
    }
  } else {
    tree_structure = calc_pathlen(tree_structure)
  }
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
  # lidar_tree must have $tree_structure_path_len, $DBH, $POM, $mean_path_len fields
  # sapwood_depth is depth at DBH (1.3m).  sapwood_depth in cm, dbh in m, pom in m, path lengths in m.
  sapwood_depth = sapwood_depth/100 # work in meters
  mean_path_len = lidar_tree$mean_path_len
  r = lidar_tree$DBH/2
  pom = lidar_tree$POM
  sapwood_area = pi*r^2 - pi*(r - sapwood_depth)^2
  sapwood_vol_area_pres = sapwood_area * mean_path_len
  lidar_tree$sapwood_vol = list("area_pres" = sapwood_vol_area_pres)

  tree_structure = lidar_tree$tree_structure_path_len
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

  lidar_tree$tree_structure = tree_structure
  lidar_tree$sapwood_vol[["sw_vol_const_total"]] = sw_vol_const_total
  lidar_tree$sapwood_vol[["sw_vol_const_per_order"]] = sw_vol_const_per_order
  lidar_tree$sapwood_vol[["sw_vol_const_per_diam_class"]] = sw_vol_const_per_diam_class

  return(lidar_tree)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cyl_file PARAM_DESCRIPTION
#' @param calc_sapwood PARAM_DESCRIPTION, Default: T
#' @param sapwood_depth PARAM_DESCRIPTION, Default: 2
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname analyze_cyl_file
#' @export
analyze_cyl_file <- function(cyl_file, calc_sapwood = T, sapwood_depth = 2) {
  # sapwood depth in cm
  library(plyr)
  library(dplyr)

  sapwood_depth = sapwood_depth/100 # everything is in meters

  tree_structure <-read.table(cyl_file,header=FALSE,sep="\t")
  colnames(tree_structure) <- c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                                "daughter_row","branch_data_row","branch_order","index_num","added_after")

  tree_structure$tree = basename(cyl_file)

  tree_structure$z_corr = tree_structure$z_start - min(tree_structure$z_start) # start z at 0 if it doesn't already
  tree_structure$surf_area = 2*pi*tree_structure$rad * tree_structure$len

  surf_area_total = sum(tree_structure$surf_area)

  surf_area_per_order = ddply(tree_structure, .(branch_order), summarize, surf_area = sum(surf_area))
  surf_area_per_order = mutate(surf_area_per_order, cum_sa_trunk_upwards = cumsum(surf_area),
                               cum_sa_tips_downward = rev(cumsum(rev(surf_area))))
  surf_area_per_order$tree = basename(cyl_file)

  tree_structure$size_class = cut(tree_structure$rad*2*100, include.lowest = T, ordered_result = T,
                                  labels = c(0, 1,2,5,10,20,30,40,50,60,70,80,90,100,110,120,130),
                                  breaks = c(0,1,2,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140))
  tree_structure$size_class = as.integer(as.character(tree_structure$size_class))
  surf_area_per_diam_class = ddply(tree_structure, .(size_class), summarize, surf_area = sum(surf_area))

  # Calculate trunk-upwards and tips-downwards cumulative surface areas
  surf_area_per_diam_class = mutate(surf_area_per_diam_class, cum_sa_tips_downward = cumsum(surf_area),
                                    cum_sa_trunk_upwards = rev(cumsum(rev(surf_area))))

  surf_area_per_diam_class$tree = basename(cyl_file)

  # calculate path lengths
  tree_structure_path_len = pathlengths(tree_structure)

  # just get the row with z_start closest to 1.4 for now; you could get more sophisticated in the future if necessary
  # also make sure that you're getting a cylinder from the main stem (drooping branches or otherwise can cross the 1.4m line)
  main_stem = subset(tree_structure, branch_order == 0)
  dbh_row = which(abs(main_stem$z_corr - 1.4) == min(abs(main_stem$z_corr - 1.4)))
  DBH = main_stem[dbh_row,]$rad*2
  POM = main_stem[dbh_row,]$z_corr

  height = max(tree_structure$z_corr)

  # From GEM manual, DBH in cm, resulting surface area in m^2
  DBH_cm = DBH*100
  Astem_lidar_chambers_2004 = 10^(-0.105-0.686*log10(DBH_cm)+2.208*(log10(DBH_cm))^2-0.627*(log10(DBH_cm))^3)

  ret = list(lidar_tree = basename(cyl_file),
             surf_area_total = surf_area_total,
             surf_area_per_order = surf_area_per_order,
             surf_area_per_diam_class = surf_area_per_diam_class,
             Astem_lidar_chambers_2004 = c(basename(cyl_file), Astem_lidar_chambers_2004),
             tree_structure_path_len = tree_structure_path_len,
             mean_path_len = mean(tree_structure_path_len$path_len, na.rm=T),
             DBH = DBH,
             POM = POM,
             height = height)

  if (calc_sapwood) {
    ret = sapwood_volume(ret)
  }

  return(ret)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param qsm_path PARAM_DESCRIPTION
#' @param parallel_process PARAM_DESCRIPTION, Default: T
#' @param cyl_file_pat PARAM_DESCRIPTION, Default: 'cyl.*.txt'
#' @param recursive PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname process_qsm_dir
#' @export
process_qsm_dir <- function(qsm_path, parallel_process = T, cyl_file_pat = "cyl.*.txt", recursive = F) {
  # rip through a directory, run analyze_cyl_file on all the qsm's

  #TODO add progress bar
  lidar_trees = list()
  cyl_files = list.files(qsm_path, pattern = cyl_file_pat, full.names = T, recursive = recursive)

  if (parallel_process) {
    clust <- makeCluster(detectCores() - 2)
    registerDoParallel(clust)
    lidar_trees = foreach (i = 1:length(cyl_files)) %dopar% {
      source("TLS_functions.r") # seems that dopar doesn't have access to the global environment when running, so need to source the functions to see analyze_cyl_file
      analyze_cyl_file(cyl_files[i])
    }
    names(lidar_trees) = basename(cyl_files)
    stopCluster(clust)
  } else {
    # non-parallel processing
    for (i in cyl_files) {
      lidar_trees[i] = analyze_cyl_file(i)
    }
  }
  return(lidar_trees)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
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
#' @export
concat_lidar_trees <- function(lidar_trees, pick_qsm = F, best_qsm = "") {

  # tree_regex = c("^(\\d+).*", "\\1"), trial_regex = c("^.*?-",""), rep_regex = c(".*-(\\d+).*","\\1"), qsm_regex = c("^.*?-(.*)-.*?$","\\1")) {

  # Concatenate and further process analyzed QSMs
  # lidar_trees is a list, and must have the following fields populated in each element:
  # tree, trial, replicate, and qsm_setting

  first_run = T
  for (lidar_tree in names(lidar_trees)) {
    this_tree = lidar_trees[[lidar_tree]]
    if (first_run) {
      first_run = F
      surf_area_per_order_tot = this_tree[["surf_area_per_order"]]
      surf_area_per_order_tot$lidar_tree = lidar_tree
      surf_area_per_order_tot$tree = this_tree$tree
      surf_area_per_order_tot$trial = this_tree$trial
      surf_area_per_order_tot$replicate = this_tree$replicate
      surf_area_per_order_tot$qsm_setting = this_tree$qsm_setting

      surf_area_per_diam_class_tot = lidar_trees[[lidar_tree]][["surf_area_per_diam_class"]]
      surf_area_per_diam_class_tot$lidar_tree = lidar_tree
      surf_area_per_diam_class_tot$tree = this_tree$tree
      surf_area_per_diam_class_tot$trial = this_tree$trial
      surf_area_per_diam_class_tot$replicate = this_tree$replicate
      surf_area_per_diam_class_tot$qsm_setting = this_tree$qsm_setting

      Astem_lidar_chambers_2004_tot = data.frame("lidar_tree" = lidar_tree,
                                                 "tree" = this_tree$tree,
                                                 "Astem" = as.numeric(this_tree[["Astem_lidar_chambers_2004"]][2]),   # TODO fix this (??)
                                                 "trial" = this_tree$trial,
                                                 "replicate" = this_tree$replicate,
                                                 "qsm_setting" = this_tree$qsm_setting)

      mean_path_len_tot = data.frame("lidar_tree" = lidar_tree,
                                     "tree" = this_tree$tree,
                                     "mean_path_len" = this_tree$mean_path_len,
                                     "trial" = this_tree$trial,
                                     "replicate" = this_tree$replicate,
                                     "qsm_setting" = this_tree$qsm_setting)

      sw_vol = data.frame("lidar_tree" = lidar_tree,
                          "tree" = this_tree$tree,
                          "sw_vol_area_pres" = this_tree[["sapwood_vol"]]$area_pres,
                          "sw_vol_const_total" = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                          "trial" = this_tree$trial,
                          "replicate" = this_tree$replicate,
                          "qsm_setting" = this_tree$qsm_setting)

      sw_vol_const_per_order_tot = this_tree[["sapwood_vol"]][["sw_vol_const_per_order"]]
      sw_vol_const_per_order_tot$lidar_tree = lidar_tree
      sw_vol_const_per_order_tot$tree = this_tree$tree
      sw_vol_const_per_order_tot$trial = this_tree$trial
      sw_vol_const_per_order_tot$replicate = this_tree$replicate
      sw_vol_const_per_order_tot$qsm_setting = this_tree$qsm_setting

      sw_vol_const_per_diam_class_tot = this_tree[["sapwood_vol"]][["sw_vol_const_per_diam_class"]]
      sw_vol_const_per_diam_class_tot$lidar_tree = lidar_tree
      sw_vol_const_per_diam_class_tot$tree = this_tree$tree
      sw_vol_const_per_diam_class_tot$trial = this_tree$trial
      sw_vol_const_per_diam_class_tot$replicate = this_tree$replicate
      sw_vol_const_per_diam_class_tot$qsm_setting = this_tree$qsm_setting

    } else {

      this_surf_area = this_tree[["surf_area_per_order"]]
      this_surf_area$lidar_tree = lidar_tree
      this_surf_area$tree = this_tree$tree
      this_surf_area$trial = this_tree$trial
      this_surf_area$replicate = this_tree$replicate
      this_surf_area$qsm_setting = this_tree$qsm_setting
      surf_area_per_order_tot = rbind(surf_area_per_order_tot, this_surf_area)

      this_surf_area = this_tree[["surf_area_per_diam_class"]]
      this_surf_area$lidar_tree = lidar_tree
      this_surf_area$tree = this_tree$tree
      this_surf_area$trial = this_tree$trial
      this_surf_area$replicate = this_tree$replicate
      this_surf_area$qsm_setting = this_tree$qsm_setting
      surf_area_per_diam_class_tot = rbind(surf_area_per_diam_class_tot, this_surf_area)

      this_Astem = data.frame("lidar_tree" = lidar_tree,
                              "tree" = this_tree$tree,
                              "Astem" = as.numeric(this_tree[["Astem_lidar_chambers_2004"]][2]),
                              "trial" = this_tree$trial,
                              "replicate" = this_tree$replicate,
                              "qsm_setting" = this_tree$qsm_setting)
      Astem_lidar_chambers_2004_tot = rbind(Astem_lidar_chambers_2004_tot, this_Astem)

      this_mean_path_len = data.frame("lidar_tree" = lidar_tree,
                                      "tree" = this_tree$tree,
                                      "mean_path_len" = this_tree$mean_path_len,
                                      "trial" = this_tree$trial,
                                      "replicate" = this_tree$replicate,
                                      "qsm_setting" = this_tree$qsm_setting)
      mean_path_len_tot = rbind(mean_path_len_tot, this_mean_path_len)

      this_sw_vol_area_pres = data.frame("lidar_tree" = lidar_tree,
                                         "tree" = this_tree$tree,
                                         "sw_vol_area_pres" = this_tree[["sapwood_vol"]]$area_pres,
                                         "sw_vol_const_total" = this_tree[["sapwood_vol"]]$sw_vol_const_total,
                                         "trial" = this_tree$trial,
                                         "replicate" = this_tree$replicate,
                                         "qsm_setting" = this_tree$qsm_setting)
      sw_vol = rbind(sw_vol, this_sw_vol_area_pres)

      this_sw_vol_const = this_tree[["sapwood_vol"]][["sw_vol_const_per_order"]]
      this_sw_vol_const$lidar_tree = lidar_tree
      this_sw_vol_const$tree = this_tree$tree
      this_sw_vol_const$trial = this_tree$trial
      this_sw_vol_const$replicate = this_tree$replicate
      this_sw_vol_const$qsm_setting = this_tree$qsm_setting
      sw_vol_const_per_order_tot = rbind(sw_vol_const_per_order_tot, this_sw_vol_const)

      this_sw_vol_const = this_tree[["sapwood_vol"]][["sw_vol_const_per_diam_class"]]
      this_sw_vol_const$lidar_tree = lidar_tree
      this_sw_vol_const$tree = this_tree$tree
      this_sw_vol_const$trial = this_tree$trial
      this_sw_vol_const$replicate = this_tree$replicate
      this_sw_vol_const$qsm_setting = this_tree$qsm_setting
      sw_vol_const_per_diam_class_tot = rbind(sw_vol_const_per_diam_class_tot, this_sw_vol_const)
    }
  }

  # surf_area_per_order_tot$trial = sub("^.*-","",surf_area_per_order_tot$tree)
  # surf_area_per_diam_class_tot$trial = sub("^.*-","",surf_area_per_diam_class_tot$tree)
  # Astem_lidar_chambers_2004_tot$trial = sub("^.*-","",Astem_lidar_chambers_2004_tot$tree)
  # mean_path_len_tot$trial = sub("^.*-","",mean_path_len_tot$tree)
  # sw_vol$trial = sub("^.*-","",sw_vol$tree)
  #
  # surf_area_per_order_tot$tree = str_extract(surf_area_per_order_tot$tree, "^M\\d+")
  # surf_area_per_diam_class_tot$tree = str_extract(surf_area_per_diam_class_tot$tree, "^M\\d+")
  # Astem_lidar_chambers_2004_tot$tree = str_extract(Astem_lidar_chambers_2004_tot$tree, "^M\\d+")
  # mean_path_len_tot$tree = str_extract(mean_path_len_tot$tree, "^M\\d+")
  # sw_vol$tree = str_extract(sw_vol$tree, "^M\\d+")

  surf_area_per_order_tot_avg = ddply(surf_area_per_order_tot, .(tree, branch_order), summarize, cum_sa_trunk_upwards_mean = mean(cum_sa_trunk_upwards),
                                      se = sd(cum_sa_trunk_upwards, na.rm=T)/sqrt(length(cum_sa_trunk_upwards)))


  surf_area_per_diam_class_tot_avg = ddply(surf_area_per_diam_class_tot, .(tree, size_class), summarize, cum_sa_trunk_upwards_mean = mean(cum_sa_trunk_upwards),
                                           se = sd(cum_sa_trunk_upwards, na.rm=T)/sqrt(length(cum_sa_trunk_upwards)))


  Astem_lidar_chambers_2004_tot_avg = ddply(Astem_lidar_chambers_2004_tot, .(tree), summarize, Astem_mean = mean(Astem),
                                            se = sd(Astem, na.rm=T)/sqrt(length(Astem)))

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
              "Astem_lidar_chambers_2004_tot_avg" = Astem_lidar_chambers_2004_tot_avg))

}

