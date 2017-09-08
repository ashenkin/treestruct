
recurse_pathlen <- function(tree_structure, this_rownum){
  this_row = tree_structure[this_rownum,]
  # have we reached the ground?
  if (this_row$parent_row == 0) {
    # reached the ground
    return(this_row$len)
  } else {
    # accumulate lengths as we go back up out of the recursion
    return(this_row$len + recurse_pathlen(tree_structure, this_row$parent_row))
  }
}

pathlengths <- function(tree_structure, start_order) {  # TODO implement start_order parameter
  # Calculate mean path length, return a tree_structure dataframe with a new path_len column
  # start_order defines the branch order to consider the "tip".  Defaults to the maximum branch order for each tip.
  # Start at each tip and go downwards to trunk, adding lengths along the way
  # Accepts a tree structre as created in the analyze_cyl_file function
  # First, find all tips
  tree_structure$path_len = NA
  tree_structure$tip = (tree_structure$daughter_row == 0)
  for (tiprow in which(tree_structure$tip)) {
    # navigate from tip to trunk
    tree_structure[tiprow,]$path_len = recurse_pathlen(tree_structure, tiprow)
  }
  return(tree_structure)
}

sapwood_volume <- function(lidar_tree, sapwood_depth = 2) {
  # Calcuate sapwood volume given a lidar tree with pathlengths already calculated
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

analyze_cyl_file <- function(cyl_file, sapwood_depth = 2) {
  # sapwood depth in cm
  library("plyr")

  sapwood_depth = sapwood_depth/100 # everything is in meters

  tree_structure <-read.table(cyl_file,header=FALSE,sep="\t")
  colnames(tree_structure) <- c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                                "daughter_row","branch_data_row","branch_order","index_num","added_after")

  tree_structure$tree = cyl_file

  tree_structure$z_corr = tree_structure$z_start - min(tree_structure$z_start) # start z at 0 if it doesn't already
  tree_structure$surf_area = 2*pi*tree_structure$rad * tree_structure$len

  surf_area_total = sum(tree_structure$surf_area)

  surf_area_per_order = ddply(tree_structure, .(branch_order), summarize, surf_area = sum(surf_area))
  surf_area_per_order = mutate(surf_area_per_order, cum_sa_trunk_upwards = cumsum(surf_area),
                               cum_sa_tips_downward = rev(cumsum(rev(surf_area))))
  surf_area_per_order$tree = cyl_file

  tree_structure$size_class = cut(tree_structure$rad*2*100, include.lowest = T, ordered_result = T,
                                  labels = c(0, 1,2,5,10,20,30,40,50,60,70,80,90,100,110,120,130),
                                  breaks = c(0,1,2,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140))
  tree_structure$size_class = as.integer(as.character(tree_structure$size_class))
  surf_area_per_diam_class = ddply(tree_structure, .(size_class), summarize, surf_area = sum(2*pi*rad*len))

  # Calculate trunk-upwards and tips-downwards cumulative surface areas
  surf_area_per_diam_class = mutate(surf_area_per_diam_class, cum_sa_tips_downward = cumsum(surf_area),
                                    cum_sa_trunk_upwards = rev(cumsum(rev(surf_area))))

  surf_area_per_diam_class$tree = cyl_file

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

  return(list(lidar_tree = cyl_file,
              surf_area_total = surf_area_total,
              surf_area_per_order = surf_area_per_order,
              surf_area_per_diam_class = surf_area_per_diam_class,
              Astem_lidar_chambers_2004 = c(cyl_file, Astem_lidar_chambers_2004),
              tree_structure_path_len = tree_structure_path_len,
              mean_path_len = mean(tree_structure_path_len$path_len, na.rm=T),
              DBH = DBH,
              POM = POM,
              height = height))
}
