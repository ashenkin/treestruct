#' treestruct: A package for processing and analysing tree structural data
#'
#' The treestruct package handles tree data derived from QSM models as well as
#' from hand-measured data.  The key differences between the two are the
#' organizations of the corresponding data frames.  QSM tree data includes the following
#' columns:
#' Hand-measured data (usually from branches), in contrast, contains the following
#' columns:
#' ```
#' "branch_code", "date", "internode_id", "parent_id", "n_furcation", "len",
#' "major_child_angle", "angle_internode_1", "angle_internode_2", "number_nodes_in_section",
#' "is_tip", "is_broken", "d_child", "d_parent"
#' ```
#'
#' @section treestruct functions:
#' The functions ...
#'
#' @docType package
#' @name treestruct
NULL
