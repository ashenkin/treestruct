#' @title TreeStructs
#' @description Constructor for TreeStructs S3 object.  Set column names before importing the treestruct dataframe
#' @param dataset name of dataset this object is based on
#' @param treestructs dataframe of architecture data.  required.
#' @param surface_area_total  total surface area of branch
#' @param pathlen_totaltotal pathlength of branch (sum of branch lengths)
#' @param idcol name of Tree id column; default = "Tree_code"
#' @param datecol name of date column; default = "date"
#' @param internodeid_col name internode id column; default = "internode_id"
#' @param parentid_col name of parent internode id column; default = "parent_id"
#' @param furcation_col name of furcation number column; default = "n_furcation"
#' @param length_col name of internode length column; default = "len"
#' @param angle_col name of angle column; default = "major_child_angle"
#' @param angle_internode1_col name of angle internode 1 column; default = "angle_internode_1"
#' @param angle_internode2_col name of angle internode 2 column; default = "angle_internode_2"
#' @param numnodes_col name of number of nodes in section column; default = "number_nodes_in_section"
#' @param istip_col name of is tip column; default = "is_tip"
#' @param isbroken_col name of is broken column; default = "is_broken"
#' @param d_child_col name of internode diameter @ child (top) column; default = "d_child"
#' @param d_parent_col name of internode diameter @ parent (base) column; default = "d_parent"
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname TreeStructs

TreeStructs <- function(dataset = NA, treestructs) {
    # TODO allow setting of colnames

    stopifnot(length(dataset) == 1)

    TreeDataset <-
        list(
            dataset = dataset,
            treestructs = NA,
            idcol = "file",
            radius_col = "rad",
            length_col = "len",
            start_loc_x_col = "x_start",
            start_loc_y_col = "y_start",
            start_loc_z_col = "z_start",
            cyl_axis_x_col = "x_cyl",
            cyl_axis_y_col = "y_cyl",
            cyl_axis_z_col = "z_cyl",
            parent_row_col = "parent_row",
            daughter_row_col = "daughter_row",
            added_after_col = "added_after",
            unmod_rad_col = "unmod_rad",
            branch_data_row_col = "branch_data_row",
            branch_order_col = "branch_order",
            index_col = "index_num",
            tips_set = F
        )

    TreeDataset = structure(TreeDataset, class = "TreeStructs")
    TreeDataset = setTreestruct(TreeDataset, treestructs)
    return(TreeDataset)
}



