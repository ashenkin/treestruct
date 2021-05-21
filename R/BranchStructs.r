#' @title BranchStructs
#' @description BranchStructs are hand-measured branches or trees.  This is a constructor for the BranchStructs S3 object.  Set column names before importing the treestruct dataframe
#' @param dataset name of dataset this object is based on
#' @param surface_area_total total surface area of branch
#' @param pathlen_total total pathlength of branch (sum of branch lengths)
#' @param treestructs dataframe of architecture data.  required.
#' @param idcol name of branch id column; default = "branch_code"
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
#' @details To create a branchstruct, pass a dataset with the following columns (unless you specify different column names, which is not recommended):
#' branch_code
#' date
#' internode_id
#' parent_id
#' n_furcation (optional)
#' len
#' major_child_angle (optional)
#' angle_internode_1 (optional)
#' angle_internode_2 (optional)
#' number_nodes_in_section (optional)
#' is_tip
#' is_broken
#' d_child
#' d_parent
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname BranchStructs

BranchStructs <- function(dataset = NA, treestructs, convert_to_meters = NA, has_topology = T) {
    # TODO allow setting of colnames

    stopifnot(length(dataset) == 1)

    branchDataset <-
        list(
            dataset = dataset,
            treestructs = NA,
            idcol = "branch_code",
            datecol = "date",
            internodeid_col = "internode_id",
            parentid_col = "parent_id",
            furcation_col = "n_furcation",
            length_col = "len",
            angle_col = "major_child_angle",
            angle_internode1_col = "angle_internode_1",
            angle_internode2_col = "angle_internode_2",
            numnodes_col = "number_nodes_in_section",
            istip_col = "is_tip",
            isbroken_col = "is_broken",
            d_child_col = "d_child", # top of cylinder (diameter at children)
            d_parent_col = "d_parent", # bottom of cylinder (diameter at parents)
            ignore_error_col = "ignore_error",
            radius_col = "rad", # computed from diameters, for compatibility with TreeStructs
            has_topology = has_topology,
            tips_set = F,
            internodes_reordered = F,
            furcations_corrected = F,
            branchnums_assigned = F
        )

    branchDataset = structure(branchDataset, class = "BranchStructs")
    branchDataset = setTreestruct(branchDataset, treestructs, convert_to_meters)
    if (! branchDataset$tips_set & branchDataset$has_topology) branchDataset = setTips(branchDataset)
    return(branchDataset)
}



