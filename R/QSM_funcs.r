# functions for interacting with QSM objects

#' @title readQSM
#' @description reads in a QSM file by dispatching to readQSM.mat or readQSM.treegraph
#' @param qsmfile path to file
#' @param qsmtype program used to create QSM - currently either QSMtree (.mat) or treegraph (.json). Default: "by_ext" (i.e., .mat = QSMtree, .json = treegraph)
#' @param qsmver version of treeQSM used to make the .mat file, Default: "by_name"
#' @return a 3-element list of file (filename), CylData, and BranchData
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname readQSM
#'
#'

readQSM <- function(qsmfile, qsmtype = c("by_ext", "QSMtree", "treegraph"), qsmver = "by_name") {
    qsmtype = match.arg(qsmtype)
    if (qsmtype == "by_ext") {
        require("tools")
        file_ext = tools::file_ext(qsmfile)
        if (file_ext == "mat") qsmtype = "QSMtree"
        else if (file_ext == "json") qsmtype = "treegraph"
        else stop(paste ("file extension", file_ext, "not recognized."))
    }
    if (qsmtype == "QSMtree") return(readQSM.mat(qsmfile, qsmver))
    if (qsmtype == "treegraph") return(readQSM.treegraph(qsmfile, qsmver))
}

#' @title readQSM.mat
#' @description reads in a QSM file in matlab format and passes back cylinder and branch data
#' @param qsmfile path to .mat file
#' @param qsmver version of treeQSM used to make the .mat file, Default: 2.3
#' @return a 3-element list of file (filename), CylData, and BranchData
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname readQSM.mat
#'
#'
readQSM.mat <- function(qsmfile, qsmver = "by_name") {
    if (!requireNamespace("rmatio", quietly = TRUE)) {
        stop("Package \"rmatio\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    QSMmat = rmatio::read.mat(qsmfile)
    # QSM files have 4 elements: treedata, inputs, qsm, and models

    #names(QSMmat[["qsm"]]) = c("cylinder", "branch", "treedata", "rundata", "pmdistance", "triangulation")  # names below 1st level in hierarchy dropped by R.matlab it seems

    cyldata = QSMmat[["qsm"]]$cylinder
    cyldata = cyldata %>% purrr::map(`[[`, 1) %>% as.data.frame() %>% dplyr::bind_cols()
    branchdata = QSMmat[["qsm"]]$branch
    branchdata = branchdata %>% purrr::map(`[[`, 1) %>% as.data.frame() %>% dplyr::bind_cols()

    cyl_col_names = c("rad" = "radius",
                      "len" = "length",
                      "x_start" = "start.1",
                      "y_start" = "start.2",
                      "z_start" = "start.3",
                      "x_cyl" = "axis.1",
                      "y_cyl" = "axis.2",
                      "z_cyl" = "axis.3",
                      "parent_row" = "parent",
                      "daughter_row" = "extension",
                      "added_after" = "added",
                      "branch_data_row" = "branch",
                      "branch_order" = "BranchOrder",
                      "index_num" = "PositionInBranch")

    branch_col_names = c("bord" = "order",
                         "bpar" = "parent",
                         "bvol" = "volume",
                         "blen" = "length",
                         "bang" = "angle",
                         "bheight" = "height",
                         "baz" = "azimuth",
                         "bdiam" = "diameter")

    # reorder and rename columns

    if (is.numeric(qsmver) & qsmver >= 2.3) {
        testit::assert("Length of cyl columns not right", length(names(cyldata)) == 15, length(names(branchdata)) == 5)
        names(cyldata) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row","branch_data_row","branch_order","index_num","added_after", "rad0")
        names(branchdata) = c("bord", "bpar", "bvol", "blen", "bang")
    } else if (qsmver == "UCL") {
        testit::assert("Length of cyl columns not right", length(names(cyldata)) == 15, length(names(branchdata)) == 8)
        names(cyldata) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row", "added_after", "unmod_rad", "branch_data_row","branch_order","index_num")
        names(branchdata) = c("bord", "bpar", "bvol", "blen", "bang", "bheight", "baz", "bdiam")
    } else if (is.numeric(qsmver) & qsmver < 2.3) {
        testit::assert("Length of cyl columns not right", length(names(cyldata)) == 14, length(names(branchdata)) == 5)
        names(cyldata) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row","branch_data_row","branch_order","index_num","added_after")
        names(branchdata) = c("bord", "bpar", "bvol", "blen", "bang")
    } else if (qsmver == "by_name") { # match by name
        cyldata = cyldata %>% dplyr::rename(!!cyl_col_names)
        branchdata = branchdata %>% dplyr::rename(!!branch_col_names)
    } else {
        stop("must specify a qsm version")
    }

    return(list(file = basename(qsmfile), CylData = cyldata, BranchData = branchdata))
}


#' @title readQSM.treegraph
#' @description reads in a QSM file in treegraph format and passes back cylinder and branch data
#' @param qsmfile path to treegraph file
#' @param qsmver version of treegraph used to make the treegraph file, Default: by_name
#' @return a 3-element list of file (filename), CylData, and BranchData
#' @details Treegraphs are described by nodes that connect cylinders, as opposed to internodes as in treeQSM.
#' In addition, treegraphs sometimes have furcations at the most basal node (think of a tree that branches at the soil
#' surface).  In these cases, treestruct's use of internodes instead of nodes as the most basic element breaks down,
#' as there is no way to indicate that the two (or more) most basal cylinders should be connected.  To get around this issue,
#' we add a 0-length, 0-width cylinder at the base of treegraph objects with basal furcations to enable that connectivity.  We must make
#' sure that analyses that could be thrown off by this extra cylinder make accomodations for it.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname readQSM.treegraph
#'
#'
readQSM.treegraph <- function(qsmfile, qsmver = "by_name") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package \"jsonlite\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    QSMgraph = jsonlite::fromJSON(qsmfile)
    # treegraph metadata:
        # cylinder data is in the "cyls" key
            # all units in meters
            # within cyls:
            #     p1 - end node number (same as internode number)
            # p2 - start node number (same as parent internode number)
            # sx/sy/sz - center of circular endcap at the beginning of the cylinder
            # ax/ay/az - unit vector of orientation of cylinder
            # radius
            # length
            # vol
            # surface_area
            # point_density
            # nbranch - number of branch this cyl has been assigned to
            # ninternode
            # ncyl
            # is_tip - logical, if cylinder is a tip
            # branch_order
            # is_tipI - integer (1/0), if cylinder is a tip

        # branch data is in the "centres" key
            # "index"
            # "centre_id"
            # "cx"
            # "cy"
            # "cz"
            # "dist2fur"
            # "distance_from_base"
            # "is_tip"
            # "n_furcation"
            # "n_points"
            # "nbranch"
            # "ncyl"
            # "ninternode"
            # "node_id"
            # "parent"
            # "parent_node"
            # "slice_id"
            # "sf_radius"
            # "sf_error"
            # "m_radius"

    unpack_json_cyl = function(x) {
        x[sapply(x, is.null)] <- NA  # NULL entries get dropped by unlist.  We don't want this.
        return( as.data.frame( unlist( unlist(x))))
    }

    cyldata = jsonlite::fromJSON(txt = QSMgraph$cyls)
    df_names = names(cyldata)
    cyldata = suppressMessages(purrr::map_dfc(cyldata, unpack_json_cyl))
    names(cyldata) = df_names

    branchdata = jsonlite::fromJSON(txt = QSMgraph$centres)
    df_names = names(branchdata)
    branchdata = suppressMessages(purrr::map_dfc(branchdata, unpack_json_cyl))
    names(branchdata) = df_names

    cyl_col_names = c("rad" = "radius",
                      "len" = "length",
                      "x_start" = "sx",
                      "y_start" = "sy",
                      "z_start" = "sz",
                      "x_cyl" = "ax",
                      "y_cyl" = "ay",
                      "z_cyl" = "az",
                      "parent_id" = "p2",
                      "internode_id" = "p1",
                      # this will be different from QSM in that it ID's the parent id, not the parent row.
                      #"parent_row" = "parent",
                      #"daughter_row" = "extension",
                      # "added_after" = NA,
                      "branch_data_row" = "nbranch", #??
                      "branch_order" = "branch_order",
                      "index_num" = "ncyl")

    # so many differences from QSM, not worth importing
    # branch_col_names = c("bord" = NA, # in cyl data
    #                      "bpar" = "parent",
    #                      "bvol" = NA, # in cyl data
    #                      "blen" = NA, # in cyl data
    #                      "bang" = "angle",
    #                      "bheight" = "height",
    #                      "baz" = "azimuth",
    #                      "bdiam" = "diameter")

    # reorder and rename columns

    if (qsmver == "by_name") { # match by name
        cyldata = cyldata %>%
            dplyr::rename(!!cyl_col_names) %>%
            dplyr::select(names(cyl_col_names)) # only keep columns that match what we need
        # add a 0-width, 0-length cylinder to the base of the QSM to enable connectivity for basal furcations
        # assume that an internode with NA parent indicates a basal internode.
        # we'll just add this extra cylinder if there is more than one basal internode (i.e. basal furcation)
        # TODO just add the basal node when it's really needed, not just when sum(NA) > 1...
        # had been assuming all na's were due to missing basal node, but turns out that's not the case.
        cyldata$parent_row = treestruct::parent_row(cyldata$parent_id, cyldata$internode_id) # set parent row if needed in other treestruct functions

        if (sum(is.na(cyldata$parent_row)) > 1) {
            base_cyl = cyldata %>% dplyr::filter(is.na(parent_row)) %>%
                slice(1L) %>%
                dplyr::mutate(rad = 0, len = 0, internode_id = parent_id, parent_id = NA, parent_row = NA, branch_data_row = NA, branch_order = NA, index_num = NA)
            cyldata = base_cyl %>% dplyr::bind_rows(cyldata)
            # re-run parent row with basal internode added
            cyldata$parent_row = treestruct::parent_row(cyldata$parent_id, cyldata$internode_id)
        }

        # branchdata = branchdata %>% dplyr::rename(!!branch_col_names)
    } else {
        stop("must specify a qsm version")
    }

    cyldata$parent_row = treestruct::parent_row(cyldata$parent_id, cyldata$internode_id) # set parent row if needed in other treestruct functions

    return(list(file = basename(qsmfile), CylData = cyldata, BranchData = NA))
}

# OLD VERSION - keeping here for now for easy reference.  not sure how this ever worked... strange mat files?
# readQSM.mat <- function(qsmfile, qsmver = 2.3) {
#     if (!requireNamespace("R.matlab", quietly = TRUE)) {
#         stop("Package \"R.matlab\" needed for this function to work. Please install it.",
#              call. = FALSE)
#     }
#     QSMmat = R.matlab::readMat(qsmfile)
#
#     CylData_cols = c("rad", "len", "sta", "axe", "cpar", "cext", "boc", "added", "unmodradius")
#     BranchData_cols = c("bord", "bpar", "bvol", "blen", "bang")
#
#     cyldata_ele2col_order = order(match(tolower(names(QSMmat)), CylData_cols), na.last = NA)
#     branchdata_ele2col_order = order(match(tolower(names(QSMmat)), BranchData_cols), na.last = NA)
#
#     #cydata_elements = grepl(paste0("^", paste(CylData_cols, collapse="$|^"), "$"), names(QSMmat), ignore.case = T)
#     #branchdata_elements = grepl(paste0("^", paste(BranchData_cols, collapse="$|^"), "$"), names(QSMmat), ignore.case = T)
#
#     CylData = data.frame(do.call(cbind, QSMmat[ cyldata_ele2col_order ]))
#     BranchData = data.frame(do.call(cbind, QSMmat[ branchdata_ele2col_order ]))
#
#     #names(CylData) =
#
#     # reorder and rename columns
#     names(BranchData) = c("bord", "bpar", "bvol", "blen", "bang")
#
#     if (qsmver >= 2.3) {
#         names(CylData) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
#                            "daughter_row","branch_data_row","branch_order","index_num","added_after", "rad0")
#     } else {
#         names(CylData) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
#                            "daughter_row","branch_data_row","branch_order","index_num","added_after")
#     }
#
#     return(list(file = basename(qsmfile), CylData = CylData, BranchData = BranchData))
# }

#' @title import_treestructs_from_dir
#' @description utility function to create treestruct dataframe that can be used to create a TreeStructs object
#' @param qsm_path PARAM_DESCRIPTION, Default: '.'
#' @param recursive PARAM_DESCRIPTION, Default: F
#' @param qsmver PARAM_DESCRIPTION, Default: "by_name"
#' @param qsmtype program used to create QSM - currently either QSMtree (.mat) or treegraph (.json). Default: "by_ext" (i.e., .mat = QSMtree, .json = treegraph)
#' @param filematch_pattern regex filematch pattern (used in list.files), Default: "(\\.mat)|(\\.json)"
#' @param nest nest cylinder dataframes in a column, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details this function extracts the CylData elements of all the "*.mat" files in the indicated directory
#' @examples
#' \dontrun{
#' if(interactive()){
#'  my_3d_trees = TreeStructs(dataset = "my_trees", treestructs = import_treestructs_from_dir("./qsm_files"))
#'  }
#' }
#' @export
#' @rdname import_treestructs_from_dir
#'
#' @import purrr
#' @import tidyr
import_treestructs_from_dir <- function(qsm_path = ".", recursive = F, qsmver = "by_name", qsmtype = c("by_ext", "QSMtree", "treegraph"), filematch_pattern = "(\\.mat)|(\\.json)", nest = F) {
    qsmtype = match.arg(qsmtype)
    verbose <- getOption("treestruct_verbose")
    if(is.null(verbose)) verbose <- FALSE
    qsm_mats = list.files(path = qsm_path,
                          pattern = filematch_pattern,
                          recursive = recursive, ignore.case = T, full.names = T)
    if (length(qsm_mats) == 0) {
        warning("No qsm files found in ", qsm_path)
        return(NA)
    }
    pb <- txtProgressBar(max = length(qsm_mats), style = 3)
    setTxtProgressBar(pb, 0)
    qsms = list()
    for (qsmfile in qsm_mats) {
        if (verbose) print(paste("Reading", qsmfile))
        # assume for now that .mat files are QSMtree, and .json files are treegraph
        qsms[[qsmfile]] = readQSM(qsmfile, qsmtype = qsmtype, qsmver = qsmver)
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    }

    close(pb)

    treestructs = qsms %>% purrr::map_dfr(~data.frame(file = .$file, .$CylData))

    if (nest) treestructs = treestructs %>% group_by(file) %>% tidyr::nest(.key = "cyldata")

    return(treestructs)
}

#' @title write_files_from_mats_dir
#' @description extracts cyl and branch data from .mat files made by treeQSM and writes txt files used by this package.
#' @param qsm_path path to directory containing treeQSM .mat files
#' @param recursive process sub-directories recursively?  Default: FALSE
#' @param qsmver version of treeQSM used to make the .mat file, Default: 2.3
#' @return vector of .mat filenames processed
#' @details If you get an error about names not being the same length as the vector, try a different qsm version.  This is probably slower than the equivalent server function, but also more reliable.  This functions does not use matlab at all.  Warning: this will overwrite existing files.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname write_files_from_mats_dir
#'
write_files_from_mats_dir <- function(qsm_path = ".", recursive = F, qsmver = 2.3) {

    qsm_mats = list.files(qsm_path, pattern = "\\.mat", full.names = T, recursive = recursive)
    pb <- txtProgressBar(max = length(qsm_mats), style = 3)
    setTxtProgressBar(pb, 0)
    for (this_qsm_mat in qsm_mats) {
        this_qsm_r = readQSM.mat(this_qsm_mat, qsmver = qsmver)
        this_dir = dirname(this_qsm_mat)
        this_file_basename = basename(tools::file_path_sans_ext(this_qsm_mat))
        cyl_filename = paste0("cyl_data_", this_file_basename, ".txt")
        branch_filename = paste0("branch_data_", this_file_basename, ".txt")
        write.table(this_qsm_r$CylData, file = file.path(this_dir, cyl_filename), quote = F, sep = "\t", col.names = F, row.names = F)
        write.table(this_qsm_r$BranchData, file = file.path(this_dir, branch_filename), quote = F, sep = "\t", col.names = F, row.names = F)
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    }
    close(pb)
    return(qsm_mats)
}


#' @title save_model_text_matlab_server
#' @description Experimental and not tested!  Uses connection to matlab server to process .mat QSM files.
#' @param qsm_path PARAM_DESCRIPTION, Default: '.'
#' @param recursive PARAM_DESCRIPTION, Default: F
#' @param treeqsm_path PARAM_DESCRIPTION, Default: ''
#' @return OUTPUT_DESCRIPTION
#' @details this may not work if your version of treeqsm doesn't match the version used when making the .mat files.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname save_model_text_matlab_server
#'
#' @import tools
save_model_text_matlab_server <- function(qsm_path = ".", recursive = F, treeqsm_path = "") {
    # this may not work if your version of treeqsm doesn't match the version used when making the .mat files
    if (!requireNamespace("R.matlab", quietly = TRUE)) {
        stop("Package \"R.matlab\" and Matlab > v6 needed for this function to work. Please install it.",
             call. = FALSE)
    }

    # start and open connection to matlab server
    R.matlab::Matlab$startServer()
    matlab <- R.matlab::Matlab()
    R.matlab::open(matlab)

        # add treeqsm funcs to matlab path
    if (treeqsm_path != "") {
        R.matlab::evaluate(matlab, paste("addpath(genpath('",treeqsm_path,"'))"))
    }

    qsm_mats = list.files(qsm_path, pattern = "\\.mat", full.names = T, recursive = recursive)
    pb <- txtProgressBar(max = length(qsm_mats), style = 3)
    setTxtProgressBar(pb, 0)
    for (this_qsm_mat in qsm_mats) {
        this_dir = dirname(this_qsm_mat)
        this_file_basename = basename(tools::file_path_sans_ext(this_qsm_mat))

        R.matlab::evaluate(matlab, paste0("save_model_text('", this_qsm_mat,"', '", this_file_basename,"');"))

        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    }
    close(pb)
    close(matlab)
    return(qsm_mats)
}


#' @title Get children of internode
#'
#' @param ts treestruct dataframe (not object)
#' @param internode_id if this is provided, the id of the internode for which children are sought
#' @param row if this is provided, the row number of the internode for which children are sought
#'
#' @return vector of rows of internode children, including the row of the specified internode
#' @export
#' @details only one of internode_id and row arguments are required.
#' @examples
get_child_rows <- function(ts, internode_id, row) {

    if (!missing(internode_id)) thisrow = which(internode_id == ts$internode_id)
    else if (!missing(row)) thisrow = row
    else stop("One of internode_id and row must be specified.")
    return(get_children(thisrow, ts$parent_row))
}


#' @title Get treestruct of children of internode
#'
#' @param ts treestruct dataframe (not object)
#' @param internode_id if this is provided, the id of the internode for which children are sought
#' @param row if this is provided, the row number of the internode for which children are sought
#'
#' @return treestruct comprised of of rows of internode children, including the row of the specified internode.  Parent_rows are re-computed.
#' @export
#' @details only one of internode_id and row arguments are required.
#' @examples
get_branch <- function(ts, internode_id, row) {

    if (!missing(internode_id)) thisrow = which(internode_id == ts$internode_id)
    else if (!missing(row)) thisrow = row
    else stop("One of internode_id and row must be specified.")

    child_ts = ts[get_child_rows(ts, row = thisrow),]
    child_ts$parent_row = parent_row(parent_id = child_ts$parent_id, internode_id = child_ts$internode_id)
    return(child_ts)
}

#' @title Prune a branch from a treestruct dataframe
#'
#' @param ts treestruct dataframe (not object)
#' @param internode_id if this is provided, the id of the internode for which children are sought
#' @param row if this is provided, the row number of the internode for which children are sought
#'
#' @return treestruct pruned of branch starting with the specified internode.  Parent_rows are re-computed.
#' @export
#' @details only one of internode_id and row arguments are required.
#' @examples
prune_branch <- function(ts, internode_id, row) {

    if (!missing(internode_id)) thisrow = which(internode_id == ts$internode_id)
    else if (!missing(row)) thisrow = row
    else stop("One of internode_id and row must be specified.")

    branch_rows = get_child_rows(ts, row = thisrow)
    ts = ts[-branch_rows,]
    ts$parent_row = parent_row(parent_id = ts$parent_id, internode_id = ts$internode_id)
    return(ts)
}
