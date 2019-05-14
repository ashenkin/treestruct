# functions for interacting with matlab QSM objects

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
readQSM.mat <- function(qsmfile, qsmver = 2.3) {
    if (!requireNamespace("rmatio", quietly = TRUE)) {
        stop("Package \"rmatio\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    QSMmat = rmatio::read.mat(qsmfile)
    # QSM files have 4 elements: treedata, inputs, qsm, and models

    #names(QSMmat[["qsm"]]) = c("cylinder", "branch", "treedata", "rundata", "pmdistance", "triangulation")  # names below 1st level in hierarchy dropped by R.matlab it seems

    cyldata = QSMmat$qsm$cylinder
    cyldata = cyldata %>% purrr::map(`[[`, 1) %>% as.data.frame() %>% dplyr::bind_cols()
    branchdata = QSMmat$qsm$branch
    branchdata = branchdata %>% purrr::map(`[[`, 1) %>% as.data.frame() %>% dplyr::bind_cols()

    # CylData_cols = c("rad", "len", "sta", "axe", "cpar", "cext", "boc", "added", "unmodradius")
    # BranchData_cols = c("bord", "bpar", "bvol", "blen", "bang")
    #
    # cyldata_ele2col_order = order(match(tolower(names(QSMmat)), CylData_cols), na.last = NA)
    # branchdata_ele2col_order = order(match(tolower(names(QSMmat)), BranchData_cols), na.last = NA)

    #cydata_elements = grepl(paste0("^", paste(CylData_cols, collapse="$|^"), "$"), names(QSMmat), ignore.case = T)
    #branchdata_elements = grepl(paste0("^", paste(BranchData_cols, collapse="$|^"), "$"), names(QSMmat), ignore.case = T)

    # CylData = data.frame(do.call(cbind, QSMmat[ cyldata_ele2col_order ]))
    # BranchData = data.frame(do.call(cbind, QSMmat[ branchdata_ele2col_order ]))

    #names(CylData) =

    # reorder and rename columns


    if (is.numeric(qsmver) & qsmver >= 2.3) {
        names(cyldata) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row","branch_data_row","branch_order","index_num","added_after", "rad0")
        names(branchdata) = c("bord", "bpar", "bvol", "blen", "bang")
    } else if (qsmver == "UCL") {
        names(cyldata) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row", "added_after", "unmod_rad", "branch_data_row","branch_order","index_num")
        names(branchdata) = c("bord", "bpar", "bvol", "blen", "bang", "bheight", "baz", "bdiam")
    } else {
        names(cyldata) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row","branch_data_row","branch_order","index_num","added_after")
        names(branchdata) = c("bord", "bpar", "bvol", "blen", "bang")
    }

    return(list(file = basename(qsmfile), CylData = cyldata, BranchData = branchdata))
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
#' @param qsmver PARAM_DESCRIPTION, Default: 2.3
#' @param filematch_pattern regex filematch pattern (used in list.files), Default: ".mat"
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
import_treestructs_from_dir <- function(qsm_path = ".", recursive = F, qsmver = 2.3, filematch_pattern = ".mat", nest = F) {
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
        qsms[[qsmfile]] = readQSM.mat(qsmfile, qsmver = qsmver)
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

