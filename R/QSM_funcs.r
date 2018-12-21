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
readQSM.mat <- function(qsmfile, qsmver = 2.3) {
    if (!requireNamespace("R.matlab", quietly = TRUE)) {
        stop("Package \"R.matlab\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    QSMmat = R.matlab::readMat(qsmfile)

    CylData_cols = c("rad", "len", "sta", "axe", "cpar", "cext", "boc", "added", "unmodradius")
    BranchData_cols = c("bord", "bpar", "bvol", "blen", "bang")

    cyldata_ele2col_order = order(match(tolower(names(QSMmat)), CylData_cols), na.last = NA)
    branchdata_ele2col_order = order(match(tolower(names(QSMmat)), BranchData_cols), na.last = NA)

    #cydata_elements = grepl(paste0("^", paste(CylData_cols, collapse="$|^"), "$"), names(QSMmat), ignore.case = T)
    #branchdata_elements = grepl(paste0("^", paste(BranchData_cols, collapse="$|^"), "$"), names(QSMmat), ignore.case = T)

    CylData = data.frame(do.call(cbind, QSMmat[ cyldata_ele2col_order ]))
    BranchData = data.frame(do.call(cbind, QSMmat[ branchdata_ele2col_order ]))

    #names(CylData) =

    # reorder and rename columns
    names(BranchData) = c("bord", "bpar", "bvol", "blen", "bang")

    if (qsmver >= 2.3) {
        names(CylData) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row","branch_data_row","branch_order","index_num","added_after", "rad0")
    } else {
        names(CylData) = c("rad","len","x_start","y_start","z_start","x_cyl","y_cyl","z_cyl","parent_row",
                           "daughter_row","branch_data_row","branch_order","index_num","added_after")
    }

    return(list(file = basename(qsmfile), CylData = CylData, BranchData = BranchData))
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

