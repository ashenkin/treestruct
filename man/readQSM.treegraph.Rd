% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QSM_funcs.r
\name{readQSM.treegraph}
\alias{readQSM.treegraph}
\title{readQSM.treegraph}
\usage{
readQSM.treegraph(qsmfile, qsmver = "by_name")
}
\arguments{
\item{qsmfile}{path to treegraph file}

\item{qsmver}{version of treegraph used to make the treegraph file, Default: by_name}
}
\value{
a 3-element list of file (filename), CylData, and BranchData
}
\description{
reads in a QSM file in treegraph format and passes back cylinder and branch data
}
\details{
Treegraphs are described by nodes that connect cylinders, as opposed to internodes as in treeQSM.
In addition, treegraphs sometimes have furcations at the most basal node (think of a tree that branches at the soil
surface).  In these cases, treestruct's use of internodes instead of nodes as the most basic element breaks down,
as there is no way to indicate that the two (or more) most basal cylinders should be connected.  To get around this issue,
we add a 0-length, 0-width cylinder at the base of treegraph objects with basal furcations to enable that connectivity.  We must make
sure that analyses that could be thrown off by this extra cylinder make accomodations for it.
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
