% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeStructs_methods.r
\name{visNetwork.TreeStructs}
\alias{visNetwork.TreeStructs}
\title{visNetwork.TreeStructs}
\usage{
\method{visNetwork}{TreeStructs}(ts, index, hierarchical = T, width_factor = 100, length_factor = 10)
}
\arguments{
\item{index}{\code{character} or \code{integer} char matching TreeStructs id column, or index of branch in TreeStructs object (required)}

\item{hierarchical}{\code{logical}, plot in hierarchical layout. Default: T}

\item{width_factor}{\code{numeric} edge width factor, Default: 100}

\item{length_factor}{\code{numeric} edge length factor, Default: 10}

\item{bss}{TreeStructs object (requred)}
}
\value{
visNetwork object
}
\description{
visualize tree network.  This overwrites visNetwork for TreeStructs objects.
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
 visNetwork(bss_obj) \%>\% visEdges(arrows = "middle")
 }
}
}
\seealso{
\code{\link[visNetwork]{visNetwork-igraph}},\code{\link[visNetwork]{visNetwork}}
}
