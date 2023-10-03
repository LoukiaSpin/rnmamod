% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.measures.plot_function.R
\name{internal_measures_plot}
\alias{internal_measures_plot}
\title{Internal measures for cluster validation
(Comparisons' comparability for transitivity evaluation)}
\usage{
internal_measures_plot(input, num_neighb, optimal_link)
}
\arguments{
\item{input}{An object of 'dist' class. It is a lower off-diagonal matrix
with the dissimilarities of all pairs of comparisons.}

\item{num_neighb}{A positive integer for the number of neighbouring
comparisons. It takes values from two to the number of comparisons minus
one. The default argument equals half the number of comparisons.}

\item{optimal_link}{A character string with values \code{"ward.D"},
\code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
\code{"mcquitty"}, \code{"median"}, or \code{"centroid"} for the optimal
linkage method, corresponding to the highest cophenetic correlation
coefficient value.}
}
\value{
\code{internal_measures_plot} returns the following list of elements:
  \item{Table_internal_measures}{A data-frame of the connectivity index,
  silhouette width, and Dunn index for a range of 2 to P-1 clusters, with P
  being the number of comparisons.}
  \item{Internal_measures_panel}{A panel of profile plots on the connectivity
  index, silhouette width, and Dunn index for a range of 2 to P-1 clusters,
  with P being the number of comparisons. The candidate optimal number of
  clusters is indicated with a red point directly on the line.}
}
\description{
\code{internal_measures_plot} prepares the table with the results of three
  internal measures (connectivity index, silhouette width, and Dunn index)
  for a range of clusters, and visualises the results using a profile plot
  for each internal measure.
}
\details{
\code{internal_measures_plot} call the functions
  \code{\link{connectivity_index}}, \code{\link{silhouette_index}} and
  \code{\link{dunn_index}} to calculate the corresponding internal measures.
  \code{internal_measures_plot} also call the function
  \code{\link{comp_clustering}} to define the argument \code{optimal_link} to
  create the silhouette plot for the selected number of clusters.

  If the network has three observed comparisons,
  \code{internal_measures_plot} will return only the
  \code{Table_internal_measures}. This is because with only three observed
  comparisons, only two clusters can be considered by the internal measures.

  \code{internal_measures_plot} is integrated in the function
  \code{\link{comp_clustering}}.
}
\references{
Handl J, Knowles J, Kell DB. Computational cluster validation in post-genomic
data analysis. \emph{Biometrics} 2005;\bold{21}(15):3201--120.
doi: 10.1093/bioinformatics/bti517

Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
validation of cluster analysis.
\emph{J Comput Appl Math} 1987;\bold{20}:53--65.

Dunn J. Well-separated clusters and optimal fuzzy partitions.
\emph{J Cybern} 1974;\bold{4}(1):95--104.
}
\seealso{
\code{\link{comp_clustering}}, \code{\link{connectivity_index}},
 \code{\link{dunn_index}}, \code{\link{silhouette_index}}
}
\author{
{Loukia M. Spineli}
}