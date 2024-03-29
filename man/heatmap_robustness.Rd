% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.robustness_function.R
\name{heatmap_robustness}
\alias{heatmap_robustness}
\title{Heatmap of robustness}
\usage{
heatmap_robustness(robust, drug_names)
}
\arguments{
\item{robust}{An object of S3 class \code{\link{robustness_index}} and
\code{\link{robustness_index_user}}.
See 'Value' in \code{\link{robustness_index}} and
\code{\link{robustness_index_user}}.}

\item{drug_names}{A vector of labels with the name of the interventions in
the order they appear in the argument \code{data} of
\code{\link{run_model}}. If \code{drug_names} is not defined,
the order of the interventions as they appear in \code{data} is used,
instead.}
}
\value{
\code{heatmap_robustness} first prints on the R console a message on
  the threshold of robustness determined by the user in
  \code{\link{robustness_index}} and \code{\link{robustness_index_user}}.
  Then, it returns a lower triangular heatmap matrix with the robustness
  index value of all possible pairwise comparisons.
}
\description{
Facilitates the detection of comparisons that are associated
  with a lack of robustness in the context of a sensitivity analysis.
}
\details{
The heatmap illustrates the robustness index for each possible
  pairwise comparison in the network. The pairwise comparisons are read from
  left to right. Comparisons highlighted with green or red colour imply
  robust or frail conclusions for the primary analysis, respectively.
  This corresponds to robustness index below or at least the selected
  threshold of robustness. \code{heatmap_robustness} inherits the threshold
  of robustness selected in the \code{\link{robustness_index}} or
  \code{\link{robustness_index_user}} function.
  The robustness index of each pairwise comparison also appears in the
  corresponding cell.
  When there is at least one comparison with frail conclusions, the primary
  analysis results may be questionable for the whole network
  (Spineli et al., 2021).

  \code{heatmap_robustness} is \emph{not} restricted to the sensitivity
  analysis concerning the impact of missing participant outcome data.

  \code{heatmap_robustness} can be used only for a network of interventions.
  Otherwise, the execution of the function will be stopped and an
  error message will be printed on the R console.
}
\examples{
data("nma.baker2009")

# Read results from 'run_sensitivity' (using the default arguments)
res_sens <- readRDS(system.file('extdata/res_sens_baker.rds',
                    package = 'rnmamod'))

# Calculate the robustness index
robust <- robustness_index(sens = res_sens,
                           threshold = 0.28)

# The names of the interventions in the order they appear in the dataset
interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
                  "fluticasone", "fluticasone plus salmeterol",
                  "formoterol", "salmeterol", "tiotropium")

# Create the heatmap of robustness
heatmap_robustness(robust = robust,
                   drug_names = interv_names)

}
\references{
Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
primary analysis results: A case study on missing outcome data in pairwise
and network meta-analysis.
\emph{Res Synth Methods} 2021;\bold{12}(4):475--90. doi: 10.1002/jrsm.1478
}
\seealso{
\code{\link{robustness_index}}, \code{\link{robustness_index_user}},
  \code{\link{run_model}}
}
\author{
{Loukia M. Spineli}
}
