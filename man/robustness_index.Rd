% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustness.index_function.R
\name{robustness_index}
\alias{robustness_index}
\title{Robustness index}
\usage{
robustness_index(sens, prediction = FALSE, threshold)
}
\arguments{
\item{sens}{An object of S3 class \code{\link{run_sensitivity}} when
sensitivity analysis refers to different scenarios about the average
missingness parameter. See 'Value' in \code{\link{run_sensitivity}}. For a
\bold{general} sensitivity analysis, insert a list of at least two objects
of S3 class \code{\link{run_model}} or \code{\link{run_metareg}} indicating
different re-analyses: the first object in the list should refer to
the primary analysis.}

\item{prediction}{Logical character on whether to consider the prediction
(\code{TRUE}) or estimation of the summary treatment effects
(\code{FALSE}). This is only relevant for a random-effects model and the
default argument is \code{FALSE} (estimation).}

\item{threshold}{A number indicating the threshold of robustness, that is,
the minimally allowed deviation between the primary analysis and
re-analysis results. See 'Details' below.}
}
\value{
\code{robustness_index} prints on the R console a message in green
  text on the threshold of robustness determined by the user.
  Then, the function returns the following list of elements:
  \item{robust_index}{A numeric scalar or vector on the robustness
  index values. In the case of a pairwise meta-analysis,
  \code{robust_index} is scalar as only one summary effect size is obtained.
  In the case of network meta-analysis, \code{robust_index} is a vector with
  length equal to the number of possible pairwise comparisons;
  one robustness index per possible comparison.}
  \item{robust}{A character or character vector (of same length with
  \code{robust_index}) on whether the primary analysis results are
  \emph{robust} or \emph{frail} to the different re-analyses.}
  \item{kld}{A vector or matrix on the Kullback-Leibler divergence
  (KLD) measure in the summary effect size from a subsequent re-analysis to
  the primary analysis. In the case of a pairwise meta-analysis, \code{kld}
  is a vector with length equal to the number of total analyses (one KLD
  value is obtained per analysis). The number of total analyses equals the
  square of the number of scenarios indicated in the argument
  \code{mean_scenarios} of \code{\link{run_sensitivity}}, in the case of
  missing participant outcome data; otherwise, the length of the character
  vector in argument \code{sens}.
  In the case of network meta-analysis, \code{robust_index} is a matrix with
  number of rows equal to the number of total analyses and number of columns
  equal to the number of  possible pairwise comparisons; one KLD value per
  analysis and possible comparison.}
  \item{threshold}{The threshold used to be inherited by the
  \code{\link{heatmap_robustness}} function.}
  \item{scenarios}{The scenarios considered to be inherited by the
  \code{\link{kld_barplot}} function.}
}
\description{
Calculates the robustness index, a novel index that quantifies the overall
  divergence of the sensitivity analysis results from the primary analysis
  results. The robustness index considers objective decision rules to infer
  the presence or lack of robustness of the primary analysis results when
  conducting a sensitivity analysis (Spineli et al., 2021).
}
\details{
Thresholds of robustness have been proposed only for the odds ratio
  and standardised mean difference (Spineli et al., 2021).
  The user may consider the values 0.28 and 0.17 in the argument
  \code{threshold} for the odds ratio and standardised mean difference effect
  measures (the default values), respectively, or consider other plausible
  values. When the argument \code{threshold} has not been defined,
  \code{robustness_index} considers the default values 0.28 and 0.17 as
  threshold for robustness for binary and continuous outcome, respectively,
  regardless of the effect measure (the default thresholds may not be proper
  choices for other effect measures; hence, use these threshold with great
  caution in this case). Spineli et al. (2021) offers a discussion on
  specifying the \code{threshold} of robustness.

  In the case of binary outcome, \code{robustness_index} considers the
  results in the odds ratio scale to calculate the robustness index.
  This is because, the odds ratio is used as the 'best-case' effect measure
  in \code{\link{run_model}}. Then, relative risk, and risk difference are
  functions of the odds ratio and the selected baseline risk (See 'Details'
  in \code{\link{run_model}}).

  In the case of missing participant outcome data, the primary analysis is
  considered to be the middle of the numbers in the argument
  \code{mean_scenarios} of \code{\link{run_sensitivity}} (see 'Arguments'
  and 'Details' in \code{\link{run_sensitivity}}). Furhermore,
  \code{robustness_index} can be used in that context only when missing
  participant outcome data have been extracted for at least one trial.
  Otherwise, the execution of the function will be stopped and an error
  message will be printed in the R console.

  In the case of a general sensitivity analysis, the compared models should
  refer to the same effect measure and the same meta-analysis model (i.e.,
  fixed-effect or random-effects).

  In \code{robust}, the value \code{"robust"} appears when
  \code{robust_index} is less than \code{threshold}; otherwise, the value
  \code{"frail"} appears.
}
\examples{
data("nma.baker2009")

# Read results from 'run_sensitivity' (using the default arguments)
res_sens <- readRDS(system.file('extdata/res_sens_baker.rds',
                    package = 'rnmamod'))

# Calculate the robustness index
robustness_index(sens = res_sens,
                 threshold = 0.28)

}
\references{
Kullback S, Leibler RA. On information and sufficiency.
\emph{Ann Math Stat} 1951;\bold{22}(1):79--86. doi: 10.1214/aoms/1177729694

Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
primary analysis results: A case study on missing outcome data in pairwise
and network meta-analysis.
\emph{Res Synth Methods} 2021;\bold{12}(4):475--90. doi: 10.1002/jrsm.1478
}
\seealso{
\code{\link{heatmap_robustness}}, \code{\link{kld_barplot}},
  \code{\link{kld_measure}}, \code{\link{run_metareg}},
  \code{\link{run_model}}, \code{\link{run_sensitivity}}
}
\author{
{Loukia M. Spineli}
}
