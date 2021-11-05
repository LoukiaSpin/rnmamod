#' rnmamod: Bayesian Network Meta-analysis with Missing Outcome Data in R
#'
#' @description
#'   An R package for performing Bayesian network meta-analysis while handling
#'   missing participant outcome data properly.
#'
#' @details
#'   R-package \bold{rnmamod} is built upon the WinBUGS program code found in
#'   the series of tutorial papers on evidence synthesis methods for decision
#'   making (Dias et al., 2013a; Dias et al., 2013b; Dias et al., 2013c) and
#'   Dias et al., (2010) that introduces the node-spliting approach.
#'   All models comprise Bayesian hierarchical models for one-stage network
#'   meta-analysis and they are implemented in JAGS through the R-package
#'   \bold{R2jags}.
#'
#'   \bold{rnmamod} comprises a suite of core models implemented in a
#'   systematic review with multiple interventions: fixed-effect and
#'   random-effects network meta-analysis, network meta-regression models for
#'   trial characteristics, and assessment of the consistency assumption
#'   globally and locally.
#'   \bold{rnmamod} also includes a rich suite of visualisation tools to aid in
#'   interpretation of the results and preparation of the manuscript for
#'   submission.
#'
#'   Missing participant outcome data are addressed in all models of the package
#'   after extending the code to incorporating the pattern-mixture model
#'   (Spineli, 2019; Spineli et al., (2021)).
#'
#'   The source for \bold{rnmamod} is available on
#'   [GitHub](https://github.com/LoukiaSpin/rnmamod) under the GPL-3.0 License.
#'
#'
#' @author {Loukia M. Spineli}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
#' treatment comparison meta-analysis.
#' \emph{Stat Med} 2010;\bold{29}(7-8):932--44.
#' \doi{10.1002/sim.3767}
#'
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013a;\bold{33}(5):641--56.
#' \doi{10.1177/0272989X12455847}
#'
#' Dias S, Sutton AJ, Welton NJ, Ades AE. Evidence synthesis for decision making
#' 3: heterogeneity--subgroups, meta-regression, bias, and bias-adjustment.
#' \emph{Med Decis Making} 2013b;\bold{33}(5):618--40.
#' \doi{10.1177/0272989X13485157}
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013c;\bold{33}(5):607--617. \doi{10.1177/0272989X12458724}
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' \doi{10.1186/s12874-019-0731-y}
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021. \doi{10.1177/0962280220983544}
#'
#' @export
