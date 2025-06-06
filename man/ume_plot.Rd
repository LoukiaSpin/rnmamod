% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UME.plot_function.R
\name{ume_plot}
\alias{ume_plot}
\title{End-user-ready results for the unrelated mean effects model}
\usage{
ume_plot(full, ume, drug_names, save_xls)
}
\arguments{
\item{full}{An object of S3 class \code{\link{run_model}}. See 'Value' in
\code{\link{run_model}}.}

\item{ume}{An object of S3 class \code{\link{run_ume}}. See 'Value' in
\code{\link{run_ume}}.}

\item{drug_names}{A vector of labels with the name of the interventions in
the order they appear in the argument \code{data} of
\code{\link{run_model}}. If \code{drug_names} is not defined,
the order of the interventions as they appear in \code{data} is used,
instead.}

\item{save_xls}{Logical to indicate whether to export the tabulated results
to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
function of the R-package
\href{https://CRAN.R-project.org/package=writexl}{writexl}) to the working
directory of the user. The default is \code{FALSE} (do not export).}
}
\value{
\code{ume_plot} prints on the R console a message on the most
  parsimonious model (if any) based on the DIC (red text). Then, the function
  returns the following list of elements:
  \item{table_effect_size}{The posterior median, posterior standard
  deviation, and 95\% credible interval of the summary effect size for each
  pairwise comparison observed in the network under the consistency model
  and the unrelated mean effects model.}
  \item{table_model_assessment}{The DIC, number of effective
  parameters, and total residual deviance under the consistency model and
  the unrelated mean effects model (Spiegelhalter et al., 2002).}
  \item{table_tau}{The posterior median and 95\% credible interval of
  \emph{tau} under the consistency model and the unrelated mean effects
  model. When a fixed-effect model has been performed, \code{ume_plot} does
  not return this element.}
  \item{scatterplots}{The scatterplot and the Bland-Altman plot on the
  posterior mean deviance contribution of the individual data points under
  the consistency model and the unrelated mean effects model. See 'Details'
  and 'Value' in \code{\link{scatterplots_dev}} and
  \code{\link{bland_altman_plot}}, respectively.}
  \item{levarage_plots}{The leverage plot under the consistency model
  and the unrelated mean effects model, separately. See 'Details' and
  'Value' in \code{\link{leverage_plot}}.}
  \item{intervalplots}{A panel of interval plots on the summary effect
  size under the consistency model and the unrelated mean effects model for
  each pairwise comparison observed in the network. See 'Details' and
  'Value' in \code{\link{intervalplot_panel_ume}}.}
}
\description{
\code{ume_plot} hosts a toolkit of functions that facilitates
  the comparison of the consistency model (via \code{\link{run_model}}) with
  the unrelated mean effects model (via \code{\link{run_ume}}) regarding the
  posterior summaries of the summary effect size for the pairwise comparisons
  observed in the network, the between-trial standard deviation (\emph{tau})
  and model assessment parameters.
}
\details{
The deviance information criterion (DIC) of the consistency model is
  compared with the DIC of the unrelated mean effects model
  (Dias et al., 2013). If the difference in DIC exceeds 5, the unrelated mean
  effects model is preferred. If the difference in DIC is less than -5, the
  consistency is preferred; otherwise, there is little to choose between the
  compared models.

  For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
  (risk difference) in \code{\link{run_model}}, \code{ume_plot} currently
  presents the results from network meta-analysis and unrelated mean effects
  in the odds ratio for being the \strong{base-case} effect measure in
  \code{\link{run_model}} for a binary outcome (see also 'Details' in
  \code{\link{run_model}}).

  Furthermore, \code{ume_plot} exports \code{table_effect_size} and
  \code{table_model_assessment} to separate 'xlsx' files (via the
  \code{\link[writexl:write_xlsx]{write_xlsx}} function) to the working
  directory of the user.

  \code{ume_plot} can be used only for a network of interventions. In the
  case of two interventions, the execution of the function will be stopped
  and an error message will be printed on the R console.
}
\examples{
data("nma.liu2013")

\donttest{

# Read results from 'run_model' (using the default arguments)
res <- readRDS(system.file('extdata/res_liu.rds', package = 'rnmamod'))

# Read results from 'run_ume' (using the default arguments)
ume <- readRDS(system.file('extdata/ume_liu.rds', package = 'rnmamod'))

# The names of the interventions in the order they appear in the dataset
interv_names <- c("placebo", "pramipexole",
                  "serotonin norepinephrine reuptake inhibitor",
                  "serotonin reuptake inhibitor",
                  "tricyclic antidepressant", "pergolide")

# Plot the results from both models
ume_plot(full = res,
         ume = ume,
         drug_names = interv_names)
}

}
\references{
Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
for decision making 4: inconsistency in networks of evidence based on
randomized controlled trials.
\emph{Med Decis Making} 2013;\bold{33}(5):641--56.
doi: 10.1177/0272989X12455847

Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of
model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}(4):583--396.
doi: 10.1111/1467-9868.00353
}
\seealso{
\code{\link{bland_altman_plot}},
  \code{\link{intervalplot_panel_ume}}, \code{\link{leverage_plot}},
  \code{\link{run_model}}, \code{\link{run_ume}},
  \code{\link[writexl:write_xlsx]{write_xlsx}}
}
\author{
{Loukia M. Spineli}
}
