% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scatterplot.sucra_function.R
\name{scatterplot_sucra}
\alias{scatterplot_sucra}
\title{Scatterplot of SUCRA values}
\usage{
scatterplot_sucra(full, reg, cov_name = "covariate value", drug_names)
}
\arguments{
\item{full}{An object of S3 class \code{\link{run_model}}.
See 'Value' in \code{\link{run_model}}.}

\item{reg}{An object of S3 class \code{\link{run_metareg}}. See 'Value' in
\code{\link{run_metareg}}.}

\item{cov_name}{A character or text to indicate the name of the covariate.}

\item{drug_names}{A vector of labels with the name of the interventions in
the order they appear in the argument \code{data} of
\code{\link{run_model}}. If \code{drug_names} is not defined,
the order of the interventions as they appear in \code{data} is used,
instead.}
}
\value{
A scatterplot of the SUCRA values under the network meta-analysis
  (y-axis) against the SUCRA values under the network meta-regression
  (x-axis) for a specified level or value of the investigated covariate.
}
\description{
Creates a scatterplot of the SUCRA values from the
  network meta-analysis and the network meta-regression for a specified level
  or value of the investigated covariate.
}
\details{
The names of the interventions appear above each point in the plot.
  Three coloured rectangles are drawn in the scatterplot: a red rectangle for
  SUCRA values up to 50\%, a yellow rectangular for SUCRA values between
  50\% and 80\%, and a green rectangle for SUCRA values over 80\%.
  Interventions falling at the green area are considered as the highest
  ranked interventions, whilst interventions falling at the red area are
  considered as the lowest ranked interventions.

  \code{scatterplot_sucra} is integrated in \code{\link{metareg_plot}}.

  \code{scatterplot_sucra} can be used only for a network of interventions.
  Otherwise, the execution of the function will be stopped and an error
  message will be printed on the R console.
}
\references{
Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
for presenting results from multiple-treatment meta-analysis: an overview and
tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
doi: 10.1016/j.jclinepi.2010.03.016
}
\seealso{
\code{\link{metareg_plot}}, \code{\link{run_metareg}},
  \code{\link{run_model}}
}
\author{
{Loukia M. Spineli}
}
