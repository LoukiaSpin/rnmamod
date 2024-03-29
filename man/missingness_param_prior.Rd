% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missingness.param.prior_function.R
\name{missingness_param_prior}
\alias{missingness_param_prior}
\title{Define the mean value of the normal distribution of the missingness parameter}
\usage{
missingness_param_prior(assumption, mean_misspar)
}
\arguments{
\item{assumption}{Character string indicating the structure of the
informative missingness parameter.
Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"},
\code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"},
\code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or
\code{"IND-UNCORR"}.
The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"},
\code{"HIE"}, and \code{"IND"} stand for identical, hierarchical and
independent, respectively. \code{"CORR"} and \code{"UNCORR"} stand for
correlated and uncorrelated, respectively.}

\item{mean_misspar}{A numeric value or a vector of two numeric values for the
mean of the normal distribution of the informative missingness parameter
(see 'Details'). The default argument is 0 and corresponds to the
missing-at-random assumption for \code{assumption = "IDE-ARM"}.}
}
\value{
A scalar or numeric vector to be passed to \code{\link{run_model}}.
}
\description{
Generates the mean value of the normal distribution of the missingness
  parameter in the proper format depending on the assumed structure of the
  missingness parameter. \code{\link{run_model}} inherits
  \code{missingness_param_prior} through the argument \code{mean_misspar}
  (see 'Argument' in \code{\link{run_model}}).
}
\details{
\code{\link{run_model}} considers the informative missingness odds
  ratio in the logarithmic scale for binary outcome data (Spineli, 2019a;
  Turner et al., 2015; White et al., 2008), the informative missingness
  difference of means when \code{measure} is \code{"MD"} or \code{"SMD"},
  and the informative missingness ratio of means in the logarithmic scale
  when \code{measure} is \code{"ROM"} (Spineli et al., 2021;
  Mavridis et al., 2015).

  When \code{assumption} is trial-specific (i.e., \code{"IDE-TRIAL"} or
  \code{"HIE-TRIAL"}), or independent (i.e., \code{"IND-CORR"} or
  \code{"IND-UNCORR"}), only one numeric value can be assigned to
  \code{mean_misspar} because the same missingness scenario is applied to all
  trials and trial-arms of the dataset, respectively. When \code{assumption}
  is \code{"IDE-ARM"} or \code{"HIE-ARM"}, a maximum of two
  \emph{different or identical} numeric values can be assigned as a vector to
  \code{mean_misspars}: the first value refers to the experimental arm,
  and the second value refers to the control arm of a trial.
  In the case of a network, the first value is considered for all
  non-reference interventions and the second value is considered for the
  reference intervention of the network (see 'Argument' \code{ref} in
  \code{\link{run_model}}). This is necessary to ensure transitivity in the
  assumptions for the missingness parameter across the comparisons in the
  network (Spineli, 2019b).

  Currently, there are no empirically-based prior distributions for the
  informative missingness parameters. The users may refer to
  Mavridis et al. (2015) and Spineli (2019) to determine \code{mean_misspar}
  for an informative missingness parameter.
}
\references{
Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
uncertainty due to missing continuous outcome data in pairwise and
network meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--41.
doi: 10.1002/sim.6365

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
data in network meta-analysis: a one-stage pattern-mixture model approach.
\emph{Stat Methods Med Res} 2021;\bold{30}(4):958--75.
doi: 10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis.
\emph{BMC Med Res Methodol} 2019a;\bold{19}(1):86.
doi: 10.1186/s12874-019-0731-y

Spineli LM. Modeling missing binary outcome data while preserving
transitivity assumption yielded more credible network meta-analysis
results. \emph{J Clin Epidemiol} 2019b;\bold{105}:19--26.
doi: 10.1016/j.jclinepi.2018.09.002

Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
for uncertainty due to missing binary outcome data in pairwise meta-analysis.
\emph{Stat Med} 2015;\bold{34}(12):2062--80. doi: 10.1002/sim.6475

White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
in meta-analysis--part 1: two-stage methods. \emph{Stat Med}
2008;\bold{27}(5):711--27. doi: 10.1002/sim.3008
}
\seealso{
\code{\link{run_model}}
}
\author{
{Loukia M. Spineli}
}
