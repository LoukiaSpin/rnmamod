#' Perform Bayesian pairwise or network meta-analysis
#'
#' @description
#'   Performs a one-stage pairwise or network meta-analysis while addressing
#'   aggregate binary or continuous missing participant outcome data via the
#'   pattern-mixture model.
#'
#' @param data A data-frame of the one-trial-per-row format with arm-level data.
#'   See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure. For a binary
#'   outcome, the following can be considered: \code{"OR"}, \code{"RR"} or
#'   \code{"RD"} for the odds ratio, relative risk, and risk difference,
#'   respectively. For a continuous outcome, the following can be considered:
#'   \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for mean difference,
#'   standardised mean difference and ratio of means, respectively.
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. The default argument is \code{"RE"}.
#' @param assumption Character string indicating the structure of the
#'   informative missingness parameter.
#'   Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"},
#'   \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"},
#'   \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"},
#'   or \code{"IND-UNCORR"}.
#'   The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"},
#'   \code{"HIE"}, and \code{"IND"} stand for identical, hierarchical and
#'   independent, respectively. \code{"CORR"} and \code{"UNCORR"} stand for
#'   correlated and uncorrelated, respectively.
#' @param heter_prior A list of three elements with the following order:
#'   1) a character string indicating the distribution with
#'   (currently available) values \code{"halfnormal"}, \code{"uniform"},
#'   \code{"lognormal"}, or \code{"logt"}; 2) two numeric values that refer to
#'   the parameters of the selected distribution. For \code{"lognormal"}, and
#'   \code{"logt"} these numbers refer to the mean and precision, respectively.
#'   For \code{"halfnormal"}, these numbers refer to zero and the scale
#'   parameter (equal to 4 or 1 being the corresponding precision of the scale
#'   parameter 0.5 or 1). For \code{"uniform"}, these numbers refer to the
#'   minimum and maximum value of the distribution.
#'   See 'Details' in \code{\link{heterogeneity_param_prior}}.
#' @param mean_misspar A scalar or numeric vector of two numeric values for the
#'   mean of the normal distribution of the informative missingness parameter
#'   (see 'Details'). The default argument is 0 and corresponds to the
#'   missing-at-random assumption.
#'   See also 'Details' in \code{\link{missingness_param_prior}}.
#' @param var_misspar A positive non-zero number for the variance of the
#'   normal distribution of the informative missingness parameter.
#'   When the \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"}
#'   the default argument is 1. When the \code{measure} is \code{"ROM"}
#'   the default argument is 0.04.
#' @param D A binary number for the direction of the outcome.
#'   Set \code{D = 1} for beneficial outcome and \code{D = 0} for harmful
#'   outcome.
#' @param ref An integer specifying the reference intervention. The number
#'   should match the intervention identifier under element \strong{t} in
#'   \code{data} (See 'Format').
#' @param base_risk A scalar, a vector of length three with elements sorted in
#'   ascending order, or a matrix with two columns and number of rows equal to
#'   the number of relevant trials. In the case of a scalar or vector, the
#'   elements should be in the interval (0, 1) (see 'Details'). If
#'   \code{base_risk} has not been defined, the function uses the median event
#'   risk for the reference intervention from the corresponding trials in
#'   \code{data}. This argument is only relevant for a binary outcome.
#' @param n_chains Positive integer specifying the number of chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Positive integer specifying the number of Markov chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Positive integer specifying the number of iterations to
#'   discard at the beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Positive integer specifying the thinning rate for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @format The columns of the data-frame in the argument \code{data} refer
#'   to the following elements for a continuous outcome:
#'   \tabular{ll}{
#'    \strong{t} \tab An intervention identifier in each arm.\cr
#'    \tab \cr
#'    \strong{y} \tab The observed mean value of the outcome in each arm.\cr
#'    \tab \cr
#'    \strong{sd} \tab The observed standard deviation of the outcome in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{m} \tab The number of missing participant outcome data in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{n} \tab The number of randomised participants in each arm.\cr
#'   }
#'
#'   For a binary outcome, the columns of the data-frame in the argument
#'   \code{data} refer to the following elements:
#'   \tabular{ll}{
#'    \strong{t} \tab An intervention identifier in each arm.\cr
#'    \tab \cr
#'    \strong{r} \tab The observed number of events of the outcome in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{m} \tab The number of missing participant outcome data in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{n} \tab The number of randomised participants in each arm.\cr
#'   }
#'   The number of rows in \code{data} equals the number of collected trials.
#'   Each element appears in \code{data} as many times as the maximum number of
#'   interventions compared in a trial of the dataset.
#'   In pairwise meta-analysis, the maximum number of arms is inherently two.
#'   The same holds for a network meta-analysis without multi-arm trials.
#'   In the case of network meta-analysis with multi-arm trials, the maximum
#'   number of arms exceeds two. See 'Examples' that illustrates the structure
#'   of \code{data} for a network with a maximum number of four arms.
#'   It is not a prerequisite of \code{run_model} that the multi-arm trials
#'   appear at the bottom of the dataset.
#'
#' @return A list of R2jags output on the summaries of the posterior
#'   distribution, and the Gelman-Rubin convergence diagnostic
#'   (Gelman et al., 1992) of the following monitored parameters for a
#'   fixed-effect pairwise meta-analysis:
#'   \item{EM}{The estimated summary effect measure (according to the argument
#'   \code{measure}).}
#'   \item{EM_LOR}{The estimated summary odd ratio in the logarithmic scale when
#'   \code{measure = "RR"} or \code{measure = "RD"}.}
#'   \item{dev_o}{The deviance contribution of each trial-arm based on the
#'   observed outcome.}
#'   \item{hat_par}{The fitted outcome at each trial-arm.}
#'   \item{phi}{The informative missingness parameter.}
#'
#'   For a fixed-effect network meta-analysis, the output additionally includes:
#'   \item{SUCRA}{The surface under the cumulative ranking curve for each
#'   intervention.}
#'   \item{SUCRA_LOR}{The surface under the cumulative ranking curve for each
#'   intervention under the odds ratio effect measure when \code{measure = "RR"}
#'   or \code{measure = "RD"}.}
#'   \item{effectiveneness}{The ranking probability of each intervention for
#'   every rank.}
#'
#'   For a random-effects pairwise meta-analysis, the output additionally
#'   includes the following elements:
#'   \item{EM_pred}{The predicted summary effect measure (according to the
#'   argument \code{measure}).}
#'   \item{EM_pred_LOR}{The predicted summary odds ratio in the logarithmic
#'   scale when \code{measure = "RR"} or \code{measure = "RD"}.}
#'   \item{delta}{The estimated trial-specific effect measure (according to the
#'   argument \code{measure}).}
#'   \item{tau}{The between-trial standard deviation.}
#'
#'   In network meta-analysis, \code{EM} and \code{EM_pred} refer to all
#'   possible pairwise comparisons of interventions in the network. Furthermore,
#'   \code{tau} is typically assumed to be common for all observed comparisons
#'   in the network. For a multi-arm trial, we estimate a total of \emph{T-1}
#'   \code{delta} for comparisons with the baseline intervention of the trial
#'   (found in the first column of the element \bold{t}), with \emph{T} being
#'   the number of interventions in the trial.
#'
#'   Furthermore, the output includes the following elements:
#'   \item{leverage_o}{The leverage for the observed outcome at each trial-arm.}
#'   \item{sign_dev_o}{The sign of the difference between observed and fitted
#'   outcome at each trial-arm.}
#'   \item{model_assessment}{A data-frame on the measures of model assessment:
#'   deviance information criterion, number of effective parameters, and total
#'   residual deviance.}
#'   \item{indic}{The sign of basic parameters in relation to the reference
#'   intervention as specified in argument \code{reg}}
#'   \item{jagsfit}{An object of S3 class \code{\link[R2jags:jags]{jags}} with
#'   the posterior results on all monitored parameters to be used in the
#'   \code{\link{mcmc_diagnostics}} function.}
#'
#'   The \code{run_model} function also returns the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{heter_prior},
#'   \code{mean_misspar}, \code{var_misspar}, \code{D}, \code{ref},
#'   \code{base_risk}, \code{n_chains}, \code{n_iter}, \code{n_burnin},
#'   and \code{n_thin} as specified by the user to be inherited by other
#'   functions of the package.
#'
#' @details The model runs in \code{JAGS} and the progress of the simulation
#'   appears on the R console. The output of \code{run_model} is used as an S3
#'   object by other functions of the package to be processed further and
#'   provide an end-user-ready output.
#'
#'   The \code{\link{data_preparation}} function is called to prepare the data
#'   for the Bayesian analysis. \code{\link{data_preparation}} creates the
#'   pseudo-data-frames \code{m_new}, and \code{I}, that have the same
#'   dimensions with the element \code{N}. \code{m_new} takes the zero
#'   value for the observed trial-arms with unreported missing participant
#'   outcome data (i.e., \code{m} equals \code{NA} for the corresponding
#'   trial-arms), the same value with \code{m} for the observed trial-arms with
#'   reported missing participant outcome data, and \code{NA} for the unobserved
#'   trial-arms. \code{I} is a dummy pseudo-data-frame and takes the value one
#'   for the observed trial-arms with reported missing participant outcome data,
#'   the zero value for the observed trial-arms with unreported missing
#'   participant outcome data (i.e., \code{m_new} equals zero for the
#'   corresponding trial-arms), and \code{NA} for the unobserved trial-arms.
#'   Thus, \code{I} indicates whether missing participant outcome data have been
#'   collected for the observed trial-arms. If the user has not defined the
#'   element \strong{m} in \code{data}, \code{m_new} and \code{I} take the zero
#'   value for all observed trial-arms to indicate that no missing participant
#'   outcome data have been collected for the analysed outcome. See 'Details' in
#'   \code{\link{data_preparation}}.
#'
#'   Furthermore, \code{\link{data_preparation}} sorts the interventions across
#'   the arms of each trial in an ascending order and correspondingly the
#'   remaining elements in \code{data} (see 'Format').
#'   \code{\link{data_preparation}} considers the first column in \strong{t} as
#'   being the control arm for every trial. Thus, this sorting ensures that
#'   interventions with a lower identifier are consistently treated as the
#'   control arm in each trial. This case is relevant in non-star-shaped
#'   networks.
#'
#'   To perform a Bayesian pairwise or network meta-analysis, the
#'   \code{\link{prepare_model}} function is called which contains the WinBUGS
#'   code as written by Dias et al. (2013a) for binomial and normal likelihood to
#'   analyse aggregate binary and continuous outcome data, respectively.
#'   \code{\link{prepare_model}} uses the consistency model (as described in
#'   Lu and Ades (2006)) to estimate all possible comparisons in the network.
#'   It also accounts for the multi-arm trials by assigning conditional
#'   univariate normal distributions on the underlying trial-specific effect
#'   size of comparisons with the baseline arm of the multi-arm trial
#'   (Dias et al., 2013a).
#'
#'   The code of Dias et al. (2013a) has been extended to incorporate the
#'   pattern-mixture model to adjust the underlying outcome in each arm of
#'   every trial for missing participant outcome data (Spineli et al., 2021;
#'   Spineli, 2019a; Turner et al., 2015). The assumptions about the
#'   missingness parameter are specified using the arguments \code{mean_misspar}
#'   and \code{var_misspar}. Specifically, \code{run_model} considers the
#'   informative missingness odds ratio in the logarithmic scale for binary
#'   outcome data (Spineli, 2019a; Turner et al., 2015; White et al., 2008), the
#'   informative missingness difference of means when \code{measure} is
#'   \code{"MD"} or \code{"SMD"}, and the informative missingness ratio of means
#'   in the logarithmic scale when \code{measure} is \code{"ROM"}
#'   (Spineli et al., 2021; Mavridis et al., 2015).
#'
#'   When \code{assumption} is trial-specific (i.e., \code{"IDE-TRIAL"} or
#'   \code{"HIE-TRIAL"}), or independent (i.e., \code{"IND-CORR"} or
#'   \code{"IND-UNCORR"}), only one numeric value can be assigned to
#'   \code{mean_misspar} because the same missingness scenario is applied to all
#'   trials and trial-arms of the dataset, respectively. When \code{assumption}
#'   is \code{"IDE-ARM"} or \code{"HIE-ARM"}, a maximum of two
#'   \emph{different} or \emph{identical} numeric values can be assigned as a
#'   vector to \code{mean_misspars}: the first value refers to the experimental
#'   arm, and the second value refers to the control arm of a trial.
#'   In the case of a network, the first value is considered for all
#'   non-reference interventions and the second value is considered for the
#'   reference intervention of the network (i.e., the intervention with
#'   identifier equal to \code{ref}). This is necessary to ensure transitivity
#'   in the assumptions for the missingness parameter across the network
#'   (Spineli, 2019b).
#'
#'   When there is at least one trial-arm with unreported missing participant
#'   outcome data (i.e., \code{m} equals \code{NA} for the corresponding
#'   trial-arms) or when missing participant outcome data have not been
#'   collected for the analysed outcome (i.e., \code{m} is missing in
#'   \code{data}), \code{run_model} assigns the assumption \code{"IND-UNCORR"}
#'   to \code{assumption}.
#'
#'   Currently, there are no empirically-based prior distributions for the
#'   informative missingness parameters. The user may refer to Spineli (2019),
#'   Turner et al. (2015), Mavridis et al. (2015), and White et al. (2008) to
#'   determine \code{mean_misspar} and select a proper value for
#'   \code{var_misspar}.
#'
#'   The scalar \code{base_risk} refers to a fixed baseline risk for the
#'   selected reference intervention (as specified with \code{ref}).
#'   When \code{base_risk} is a three-element vector, it refers to a random
#'   baseline risk and the elements should be sorted in ascending order as they
#'   refer to the lower bound, mean value, and upper bound of the 95\%
#'   confidence interval for the baseline risk for the selected reference
#'   intervention. The \code{\link{baseline_model}} function is called to
#'   calculate the mean and variance of the approximately normal distribution of
#'   the logit of an event for \code{ref} using these three elements
#'   (Dias et al., 2018). When \code{base_risk} is a matrix, it refers to the
#'   predicted baseline risk with first column being the number of events, and
#'   second column being the sample size of the corresponding trials on the
#'   selected reference intervention. Then the \code{\link{baseline_model}}
#'   function is called that contains the WinBUGS code as written by Dias et al.
#'   (2013b) for the hierarchical baseline model. The posterior mean and
#'   precision of the predictive distribution of the logit of an event
#'   for the selected reference intervention are plugged in the WinBUGS code for
#'   the relative effects model (via the \code{\link{prepare_model}} function).
#'   The matrix \code{base_risk} should not comprise the trials in \code{data}
#'   that include the \code{ref}, unless justified (Dias et al., 2018).
#'
#'   To obtain unique absolute risks for each intervention, the network
#'   meta-analysis model has been extended to incorporate the transitive risks
#'   framework, namely, an intervention has the same absolute risk regardless of
#'   the comparator intervention(s) in a trial (Spineli et al., 2017).
#'   The absolute risks are a function of the odds ratio (the \strong{base-case}
#'   effect measure for a binary outcome) and the selected baseline risk for the
#'   reference intervention (\code{ref}) (Appendix in Dias et al., 2013a).
#'   We advocate using the odds ratio as an effect measure for its desired
#'   mathematical properties. Then, the relative risk and risk difference can be
#'   obtained as a function of the absolute risks of the corresponding
#'   interventions in the comparison of interest. Hence, regardless of the
#'   selected \code{measure} for a binary outcome, \code{run_model} performs
#'   pairwise or network meta-analysis based on the odds ratio.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{baseline_model}}, \code{\link{data_preparation}},
#'   \code{\link{heterogeneity_param_prior}}, \code{\link[R2jags:jags]{jags}},
#'   \code{\link{missingness_param_prior}}, \code{\link{prepare_model}}
#'
#' @references
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing
#' between-study heterogeneity and inconsistency in mixed treatment
#' comparisons: Application to stroke prevention treatments in individuals
#' with non-rheumatic atrial fibrillation.
#' \emph{Stat Med} 2009;\bold{28}(14):1861--81. doi: 10.1002/sim.3594
#'
#' Dias S, Ades AE, Welton NJ, Jansen JP, Sutton AJ. Network Meta-Analysis for
#' Decision Making. Chichester (UK): Wiley; 2018.
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013a;\bold{33}(5):607--17. doi: 10.1177/0272989X12458724
#'
#' Dias S, Welton NJ, Sutton AJ, Ades AE. Evidence synthesis for decision
#' making 5: the baseline natural history model. \emph{Med Decis Making}
#' 2013b;\bold{33}(5):657--70. doi: 10.1177/0272989X13485155
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}(4):457--72.
#' doi: 10.1214/ss/1177011136
#'
#' Lu G, Ades AE. Assessing evidence inconsistency in mixed treatment
#' comparisons. \emph{J Am Stat Assoc} 2006;\bold{101}:447--59.
#' doi: 10.1198/016214505000001302
#'
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
#' uncertainty due to missing continuous outcome data in pairwise and
#' network meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--41.
#' doi: 10.1002/sim.6365
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021;\bold{30}(4):958--75.
#' doi: 10.1177/0962280220983544
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019a;\bold{19}(1):86.
#' doi: 10.1186/s12874-019-0731-y
#'
#' Spineli LM. Modeling missing binary outcome data while preserving
#' transitivity assumption yielded more credible network meta-analysis
#' results. \emph{J Clin Epidemiol} 2019b;\bold{105}:19--26.
#' doi: 10.1016/j.jclinepi.2018.09.002
#'
#' Spineli LM, Brignardello-Petersen R, Heen AF, Achille F, Brandt L,
#' Guyatt GH, et al. Obtaining absolute effect estimates to facilitate shared
#' decision making in the context of multiple-treatment comparisons.
#' Abstracts of the Global Evidence Summit, Cape Town, South Africa.
#' \emph{Cochrane Database of Systematic Reviews} 2017;\bold{9}(Suppl 1):18911.
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
#' for uncertainty due to missing binary outcome data in pairwise
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(12):2062--80.
#' doi: 10.1002/sim.6475
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
#' in meta-analysis--part 1: two-stage methods. \emph{Stat Med}
#' 2008;\bold{27}(5):711--27. doi: 10.1002/sim.3008
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Show the first six trials of the dataset
#' head(nma.baker2009)
#'
#' \donttest{
#' # Perform a random-effects network meta-analysis
#' # Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
#' run_model(data = nma.baker2009,
#'           measure = "OR",
#'           model = "RE",
#'           assumption = "IDE-ARM",
#'           heter_prior = list("halfnormal", 0, 1),
#'           mean_misspar = c(0, 0),
#'           var_misspar = 1,
#'           D = 0,
#'           ref = 1,
#'           n_chains = 3,
#'           n_iter = 1000,
#'           n_burnin = 100,
#'           n_thin = 1)
#' }
#'
#' @export
run_model <- function(data,
                      measure,
                      model,
                      assumption,
                      heter_prior,
                      mean_misspar,
                      var_misspar,
                      D,
                      ref,
                      base_risk,
                      n_chains,
                      n_iter,
                      n_burnin,
                      n_thin) {


  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)

  # Missing and default arguments
  model <- if (missing(model)) {
    "RE"
  } else if (!is.element(model, c("RE", "FE"))) {
    stop("Insert 'RE', or 'FE'.", call. = FALSE)
  } else {
    model
  }
  assumption <- if (missing(assumption) & min(na.omit(unlist(item$I))) == 1) {
    message("The 'IDE-ARM' has been used as the default.")
    "IDE-ARM"
  } else if ((missing(assumption) || assumption != "IND-UNCORR") &
             (min(na.omit(unlist(item$I))) == 0 &
              max(na.omit(unlist(item$I))) == 1)) {
    aa <- "Partially extracted missing participants:"
    message(paste(aa, "the 'IND-UNCORR' has been used."))
    "IND-UNCORR"
  } else if ((missing(assumption) || assumption != "IND-UNCORR") &
             (min(na.omit(unlist(item$I))) == 0 &
              max(na.omit(unlist(item$I))) == 0)) {
    "IND-UNCORR"
  } else {
    assumption
  }
  heterog_prior <- heterogeneity_param_prior(measure, model, heter_prior)
  mean_misspar <- if (missing(assumption)  & min(na.omit(unlist(item$I))) == 1) {
    0
  } else if ((missing(assumption) || assumption != "IND-UNCORR") &
             (min(na.omit(unlist(item$I))) == 0 &
              max(na.omit(unlist(item$I))) == 1)) {
    0
  } else if ((missing(assumption) || assumption != "IND-UNCORR") &
             (min(na.omit(unlist(item$I))) == 0 &
              max(na.omit(unlist(item$I))) == 0)) {
    0
  } else {
    missingness_param_prior(assumption, mean_misspar)
  }
  var_misspar <- if (missing(var_misspar) &
                     is.element(measure, c("OR", "RR", "RD", "MD", "SMD"))) {
    1
  } else if (missing(var_misspar) & measure == "ROM") {
    0.2^2
  } else if (var_misspar <= 0) {
    stop("The argument 'var_misspar' must be a positive non-zero number.",
         call. = FALSE)
  } else {
    var_misspar
  }
  D <- if (missing(D)) {
    stop("The argument 'D' needs to be defined.", call. = FALSE)
  } else if (D != 0 & D != 1) {
    stop("The argument 'D' must be '0' or '1'.", call. = FALSE)
  } else {
    D
  }
  ref <- if (missing(ref)) {
    1
  } else if (ref < 1 || ref > item$nt) {
    stop(paste("The argument 'ref' must be an integer from 1 to",
               paste0(item$nt, ".")),
         call. = FALSE)
  } else {
    ref
  }
  ref_base <- if (is.element(measure, c("OR", "RR", "RD")) &
                  missing(base_risk)) {
    base_risk <-
      describe_network(data = data,
                       drug_names = 1:item$nt,
                       measure = measure)$table_interventions[ref, 7]/100
    rep(log(base_risk / (1 - base_risk)), 2)
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    baseline_model(base_risk,
                   n_chains,
                   n_iter,
                   n_burnin,
                   n_thin)$ref_base
  } else if (!is.element(measure, c("OR", "RR", "RD"))) {
    NA
  }
  n_chains <- if (missing(n_chains)) {
    2
  } else if (n_chains < 1) {
    stop("The argument 'n_chains' must be a positive integer.", call. = FALSE)
  } else {
    n_chains
  }
  n_iter <- if (missing(n_iter)) {
    10000
  } else if (n_iter < 1) {
    stop("The argument 'n_iter' must be a positive integer.", call. = FALSE)
  } else {
    n_iter
  }
  n_burnin <- if (missing(n_burnin)) {
    1000
  } else if (n_burnin < 1) {
    stop("The argument 'n_burnin' must be a positive integer.", call. = FALSE)
  } else {
    n_burnin
  }
  n_thin <- if (missing(n_thin)) {
    1
  } else if (n_thin < 1) {
    stop("The argument 'n_thin' must be a positive integer.", call. = FALSE)
  } else {
    n_thin
  }

  # Sign of basic parameters in relation to 'ref'
  indic <- matrix(NA, nrow = item$ns, ncol = max(item$na))
  for (i in 1:item$ns) {
    for (k in 1:max(item$na)) {
      indic[i, k] <- if (item$t[i, k] < ref & !is.na(item$t[i, k])) {
        -1
      } else if (item$t[i, k] >= ref & !is.na(item$t[i, k])) {
        1
      } else if (is.na(item$t[i, k])) {
        NA
      }
    }
  }

  # Data in list format for R2jags
  data_jag <- list("m" = item$m,
                   "N" = item$N,
                   "t" = item$t,
                   "na" = item$na,
                   "nt" = item$nt,
                   "ns" = item$ns,
                   "ref" = ref,
                   "I" = item$I,
                   "indic" = indic,
                   "D" = D)

  data_jag <- if (is.element(measure, c("MD", "SMD", "ROM"))) {
    append(data_jag, list("y.o" = item$y0,
                          "se.o" = item$se0))
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    append(data_jag, list("r" = item$r,
                          "ref_base" = ref_base))
  }

  data_jag <- if (is.element(assumption, "IND-CORR")) {
    append(data_jag, list("M" = ifelse(!is.na(item$m), mean_misspar, NA),
                          "cov.phi" = 0.5 * var_misspar,
                          "var.phi" = var_misspar))
  } else {
    append(data_jag, list("meand.phi" = mean_misspar,
                          "precd.phi" = 1 / var_misspar))
  }

  data_jag <- if (model == "RE") {
    append(data_jag, list("heter.prior" = heterog_prior))
  } else {
    data_jag
  }

  param_jags <- c("EM",
                  "SUCRA",
                  "effectiveness",
                  "dev.o",
                  "totresdev.o",
                  "hat.par",
                  "EM.pred",
                  "tau",
                  "delta")

  param_jags <- if (is.element(assumption,
                               c("HIE-COMMON", "HIE-TRIAL", "HIE-ARM"))) {
    append(param_jags, "mean.phi")
  } else {
    append(param_jags, "phi")
  }

  param_jags <- if (model == "RE" & !is.element(measure, c("RR", "RD"))) {
    param_jags
  } else if (model == "RE" & is.element(measure, c("RR", "RD"))) {
    append(param_jags, "EM.pred.LOR")
  } else if (model == "FE") {
    param_jags[!is.element(param_jags,
                           c("EM.pred", "tau", "delta"))]
  }

  param_jags <- if (is.element(measure, c("OR", "RR", "RD"))) {
    append(param_jags, "abs_risk")
  } else {
    param_jags
  }

  param_jags <- if (is.element(measure, c("RR", "RD"))) {
    append(param_jags, c("SUCRA.LOR", "EM.LOR"))
  } else {
    param_jags
  }

  # Run the Bayesian analysis
  jagsfit <- jags(data = data_jag,
                  parameters.to.save = param_jags,
                  model.file = textConnection(
                    prepare_model(measure,
                                  model,
                                  covar_assumption = "NO",
                                  assumption)
                    ),
                  n.chains = n_chains,
                  n.iter = n_iter,
                  n.burnin = n_burnin,
                  n.thin = n_thin)

  # Turn R2jags object into a data-frame
  get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(get_results %>% dplyr::select(starts_with("EM[")))

  # Predictive effects of all unique pairwise comparisons
  EM_pred <- t(get_results %>% dplyr::select(starts_with("EM.pred[")))

  # Unique absolute risks for all interventions (only binary data)
  abs_risk <- t(get_results %>% dplyr::select(starts_with("abs_risk[")))

  # Estimated og odds ratio of all unique pairwise comparisons
  # (when RR and RD have been selected as effect measures)
  EM_LOR <- t(get_results %>% dplyr::select(starts_with("EM.LOR[")))

  # Predicted log odds ratio of all unique pairwise comparisons
  # (when RR and RD have been selected as effect measures)
  EM_pred_LOR <- t(get_results %>% dplyr::select(starts_with("EM.pred.LOR[")))

  # Between-trial standard deviation
  tau <- t(get_results %>% dplyr::select(starts_with("tau")))

  # SUrface under the Cumulative RAnking curve values
  SUCRA <- t(get_results %>% dplyr::select(starts_with("SUCRA[")))

  # SUrface under the Cumulative RAnking curve values
  # (when RR and RD have been selected as effect measures)
  SUCRA_LOR <- t(get_results %>% dplyr::select(starts_with("SUCRA.LOR[")))

  # Within-trial effects size
  delta <- t(get_results %>% dplyr::select(starts_with("delta") &
                                             !ends_with(",1]")))

  # Ranking probability of each intervention for every rank
  effectiveness <- t(get_results %>% dplyr::select(
    starts_with("effectiveness")))

  # Estimated missingness parameter
  phi <- if (min(na.omit(unlist(item$I))) == 1 &
             max(na.omit(unlist(item$I))) == 1) {
    t(get_results %>% dplyr::select(starts_with("phi") |
                                      starts_with("mean.phi") |
                                      starts_with("mean.phi[") |
                                      starts_with("phi[")))
  } else if (min(na.omit(unlist(item$I))) == 0 &
             max(na.omit(unlist(item$I))) == 1) {
    t(get_results %>% dplyr::select(starts_with("phi[")))
  } else if (min(na.omit(unlist(item$I))) == 0 &
             max(na.omit(unlist(item$I))) == 0) {
    NULL
  }

  # Trial-arm deviance contribution for observed outcome
  dev_o <- t(get_results %>% dplyr::select(starts_with("dev.o")))

  # Fitted/predicted outcome
  hat_par <- t(get_results %>% dplyr::select(starts_with("hat.par")))

  # Total residual deviance
  dev <- jagsfit$BUGSoutput$summary["totresdev.o", "mean"]

  # Calculate the deviance at posterior mean of fitted values
  # Turn 'N' and 'm' into a vector (first column, followed by second, etc)
  m_new <- suppressMessages({
    as.vector(na.omit(melt(item$m)[, 2]))
    })
  N_new <- suppressMessages({
    as.vector(na.omit(melt(item$N)[, 2]))
    })
  obs <- N_new - m_new

  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Turn 'y0', 'se0' into a vector as above
    y0_new <- suppressMessages({
      as.vector(na.omit(melt(item$y0)[, 2]))
      })
    se0_new <- suppressMessages({
      as.vector(na.omit(melt(item$se0)[, 2]))
      })

    # Deviance contribution at the posterior mean of the fitted mean outcome
    dev_post_o <- (y0_new -
                     as.vector(hat_par[, 1])) *
      (y0_new - as.vector(hat_par[, 1])) * (1 / se0_new^2)

    # Sign of the difference between observed and fitted mean outcome
    sign_dev_o <- sign(y0_new - as.vector(hat_par[, 1]))
  } else {
    # Turn 'r' and number of observed into a vector as above
    r_new <- suppressMessages({
      as.vector(na.omit(melt(item$r)[, 2]))
      })

    # Correction for zero events in trial-arm
    r0 <- ifelse(r_new == 0, r_new + 0.01,
                 ifelse(r_new == obs, r_new - 0.01, r_new))

    # Deviance contribution at the posterior mean of the fitted response
    dev_post_o <- 2 * (r0 * (log(r0) -
                               log(as.vector(hat_par[, 1]))) +
                         (obs - r0) * (log(obs - r0) -
                                         log(obs - as.vector(hat_par[, 1]))))

    # Sign of the difference between observed and fitted response
    sign_dev_o <- sign(r0 - as.vector(hat_par[, 1]))
  }

  # Number of unconstrained data-points
  n_data <- sum(item$na)

  # Obtain the leverage for observed outcomes
  leverage_o <- as.vector(dev_o[, 1]) - dev_post_o

  # Number of effective parameters
  pD <- dev - sum(dev_post_o)

  # Deviance information criterion
  DIC <- pD + dev

  # A data-frame on the measures of model assessment:
  # DIC, pD, total residual deviance, and number of unconstrained data-points
  model_assessment <- data.frame(DIC, pD, dev, n_data)

  # Return a list of results
  results <- list(EM = EM,
                  dev_o = dev_o,
                  hat_par = hat_par,
                  leverage_o = leverage_o,
                  sign_dev_o = sign_dev_o,
                  phi = phi,
                  model_assessment = model_assessment,
                  data = data,
                  measure = measure,
                  model = model,
                  assumption = assumption,
                  mean_misspar = mean_misspar,
                  var_misspar = var_misspar,
                  D = D,
                  ref = ref,
                  indic = indic,
                  jagsfit = jagsfit,
                  n_chains = n_chains,
                  n_iter = n_iter,
                  n_burnin = n_burnin,
                  n_thin = n_thin,
                  type = "nma")
  if (model == "RE" & !is.element(measure, c("OR", "RR", "RD"))) {
    ma_results <- append(results, list(EM_pred = EM_pred,
                                       tau = tau,
                                       delta = delta,
                                       heter_prior = heterog_prior))
    nma_results <- append(results, list(EM_pred = EM_pred,
                                        tau = tau,
                                        delta = delta,
                                        heter_prior = heterog_prior,
                                        SUCRA = SUCRA,
                                        effectiveness = effectiveness))
  } else if (model == "RE" & measure == "OR") {
    ma_results <- append(results, list(EM_pred = EM_pred,
                                       tau = tau,
                                       delta = delta,
                                       heter_prior = heterog_prior,
                                       abs_risk = abs_risk,
                                       base_risk = base_risk))
    nma_results <- append(results, list(EM_pred = EM_pred,
                                        tau = tau,
                                        delta = delta,
                                        heter_prior = heterog_prior,
                                        SUCRA = SUCRA,
                                        effectiveness = effectiveness,
                                        abs_risk = abs_risk,
                                        base_risk = base_risk))
  } else if (model == "RE" & is.element(measure, c("RR", "RD"))) {
    ma_results <- append(results, list(EM_pred = EM_pred,
                                       EM_LOR = EM_LOR,
                                       EM_pred_LOR = EM_pred_LOR,
                                       tau = tau,
                                       delta = delta,
                                       heter_prior = heterog_prior,
                                       abs_risk = abs_risk,
                                       base_risk = base_risk))
    nma_results <- append(results, list(EM_pred = EM_pred,
                                        EM_LOR = EM_LOR,
                                        EM_pred_LOR = EM_pred_LOR,
                                        tau = tau,
                                        delta = delta,
                                        heter_prior = heterog_prior,
                                        SUCRA = SUCRA,
                                        SUCRA_LOR = SUCRA_LOR,
                                        effectiveness = effectiveness,
                                        abs_risk = abs_risk,
                                        base_risk = base_risk))
  } else if (model == "FE" & !is.element(measure, c("OR", "RR", "RD"))) {
    ma_results <- results
    nma_results <- append(results, list(SUCRA = SUCRA,
                                        effectiveness = effectiveness))
  } else if (model == "FE" & measure == "OR") {
    ma_results <- append(results, list(abs_risk = abs_risk,
                                       base_risk = base_risk))
    nma_results <- append(results, list(SUCRA = SUCRA,
                                        effectiveness = effectiveness,
                                        abs_risk = abs_risk,
                                        base_risk = base_risk))
  } else if (model == "FE" & is.element(measure, c("RR", "RD"))) {
    ma_results <- append(results, list(EM_LOR = EM_LOR,
                                       abs_risk = abs_risk,
                                       base_risk = base_risk))
    nma_results <- append(results, list(EM_LOR = EM_LOR,
                                        SUCRA = SUCRA,
                                        SUCRA_LOR = SUCRA_LOR,
                                        effectiveness = effectiveness,
                                        abs_risk = abs_risk,
                                        base_risk = base_risk))
  }

  ifelse(item$nt > 2, return(nma_results), return(ma_results))
}
