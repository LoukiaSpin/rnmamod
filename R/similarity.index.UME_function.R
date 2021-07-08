#' Similarity index: consistency model versus unrelated mean effects model
#'
#' @description This function calculates the Kullback-Leibler divergence (KLD) measure (Kullback and Leibler, 1951) in the summary effect size from the unrelated mean effects model
#'   to the consistency model. Then, for each pairwise comparison observed in the network, the KLD value is compared with a threshold that reflects
#'   the minimally allowed deviation between the compared models.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param ume An object of S3 class \code{\link{run.UME}}. See 'Value' in \code{\link{run.UME}}.
#' @param threshold A number indicating the threshold of similarity. See 'Details' below.
#'
#' @return A list of the following two vectors:
#' \tabular{ll}{
#'  \code{KLD} \tab A numeric vector on the KLD values; one KLD value per observed comparison.\cr
#'  \tab \cr
#'  \code{robust} \tab A character vector on whether the divergence is acceptable or substantial.\cr
#' }
#'
#' In \code{robust}, acceptable divergence is labeled as \code{"robust"} (i.e., \code{KLD} \eqn{<} \code{threshold}), and substantial divergence is labeled as \code{"frail"}. This is equivalent to concluding
#' that the summary effect size of a comparison is robust and sensible, respectively, to the applied model.
#'
#' @details \code{similarity.index.UME} is integrated in the \code{heatmap.similarity.UME} function.
#'   The user may consider the values 0.28 and 0.17 as \code{threshold} for binary and continuous outcome data, (the default values), respectively, or consider other plausible values.
#'   Spineli et al. (2021) offers a discussion on specifying the \code{threshold}.
#'   These thresholds have been originally developed by Spineli et al. (2021) and considered also by Spineli (2021) in the proposed framework of global evaluation of the consistency assumption.
#'
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.UME}}, \code{\link{heatmap.similarity.UME}}
#'
#' @references
#' Spineli LM. A novel framework to evaluate the consistency assumption globally in a network of interventions. \emph{submitted} 2021.
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Kullback S, Leibler RA. On information and sufficiency. \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. [\doi{10.1214/aoms/1177729694}]
#'
#' @export
similarity.index.UME <- function(full, ume, threshold){


  ## Function for the Kullback-Leibler Divergence (comparing two univariate normal distributions)
  KLD.measure.univ <- function(mean.y, sd.y, mean.x, sd.x){

    # x is the 'truth' (e.g. the MAR assumption)
    KLD.xy <- 0.5*(((sd.x/sd.y)^2) + ((mean.y - mean.x)^2)/(sd.y^2) - 1 + 2*log(sd.y/sd.x))

    return(KLD.xy = KLD.xy)
  }


  kldxy <- rep(NA, length(full[, 1]))

  for(i in 1:length(full[, 1])){ ## We are interested in all observed comparisons of the network

    ## Returns the KLD of UME when compared with NMA (NMA as 'true') for comparison i
    kldxy[i] <- KLD.measure.univ(ume[i, 1], ume[i, 2], full[i, 1], full[i, 2])

  }


  robust <- ifelse(kldxy < threshold, "robust", "frail")

  return(list(robust = robust, KLD = kldxy))
}
