design_matrices <- function (study_name,
                             base_t,
                             exp_t,
                             ref_t,
                             covar,
                             covar_assum) {


  ## Find the unique treatments and sort in ascending order
  unique_treats <- sort(unique(c(cbind(base_t, exp_t))))


  ## Default arguments
  # The covariate variable
  covar <- if(missing(covar)) {
    stop("'covar' must be defined", call. = FALSE)
  } else {
    covar
  }
  # Create a centered variable if 'covar' is metric
  covar.val <-
    if (length(unique(covar)) > 2) { # Center metric variable
      covar - mean(covar)
    } else if (length(unique(covar)) == 2) { # Binary
      covar
    }
  # Treatment-by-covariate Interaction assumption
  covar_assum <- if (missing(covar_assum)) {
    stop("'covar_assum' must be defined", call. = FALSE)
  } else if (!is.element(covar_assum,
                         c("no", "common", "exchangeable", "independent"))) {
    aa <- "'common', 'exchangeable', or 'independent'."
    stop(paste("'covar_assum' must be any of the following: 'no',", aa),
         call. = FALSE)
  } else {
    covar_assum
  }
  # The reference treatment
  ref_t <- if (!is.element(ref_t, unique_treats)) {
    stop("'ref_t' must be any of the treatments in the dataset", call. = FALSE)
  } else {
    ref_t
  }
  # Rename only multi-arm trials
  study_name_new <-
    ave(study_name, study_name,
        FUN = function(x) if (length(x) > 1) paste0(x[1], "(", seq_along(x), ")") else x[1])


  ## Turn into factors
  base <- factor(base_t, levels = unique_treats)
  exp <- factor(exp_t, levels = unique_treats)


  ## PART 1: Design matrix X (studies by basic parameters & regression coefficients)
  design_X <- model.matrix(~ exp - 1) - model.matrix(~ base - 1)
  colnames(design_X) <- unique_treats

  if (missing(ref_t)) {
    ref_t <- unique_treats[1]
    design_X[, colnames(design_X) != ref_t]
  }


  ## Define design matrix X based on the 'covar_assum' argument
  # Columns refer to basic parameters, followed by basic regression coefficients (if relevant)
  if (covar_assum == "no") {
    design_X <- design_X[, -ref_t]
    colnames(design_X) <- paste0("d", unique_treats[-ref_t], ref_t)
  } else if (covar_assum == "common") {
    design_mat_X <- cbind(design_X[, -ref_t], design_X[, -ref_t] * covar.val)
    colnames(design_mat_X) <- c(paste0("d", unique_treats[-ref_t], ref_t),
                                paste0("beta", unique_treats[-ref_t], ref_t))
    # Final design matrix X for 'common' interaction
    design_X <-
      data.frame(design_mat_X[, -c(length(unique_treats):(2*(length(unique_treats) - 1)))],
                 beta = apply(design_mat_X[, length(unique_treats):(2*(length(unique_treats) - 1))], 1, sum))

  } else if (!is.element(covar_assum, c("no", "common"))) {
    design_X <- cbind(design_X[, -ref_t], design_X[, -ref_t] * covar.val)
    colnames(design_X) <- c(paste0("d", unique_treats[-ref_t], ref_t),
                            paste0("beta", unique_treats[-ref_t], ref_t))
  }
  rownames(design_X) <- paste0("y", study_name_new)


  ## PART 2: Design matrix Z (all comparisons by basic parameters & regression coefficients)
  poss_comb <- t(combn(unique_treats, 2))

  # Turn into factors
  baseline <- factor(poss_comb[, 1], levels = unique_treats)
  comparator <- factor(poss_comb[, 2], levels = unique_treats)

  # Design part of the matrix Z
  design_mat_Z0 <- (model.matrix(~ comparator - 1) - model.matrix(~ baseline - 1))[, -ref_t]

  # Turn it into block diagonal, the design matrix Z (based on 'covar_assum')
  if (covar_assum == "no") {
    design_Z <- design_mat_Z0
    colnames(design_Z) <- paste0("d", unique_treats[-ref_t], ref_t)
    rownames(design_Z) <- paste0("d", apply(combn(unique_treats, 2), 2,
                                            function(x) paste0((x), collapse = "")))
  } else if (covar_assum == "common") {
    design_mat_Z <- as.matrix(bdiag(design_mat_Z0, design_mat_Z0))
    colnames(design_mat_Z) <- c(paste0("d", unique_treats[-ref_t], ref_t),
                                paste0("beta", unique_treats[-ref_t], ref_t))
    rownames(design_mat_Z) <- c(paste0("d", apply(combn(unique_treats, 2), 2,
                                                  function(x) paste0((x), collapse = ""))),
                                paste0("beta", apply(combn(unique_treats, 2), 2,
                                                     function(x) paste0((x), collapse = ""))))
    # Final design matrix Z for 'common' interaction
    design_Z <-
      design_mat_Z[1:(dim(combn(length(unique_treats), 2))[2] + 1),
                   -c((length(unique_treats) + 1):(2*(length(unique_treats) - 1)))]
  } else if (!is.element(covar_assum, c("no", "common"))) {
    design_Z <- as.matrix(bdiag(design_mat_Z0, design_mat_Z0))
    colnames(design_Z) <- c(paste0("d", unique_treats[-ref_t], ref_t),
                            paste0("beta", unique_treats[-ref_t], ref_t))
    rownames(design_Z) <- c(paste0("d", apply(combn(unique_treats, 2), 2,
                                              function(x) paste0((x), collapse = ""))),
                            paste0("beta", apply(combn(unique_treats, 2), 2,
                                                 function(x) paste0((x), collapse = ""))))
  }

  return(list(design_X = as.matrix(design_X),
              design_Z = as.matrix(design_Z)))
}
