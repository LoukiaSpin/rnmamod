

data <- nma.baker2009
measure <-"OR"
mean.value <- 0
var.value <- 1
rho <- 0


unrelated.effects.plot <- function(data, measure, mean.value, var.value, rho) {


  # Default arguments
  mean.value <- ifelse(missing(mean.value), 0, mean.value)  # MAR assumption
  var.value <- ifelse(missing(var.value) & (is.element(measure, c("OR", "MD", "SMD"))), 1, ifelse(missing(var.value) & measure == "ROM", 0.2^2, var.value))
  rho <- ifelse(missing(rho), 0, rho)                       # Uncorrelated within-trial missingness parameters


  item <- data.preparation(data, measure)


  ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    pairwise.observed <- pairwise(as.list(item$t),
                                  mean = as.list(item$y0),
                                  sd = as.list(item$sd0),
                                  n = as.list(item$N),
                                  data = data,
                                  studlab = 1:item$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    colnames(pairwise.observed) <- c("study", "arm1", "arm2", "y1", "y2", "sd1", "sd2", "n1", "n2")

    pairwise.mod <- pairwise(as.list(item$t),
                             mean = as.list(item$y0),
                             sd = as.list(item$sd0),
                             n = as.list(item$m),
                             data = data,
                             studlab = 1:item$ns)[, c(6, 9)]
    colnames(pairwise.mod) <- c("m1", "m2")
  } else {
    pairwise.observed <- pairwise(as.list(item$t),
                                  event = as.list(item$r),
                                  n = as.list(item$N),
                                  data = data,
                                  studlab = 1:item$ns)[, c(3:6, 8, 7, 9)]
    colnames(pairwise.observed) <- c("study", "arm1", "arm2", "r1", "r2", "n1", "n2")

    pairwise.mod <- pairwise(as.list(item$t),
                             event = as.list(item$m),
                             n = as.list(item$N),
                             data = data,
                             studlab = 1:item$ns)[, c(6, 8)]
    colnames(pairwise.mod) <- c("m1", "m2")
  }


  ## The dataset to perform the fixed-effects analysis (i.e. unrelated trial effects)
  pairwise.data <- data.frame(pairwise.observed, pairwise.mod)


  contrast0 <- if (is.element(measure, c("MD", "SMD", "ROM"))) {
    Taylor.IMDoM.IMRoM(pairwise.data, measure, mean.value, var.value, rho)
  } else {
    Taylor.IMOR(pairwise.data, measure, mean.value, var.value, rho)
  }

  contrast <- contrast0[, c(1, 16:19)]
  contrast$LROM <- ifelse(contrast$t2 == "Control", (-1)*contrast$LROM,  contrast$LROM)

  (contrast$lower <- contrast$LROM - 1.95*contrast$SE.LROM)
  (contrast$upper <- contrast$LROM + 1.95*contrast$SE.LROM)
  (contrast$studlab <- sub("_", " ", rep(primary[, 1], na)))
  contrast$ROM <- exp(ifelse(contrast$t2 == "Control", (-1)*contrast$LROM,  contrast$LROM)) # EDW!!!!
  (contrast$ROM.lower <- exp(contrast$LROM - 1.95*contrast$SE.LROM))
  (contrast$ROM.upper <- exp(contrast$LROM + 1.95*contrast$SE.LROM))
  (contrast$comp <- paste(contrast$t2, "versus", contrast$t1))
  (contrast$design <- rep(primary[, 2], na))
  (contrast$bcount <- rep(primary[, 4], na))
  (contrast$rob <- rep(primary[, 37], na))




}


ggplot(contrast, aes(x = exp(LROM), y = studlab, xmin = exp(lower), xmax = exp(upper), color = design, linetype = bcount, shape = rob)) +
  geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, lty = 1, col = "grey") +
  geom_point(size = 3, color = "black") +
  geom_text(aes(x = exp(LROM), y = studlab, label = round(exp(LROM), 2)), color = "black", hjust = -0.3, vjust = -0.1, size = 3.5,
            check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
  scale_color_manual(breaks = c("NRCT", "RCT"), values = c("grey", "grey47")) +
  facet_wrap(vars(comp), scales = "free_x") +
  labs(y = "", x = "Ratio of ratio of means in total bacterial", color = "Design", linetype = "Bacterial count", shape = " Risk of bias") +
  theme_classic() +
  scale_x_continuous(trans = 'log10') +
  theme(axis.title.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 12),
        strip.text = element_text(size = 12), legend.position = "bottom", legend.text = element_text(color = "black", size = 12), legend.title = element_text(color = "black", size = 12, face = "bold"))
