% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kld.inconsistency.user_function.R
\name{kld_inconsistency_user}
\alias{kld_inconsistency_user}
\title{Density plots of local inconsistency results and Kullback-Leibler divergence
 (When dataset is created by the user)}
\usage{
kld_inconsistency_user(
  dataset,
  threshold = 1e-05,
  level = 0.05,
  outcome = NULL,
  scales = "free",
  show_incons = TRUE,
  y_axis_name = TRUE,
  title_name = NULL,
  axis_title_size = 13,
  axis_text_size = 13,
  text_size = 3.5,
  strip_text_size = 13,
  legend_title_size = 13,
  legend_text_size = 13,
  str_wrap_width = 10
)
}
\arguments{
\item{dataset}{A data-frame of seven columns and as many rows as the split
nodes. The first column contains the names of the split nodes, and the
remaining columns have the point estimate and standard error of the direct,
indirect and inconsistency parameter in that order.}

\item{threshold}{A positive number indicating the threshold of not concerning
inconsistency, that is, the minimally allowed deviation between the direct
and indirect estimates for a split node that does raise concerns for
material inconsistency. The argument is optional.}

\item{level}{A number indicating the significance level. Suggested values
are 0.05 and 0.10. The default value is 0.05.}

\item{outcome}{Optional argument to describe the effect measure used (the
x-axis of the plots).}

\item{scales}{A character on whether both axes should be fixed
(\code{"fixed"}) or free (\code{"free"}) or only one of them be free
(\code{"free_x"} or \code{"free_y"}). \code{scales} determines the scales
argument found in function (\code{\link[ggplot2:facet_wrap]{facet_wrap}})
in the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}. The default is
(\code{"free"}).}

\item{show_incons}{Logical to indicate whether to present the point estimate
and 95% interval of the inconsistency parameter. The default is \code{TRUE}
(report).}

\item{y_axis_name}{Logical to indicate whether to present the title of y-axis
('Density'). The default is \code{TRUE} (report).}

\item{title_name}{Text for the title of the plot. \code{title_name}
determines the labs argument of the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.}

\item{axis_title_size}{A positive integer for the font size of axis title.
\code{axis_title_size} determines the axis.title argument found in the
theme's properties in the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
The default option is 13.}

\item{axis_text_size}{A positive integer for the font size of axis text.
\code{axis_text_size} determines the axis.text argument found in the
theme's properties in the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
The default option is 13.}

\item{text_size}{A positive integer for the font size of labels.
\code{text_size} determines the size argument found in the geom_text
function in the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
The default option is 3.5.}

\item{strip_text_size}{A positive integer for the font size of facet labels.
\code{legend_text_size} determines the legend.text argument found in
the theme's properties in the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
The default option is 13.}

\item{legend_title_size}{A positive integer for the font size of legend
title. \code{legend_text_size} determines the legend.text argument found in
the theme's properties in the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
The default option is 13.}

\item{legend_text_size}{A positive integer for the font size of legend text.
\code{legend_text_size} determines the legend.text argument found in the
theme's properties in the R-package
\href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
The default option is 13.}

\item{str_wrap_width}{A positive integer for wrapping the axis labels in the
percent stacked bar-plot. \code{str_wrap_width} determines the
\code{\link[stringr:str_wrap]{str_wrap}} function of the R-package
\href{https://CRAN.R-project.org/package=stringr}{stringr}.}
}
\value{
The first plot is a panel of density plots for each split node sorted
in ascending order of the Kullback-Leibler divergence value. Blue and black
lines refer to the direct and indirect estimates, respectively. The grey
segment refers to the (1 - \code{level})\% 'pseudo' confidence interval of
the inconsistency parameter based on the corresponding normal z-scores,
with a darker grey line  referring to the point estimate. The names of the
selected comparisons appear at the top of each plot. The mean estimate on
the scale of the selected effect measure appears at the top of each density
curve.

The Kullback-Leibler divergence value appears at the top left of each plot
in three colours: black, if no threshold has been defined (the default),
green, if the Kullback-Leibler divergence is below the specified
\code{threshold} (not concerning inconsistency) and red, if the
Kullback-Leibler divergence is at least the specified \code{threshold}
(substantial inconsistency).

The second plot is a percent stacked bar plot on the percentage contribution
of approximating direct with indirect estimate (and vice-versa) to the total
information loss for each target comparison. Total information loss is
defined as the sum of the Kullback-Leibler divergence value
when approximating the direct with indirect estimate (blue bars), and the
Kullback-Leibler divergence value when approximating the indirect
with direct estimate (black bars). Values parentheses refer to the
corresponding Kullback-Leibler divergence value. Bars are sorted in ascending
order of the average Kullback-Leibler divergence value.

The function also returns the data-frame \code{average_KLD} that includes the
split comparisons and the corresponding average Kullback-Leibler divergence
value.
}
\description{
When the user has extracted results obtained from a method of local
inconsistency evaluation (e.g., loop-specific, back-calculation or
node-splitting approaches) as reported in publication, this function
provides the same output with the function \code{\link{kld_inconsistency}}.
A panel of density plots on the direct and indirect estimates of the
selected comparisons based on approach for local inconsistency evaluation,
such as back-calculation and node-splitting approaches (Dias et al., 2010;
van Valkenhoef et al., 2016) and loop-specific approach (Bucher et al., 1997)
accompanied by the average Kullback-Leibler divergence. Additionally, stacked
bar plots on the percentage contribution of either Kullback-Leibler
divergence (from direct to indirect, and vice-versa) to the total information
loss for each selected comparison are presented (Spineli, 2024).
}
\examples{

\dontrun{
## Data are taken from Table II in Dias et al. (2010)
# Treatments compared
treat <-
c("SK", "t-PA", "Acc t-PA", "SK+t-PA", "r-PA", "TNK", "PTCA", "UK", "ASPAC")

# Baseline arm (from each selected comparison)
base <- rep(1:3, c(6, 3, 5))

# Non-baseline arm (from each selected comparison)
nonbase <- c(2, 3, 5, 7, 8, 9, 7, 8, 9, 4, 5, 7, 8, 9)

# Compared treatments with their names
treat_comp <-
mapply(function(x, y) paste(treat[x], "vs", treat[y]), base, nonbase)

# Direct results
direct_mean <- c(0.000, -0.158, -0.060, -0.666, -0.369, 0.009, -0.545,
                 -0.295, 0.006, 0.126, 0.019, -0.216, 0.143, 1.409)
direct_sd <- c(0.030, 0.048, 0.089, 0.185, 0.518, 0.037, 0.417, 0.347, 0.037,
               0.054, 0.066, 0.118, 0.356, 0.415)

# Indirect results
indirect_mean <- c(0.189, -0.247, -0.175, -0.393, -0.168, 0.424, -0.475,
                   -0.144, 0.471, 0.630, 0.135, -0.477, -0.136, 0.165)
indirect_sd <- c(0.235, 0.092, 0.081, 0.120, 0.244, 0.252, 0.108, 0.290,
                 0.241, 0.697, 0.101, 0.174, 0.288, 0.057)

# Inconsistency
incons_mean <- c(-0.190, 0.088, 0.115, -0.272, -0.207, -0.413, -0.073,
                 -0.155, -0.468, -0.506, -0.116, 0.260, 0.277, 1.239)
incons_sd <- c(0.236, 0.104, 0.121, 0.222, 0.575, 0.253, 0.432, 0.452, 0.241,
               0.696, 0.120, 0.211, 0.461, 0.420)

# Collect results in a data-frame (exactly as required from the function)
dias_results <- data.frame(treat_comp, direct_mean, direct_sd, indirect_mean,
                           indirect_sd, incons_mean, incons_sd)

# Apply the function
kld_inconsistency_user(dataset = dias_results,
                       threshold = 0.13,
                       outcome = "Odds ratio (logarithmic scale)")
}

}
\references{
Bucher HC, Guyatt GH, Griffith LE, Walter SD. The results of direct and
indirect treatment comparisons in meta-analysis of randomized controlled
trials. \emph{J Clin Epidemiol} 1997;\bold{50}(6):683--91.

Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
treatment comparison meta-analysis.
\emph{Stat Med} 2010;\bold{29}(7-8):932--44.
doi: 10.1002/sim.3767

Kullback S, Leibler RA. On information and sufficiency.
\emph{Ann Math Stat} 1951;\bold{22}(1):79--86. doi: 10.1214/aoms/1177729694

Spineli LM. Local inconsistency detection using the Kullback-Leibler
divergence measure. \emph{Syst Rev} 2024;\bold{13}(1):261.
doi: 10.1186/s13643-024-02680-4.

van Valkenhoef G, Dias S, Ades AE, Welton NJ. Automated generation of
node-splitting models for assessment of inconsistency in network
meta-analysis. \emph{Res Synth Methods} 2016;\bold{7}(1):80--93.
doi: 10.1002/jrsm.1167
}
\seealso{
\code{\link[ggplot2:facet_wrap]{facet_wrap}},
  \code{\link{kld_inconsistency}},
  \code{\link{kld_measure}}
}
\author{
{Loukia M. Spineli}
}
