# Density plots of local inconsistency results and Kullback-Leibler divergence when 'rnmamod', 'netmeta' or 'gemtc' R packages are used

A panel of density plots on the direct and indirect estimates of the
selected comparisons based on approach for local inconsistency
evaluation, such as back-calculation and node-splitting approaches (Dias
et al., 2010; van Valkenhoef et al., 2016) and loop-specific approach
(Bucher et al., 1997) accompanied by the average Kullback-Leibler
divergence. Additionally, stacked bar plots on the percentage
contribution of either Kullback-Leibler divergence (from direct to
indirect, and vice-versa) to the total information loss for each
selected comparison are presented (Spineli, 2024). The function handles
results also from the R-packages
[gemtc](https://CRAN.R-project.org/package=gemtc) and
[netmeta](https://CRAN.R-project.org/package=netmeta).

## Usage

``` r
kld_inconsistency(
  node,
  threshold = 1e-05,
  drug_names = NULL,
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
```

## Arguments

- node:

  An object of S3 class
  [`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md)
  or class
  [`mtc.nodesplit`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html)
  (see [gemtc](https://CRAN.R-project.org/package=gemtc)) or class
  [`netsplit`](https://rdrr.io/pkg/netmeta/man/netsplit.html) (see
  [netmeta](https://CRAN.R-project.org/package=netmeta)).

- threshold:

  A positive number indicating the threshold of not concerning
  inconsistency, that is, the minimally allowed deviation between the
  direct and indirect estimates for a split node that does raise
  concerns for material inconsistency. The argument is optional.

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data`. It is not relevant for
  [gemtc](https://CRAN.R-project.org/package=gemtc) and
  [netmeta](https://CRAN.R-project.org/package=netmeta).

- outcome:

  Optional argument to describe the effect measure used (the x-axis of
  the plots).

- scales:

  A character on whether both axes should be fixed (`"fixed"`) or free
  (`"free"`) or only one of them be free (`"free_x"` or `"free_y"`).
  `scales` determines the scales argument found in function
  ([`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html))
  in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2). The default is
  (`"free"`).

- show_incons:

  Logical to indicate whether to present the point estimate and 95
  (report).

- y_axis_name:

  Logical to indicate whether to present the title of y-axis
  ('Density'). The default is `TRUE` (report).

- title_name:

  Text for the title of the plot. `title_name` determines the labs
  argument of the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_title_size:

  A positive integer for the font size of axis title. `axis_title_size`
  determines the axis.title argument found in the theme's properties in
  the R-package [ggplot2](https://CRAN.R-project.org/package=ggplot2).
  The default option is 13.

- axis_text_size:

  A positive integer for the font size of axis text. `axis_text_size`
  determines the axis.text argument found in the theme's properties in
  the R-package [ggplot2](https://CRAN.R-project.org/package=ggplot2).
  The default option is 13.

- text_size:

  A positive integer for the font size of labels. `text_size` determines
  the size argument found in the geom_text function in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2). The default
  option is 3.5.

- strip_text_size:

  A positive integer for the font size of facet labels.
  `legend_text_size` determines the legend.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2). The default
  option is 13.

- legend_title_size:

  A positive integer for the font size of legend title.
  `legend_text_size` determines the legend.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2). The default
  option is 13.

- legend_text_size:

  A positive integer for the font size of legend text.
  `legend_text_size` determines the legend.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2). The default
  option is 13.

- str_wrap_width:

  A positive integer for wrapping the axis labels in the percent stacked
  bar-plot. `str_wrap_width` determines the
  [`str_wrap`](https://stringr.tidyverse.org/reference/str_wrap.html)
  function of the R-package
  [stringr](https://CRAN.R-project.org/package=stringr).

## Value

The first plot is a panel of density plots for each split node sorted in
ascending order of the Kullback-Leibler divergence value. Blue and black
lines refer to the direct and indirect estimates, respectively. The grey
segment refers to the 95% credible (confidence) interval of the
inconsistency parameter, when
[`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md)
([`netsplit`](https://rdrr.io/pkg/netmeta/man/netsplit.html)) has been
applied, with a darker grey line referring to the point estimate. When
[`mtc.nodesplit`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html) has
been employed, the 95% confidence interval has been approximated using
the Bucher's approach based on the corresponding direct and indirect
results. This was necessary because
[`mtc.nodesplit`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html)
(version 1.0-2) returns only the inconsistency p-values rather than the
posterior results on the inconsistency parameters. The mean estimate on
the scale of the selected effect measure appears at the top of each
density curve.

The Kullback-Leibler divergence value appears at the top left of each
plot in three colours: black, if no threshold has been defined (the
default), green, if the Kullback-Leibler divergence is below the
specified `threshold` (not concerning inconsistency) and red, if the
Kullback-Leibler divergence is at least the specified `threshold`
(substantial inconsistency).

The second plot is a percent stacked bar plot on the percentage
contribution of approximating direct with indirect estimate (and
vice-versa) to the total information loss for each target comparison.
Total information loss is defined as the sum of the KLD value when
approximating the direct with indirect estimate (blue bars), and the KLD
when approximating the indirect with direct estimate (black bars).
Values parentheses refer to the corresponding KLD value.

The function also returns the data-frame `average_KLD` that includes the
split comparisons and the corresponding average Kullback-Leibler
divergence value.

## References

Bucher HC, Guyatt GH, Griffith LE, Walter SD. The results of direct and
indirect treatment comparisons in meta-analysis of randomized controlled
trials. *J Clin Epidemiol* 1997;**50**(6):683–91.

Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
treatment comparison meta-analysis. *Stat Med* 2010;**29**(7-8):932–44.
doi: 10.1002/sim.3767

Kullback S, Leibler RA. On information and sufficiency. *Ann Math Stat*
1951;**22**(1):79–86. doi: 10.1214/aoms/1177729694

Spineli LM. Local inconsistency detection using the Kullback-Leibler
divergence measure. *Syst Rev* 2024;**13**(1):261. doi:
10.1186/s13643-024-02680-4.

van Valkenhoef G, Dias S, Ades AE, Welton NJ. Automated generation of
node-splitting models for assessment of inconsistency in network
meta-analysis. *Res Synth Methods* 2016;**7**(1):80–93. doi:
10.1002/jrsm.1167

## See also

[`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html),
[`mtc.nodesplit`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html),
[`kld_measure`](https://loukiaspin.github.io/rnmamod/reference/kld_measure.md),
[`netsplit`](https://rdrr.io/pkg/netmeta/man/netsplit.html),
[`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md),
[`str_wrap`](https://stringr.tidyverse.org/reference/str_wrap.html)

## Author

Loukia M. Spineli

## Examples

``` r
if (FALSE) { # \dontrun{
data("nma.baker2009")

# Read results from 'run_nodesplit' (using the default arguments)
node <- readRDS(system.file('extdata/node_baker.rds', package = 'rnmamod'))

# The names of the interventions in the order they appear in the dataset
interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
                  "fluticasone", "fluticasone plus salmeterol",
                  "formoterol", "salmeterol", "tiotropium")

# Apply the function
kld_inconsistency(node = node,
                  threshold = 0.64,
                  drug_names = interv_names,
                  outcome = "Odds ratio (logarithmic scale)",
                  str_wrap_width = 15)
} # }
```
