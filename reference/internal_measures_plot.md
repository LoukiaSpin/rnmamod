# Internal measures for cluster validation (Comparisons' comparability for transitivity evaluation)

`internal_measures_plot` currently prepares the table with the results
of the average silhouette width for a range of clusters, and visualises
the results using a profile plot.

## Usage

``` r
internal_measures_plot(
  input,
  optimal_link,
  label_size = 4,
  axis_title_size = 14,
  axis_text_size = 14
)
```

## Arguments

- input:

  An object of 'dist' class. It is a lower off-diagonal matrix with the
  dissimilarities of all pairs of comparisons.

- optimal_link:

  A character string with values `"ward.D"`, `"ward.D2"`, `"single"`,
  `"complete"`, `"average"`, `"mcquitty"`, `"median"`, or `"centroid"`
  for the optimal linkage method, corresponding to the highest
  cophenetic correlation coefficient value.

- label_size:

  A positive integer for the font size of labels in the profile plot
  with the average silhouette width per candidate cluster. `label_size`
  determines the size argument found in the geom's aesthetic properties
  in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_title_size:

  A positive integer for the font size of axis title in the profile plot
  with the average silhouette width per candidate cluster.
  `axis_title_size` determines the axis.title argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_text_size:

  A positive integer for the font size of axis text in the profile plot
  with the average silhouette width per candidate cluster.
  `axis_text_size` determines the axis.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

## Value

`internal_measures_plot` currently returns the following list of
elements:

- Table_internal_measures:

  A data-frame of the average silhouette width for a range of 2 to P-1
  clusters, with P being the number of trials

- Internal_measures_panel:

  A profile plot on the average silhouette width for a range of 2 to P-1
  clusters, with P being the number of trials The candidate optimal
  number of clusters is indicated with a red point directly on the line.

## Details

`internal_measures_plot` also calls the function
[`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md)
to define the argument `optimal_link` to create the silhouette plot for
the selected number of clusters.

`internal_measures_plot` calls the
[`silhouette`](https://rdrr.io/pkg/cluster/man/silhouette.html) function
in the R-package [cluster](https://CRAN.R-project.org/package=cluster)
to obtain the results on average silhouette for each candidate cluster.

`internal_measures_plot` is integrated in the function
[`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md).

## References

Handl J, Knowles J, Kell DB. Computational cluster validation in
post-genomic data analysis. *Biometrics* 2005;**21**(15):3201–120. doi:
10.1093/bioinformatics/bti517

Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
validation of cluster analysis. *J Comput Appl Math* 1987;**20**:53–65.

## See also

[`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md),
[`silhouette`](https://rdrr.io/pkg/cluster/man/silhouette.html)

## Author

Loukia M. Spineli
