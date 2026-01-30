# Plot Gower's disimilarity values for each study (Transitivity evaluation)

Illustrating the range of Gower's dissimilarity values for each study in
the network, as well as their between- and within-comparison
dissimilarities

## Usage

``` r
plot_study_dissimilarities(
  results,
  axis_title_size = 12,
  axis_text_size = 12,
  strip_text_size = 11,
  label_size = 3.5
)
```

## Arguments

- results:

  An object of S3 class
  [`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md).
  See 'Value' in
  [`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md).

- axis_title_size:

  A positive integer for the font size of axis title (both axes).
  `axis_title_size` determines the axis.title argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_text_size:

  A positive integer for the font size of axis text (both axes).
  `axis_text_size` determines the axis.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- strip_text_size:

  A positive integer for the font size of facet labels.
  `strip_text_size` determines the strip.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- label_size:

  A positive integer for the font size of labels appearing on each
  study-specific segment. `label_size` determines the size argument
  found in the geom's aesthetic properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

## Value

A horizontal bar plot illustrating the range of Gower's dissimilarity
values for each study with those found in other comparisons. The study
names appear on the y-axis in the order they appear in `results` and the
dissimilarity values appear on the x-axis. Red and blue points refer to
the (average) within-comparison and between-comparison dissimilarity,
respectively, for each study.

A data-frame on the (average) within-comparison and between-comparison
dissimilarities for each study alongside the study name and comparison.
The last two columns refer to the within-comparison and
between-comparison dissimilarities, respectively, after replacing with
the maximum value in the multi-arm trials. These two columns should be
used as a covariate in the function
[`study_perc_contrib`](https://loukiaspin.github.io/rnmamod/reference/study_perc_contrib.md)
to obtain the percentage contribution of each study based on the
covariate values.

## Details

The range of Gower's dissimilarity values for each study versus the
remaining studies in the network for a set of clinical and
methodological characteristics that may act as effect modifiers. Gower's
dissimilarities take values from 0 to 1, with 0 and 1 implying perfect
similarity and perfect dissimilarity, respectively.

The unique dissimilarity values appear as dotted, vertical, grey lines
on each study

## References

Gower J. General Coefficient of Similarity and Some of Its Properties.
*Biometrics* 1971;**27**(4):857–71. doi: 10.2307/2528823

## See also

[`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md),
[`study_perc_contrib`](https://loukiaspin.github.io/rnmamod/reference/study_perc_contrib.md)

## Author

Loukia M. Spineli

## Examples

``` r
# \donttest{
# Fictional dataset
data_set <- data.frame(Trial_name = paste("study", as.character(1:7)),
                      arm1 = c("1", "1", "1", "1", "1", "2", "2"),
                      arm2 = c("2", "2", "2", "3", "3", "3", "3"),
                      sample = c(140, 145, 150, 40, 45, 75, 80),
                      age = c(18, 18, 18, 48, 48, 35, 35),
                      blinding = factor(c("yes", "yes", "yes", "no", "no", "no", "no")))

# Obtain comparison dissimilarities (informative = TRUE)
res <- comp_clustering(input = data_set,
                       drug_names = c("A", "B", "C"),
                       threshold = 0.13,  # General research setting
                       informative = TRUE,
                       get_plots = TRUE)
#> - 3 observed comparisons (0 single-study comparisons)
#> - Dropped characteristics: none
#> $Trials_diss_table
#>             study 1 B-A study 2 B-A study 3 B-A study 4 C-A study 5 C-A
#> study 1 B-A       0.000          NA          NA          NA          NA
#> study 2 B-A       0.015       0.000          NA          NA          NA
#> study 3 B-A       0.030       0.015       0.000          NA          NA
#> study 4 C-A       0.970       0.985       1.000       0.000          NA
#> study 5 C-A       0.955       0.970       0.985       0.015       0.000
#> study 6 C-B       0.719       0.734       0.749       0.251       0.235
#> study 7 C-B       0.704       0.719       0.734       0.266       0.251
#>             study 6 C-B study 7 C-B
#> study 1 B-A          NA          NA
#> study 2 B-A          NA          NA
#> study 3 B-A          NA          NA
#> study 4 C-A          NA          NA
#> study 5 C-A          NA          NA
#> study 6 C-B       0.000          NA
#> study 7 C-B       0.015           0
#> 
#> $Comparisons_diss_table
#>      B-A  C-A  C-B
#> B-A 0.02   NA   NA
#> C-A 0.98 0.02   NA
#> C-B 0.73 0.25 0.02
#> 
#> $Total_dissimilarity
#>   comparison total_dissimilarity         index_type
#> 5 C-A vs C-B                0.25 Between-comparison
#> 3 B-A vs C-B                0.73 Between-comparison
#> 2 B-A vs C-A                0.98 Between-comparison
#> 1        B-A                0.02  Within-comparison
#> 4        C-A                0.02  Within-comparison
#> 6        C-B                0.02  Within-comparison
#> 
#> $Types_used
#>   characteristic    type
#> 1         sample  double
#> 2            age  double
#> 3       blinding integer
#> 
#> $Total_missing
#> [1] "0%"
#> 
#> $Within_comparison_dissimilarity

#> 
#> $Between_comparison_dissimilarity

#> 
#> $Dissimilarity_heatmap

#> 
#> attr(,"class")
#> [1] "comp_clustering"

plot_study_dissimilarities(results = res,
                           axis_title_size = 12,
                           axis_text_size = 12,
                           strip_text_size = 11,
                           label_size = 3.5)
#> [[1]]

#> 
#> $diss_values
#>   study_id   study comp within_value between_value within_multiarm
#> 1        1 study 1  B-A   0.02371708     0.8463897      0.02371708
#> 2        2 study 2  B-A   0.01500000     0.8612262      0.01500000
#> 3        3 study 3  B-A   0.02371708     0.8760682      0.02371708
#> 4        4 study 4  C-A   0.01500000     0.7803694      0.01500000
#> 5        5 study 5  C-A   0.01500000     0.7669910      0.01500000
#> 6        6 study 6  C-B   0.01500000     0.5890576      0.01500000
#> 7        7 study 7  C-B   0.01500000     0.5805325      0.01500000
#>   between_multiarm
#> 1        0.8463897
#> 2        0.8612262
#> 3        0.8760682
#> 4        0.7803694
#> 5        0.7669910
#> 6        0.5890576
#> 7        0.5805325
#> 
# }
```
