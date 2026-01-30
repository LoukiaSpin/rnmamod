# Robustness index

Calculates the robustness index, a novel index that quantifies the
overall divergence of the sensitivity analysis results from the primary
analysis results. The robustness index considers objective decision
rules to infer the presence or lack of robustness of the primary
analysis results when conducting a sensitivity analysis (Spineli et al.,
2021).

## Usage

``` r
robustness_index(sens, prediction = FALSE, threshold)
```

## Arguments

- sens:

  An object of S3 class
  [`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md)
  when sensitivity analysis refers to different scenarios about the
  average missingness parameter. See 'Value' in
  [`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md).
  For a **general** sensitivity analysis, insert a list of at least two
  objects of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
  or
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md)
  indicating different re-analyses: the first object in the list should
  refer to the primary analysis.

- prediction:

  Logical character on whether to consider the prediction (`TRUE`) or
  estimation of the summary treatment effects (`FALSE`). This is only
  relevant for a random-effects model and the default argument is
  `FALSE` (estimation).

- threshold:

  A number indicating the threshold of robustness, that is, the
  minimally allowed deviation between the primary analysis and
  re-analysis results. See 'Details' below.

## Value

`robustness_index` prints on the R console a message in green text on
the threshold of robustness determined by the user. Then, the function
returns the following list of elements:

- robust_index:

  A numeric scalar or vector on the robustness index values. In the case
  of a pairwise meta-analysis, `robust_index` is scalar as only one
  summary effect size is obtained. In the case of network meta-analysis,
  `robust_index` is a vector with length equal to the number of possible
  pairwise comparisons; one robustness index per possible comparison.

- robust:

  A character or character vector (of same length with `robust_index`)
  on whether the primary analysis results are *robust* or *frail* to the
  different re-analyses.

- kld:

  A vector or matrix on the Kullback-Leibler divergence (KLD) measure in
  the summary effect size from a subsequent re-analysis to the primary
  analysis. In the case of a pairwise meta-analysis, `kld` is a vector
  with length equal to the number of total analyses (one KLD value is
  obtained per analysis). The number of total analyses equals the square
  of the number of scenarios indicated in the argument `mean_scenarios`
  of
  [`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md),
  in the case of missing participant outcome data; otherwise, the length
  of the character vector in argument `sens`. In the case of network
  meta-analysis, `robust_index` is a matrix with number of rows equal to
  the number of total analyses and number of columns equal to the number
  of possible pairwise comparisons; one KLD value per analysis and
  possible comparison.

- threshold:

  The threshold used to be inherited by the
  [`heatmap_robustness`](https://loukiaspin.github.io/rnmamod/reference/heatmap_robustness.md)
  function.

- scenarios:

  The scenarios considered to be inherited by the
  [`kld_barplot`](https://loukiaspin.github.io/rnmamod/reference/kld_barplot.md)
  function.

## Details

Thresholds of robustness have been proposed only for the odds ratio and
standardised mean difference (Spineli et al., 2021). The user may
consider the values 0.28 and 0.17 in the argument `threshold` for the
odds ratio and standardised mean difference effect measures (the default
values), respectively, or consider other plausible values. When the
argument `threshold` has not been defined, `robustness_index` considers
the default values 0.28 and 0.17 as threshold for robustness for binary
and continuous outcome, respectively, regardless of the effect measure
(the default thresholds may not be proper choices for other effect
measures; hence, use these threshold with great caution in this case).
Spineli et al. (2021) offers a discussion on specifying the `threshold`
of robustness.

In the case of binary outcome, `robustness_index` considers the results
in the odds ratio scale to calculate the robustness index. This is
because, the odds ratio is used as the 'best-case' effect measure in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
Then, relative risk, and risk difference are functions of the odds ratio
and the selected baseline risk (See 'Details' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

In the case of missing participant outcome data, the primary analysis is
considered to be the middle of the numbers in the argument
`mean_scenarios` of
[`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md)
(see 'Arguments' and 'Details' in
[`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md)).
Furhermore, `robustness_index` can be used in that context only when
missing participant outcome data have been extracted for at least one
trial. Otherwise, the execution of the function will be stopped and an
error message will be printed in the R console.

In the case of a general sensitivity analysis, the compared models
should refer to the same effect measure and the same meta-analysis model
(i.e., fixed-effect or random-effects).

In `robust`, the value `"robust"` appears when `robust_index` is less
than `threshold`; otherwise, the value `"frail"` appears.

## References

Kullback S, Leibler RA. On information and sufficiency. *Ann Math Stat*
1951;**22**(1):79–86. doi: 10.1214/aoms/1177729694

Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness
of primary analysis results: A case study on missing outcome data in
pairwise and network meta-analysis. *Res Synth Methods*
2021;**12**(4):475–90. doi: 10.1002/jrsm.1478

## See also

[`heatmap_robustness`](https://loukiaspin.github.io/rnmamod/reference/heatmap_robustness.md),
[`kld_barplot`](https://loukiaspin.github.io/rnmamod/reference/kld_barplot.md),
[`kld_measure`](https://loukiaspin.github.io/rnmamod/reference/kld_measure.md),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md)

## Author

Loukia M. Spineli

## Examples

``` r
data("nma.baker2009")

# Read results from 'run_sensitivity' (using the default arguments)
res_sens <- readRDS(system.file('extdata/res_sens_baker.rds',
                    package = 'rnmamod'))

# Calculate the robustness index
robustness_index(sens = res_sens,
                 threshold = 0.28)
#> The value 0.28 was assigned as 'threshold' for odds ratio.
#> $robust_index
#>  [1] 0.45825757 0.40000000 0.45387223 0.51478151 0.48062459 0.31780497
#>  [7] 0.34785054 0.08366600 0.11401754 0.08366600 0.10488088 0.24698178
#> [13] 0.26457513 0.07745967 0.06324555 0.07745967 0.18973666 0.20000000
#> [19] 0.12247449 0.06324555 0.18165902 0.20000000 0.04472136 0.18973666
#> [25] 0.21213203 0.11832160 0.16124515 0.05477226
#> 
#> $robust
#>  [1] "frail"  "frail"  "frail"  "frail"  "frail"  "frail"  "frail"  "robust"
#>  [9] "robust" "robust" "robust" "robust" "robust" "robust" "robust" "robust"
#> [17] "robust" "robust" "robust" "robust" "robust" "robust" "robust" "robust"
#> [25] "robust" "robust" "robust" "robust"
#> 
#> $kld
#>             [,1]       [,2]        [,3]        [,4]        [,5]        [,6]
#>  [1,] 0.03143496 0.02796963 0.055058416 0.125452928 0.196018804 0.013912157
#>  [2,] 0.02543702 0.01920250 0.043161377 0.091735629 0.146776953 0.019573044
#>  [3,] 0.03405799 0.03059699 0.056604557 0.116167819 0.170723488 0.016309098
#>  [4,] 0.02951539 0.04548642 0.093447495 0.146961033 0.214858236 0.009077864
#>  [5,] 0.01113434 0.03216688 0.109525620 0.166811784 0.247708559 0.003877332
#>  [6,] 0.01951722 0.01086042 0.020786834 0.091750311 0.165670683 0.029723758
#>  [7,] 0.04998651 0.02082673 0.021052758 0.076821621 0.147532529 0.046796072
#>  [8,] 0.04989347 0.03309033 0.019555306 0.027614613 0.033629565 0.005643233
#>  [9,] 0.04486888 0.03923226 0.033844111 0.041112954 0.052304306 0.014171171
#> [10,] 0.02902816 0.02737892 0.021368567 0.029429527 0.035216212 0.014288434
#> [11,] 0.01450164 0.01118354 0.010032892 0.014450268 0.019383796 0.004623395
#> [12,] 0.02734489 0.02961612 0.033151278 0.061701609 0.086408868 0.008612444
#> [13,] 0.03628689 0.03747297 0.047533123 0.080519739 0.109411182 0.012209439
#> [14,] 0.03472138 0.02893725 0.019598727 0.024299059 0.027933502 0.015519021
#> [15,] 0.03039907 0.03019554 0.018994336 0.015199623 0.018396590 0.019615407
#> [16,] 0.01750883 0.01454469 0.010723813 0.008458197 0.005050863 0.009345788
#> [17,] 0.01921020 0.02019464 0.024987792 0.037464730 0.051152380 0.011409928
#> [18,] 0.02358055 0.02438248 0.031703395 0.052986855 0.067429896 0.012530294
#> [19,] 0.05688734 0.06522904 0.046346822 0.040560882 0.037222429 0.027758787
#> [20,] 0.01958700 0.02872948 0.021625476 0.013624033 0.011923788 0.006174286
#> [21,] 0.04855291 0.04361876 0.030655428 0.035677143 0.031113092 0.019055721
#> [22,] 0.05115298 0.04878439 0.041096263 0.053712343 0.048622138 0.022627111
#> [23,] 0.01142135 0.01065841 0.008975878 0.004187517 0.003444494 0.005165492
#> [24,] 0.03599252 0.04075226 0.048066877 0.042427435 0.048430174 0.016616661
#> [25,] 0.05016524 0.06033360 0.062757451 0.061811728 0.070595835 0.022314447
#> [26,] 0.02050571 0.02854462 0.053504773 0.045924247 0.056261386 0.006777018
#> [27,] 0.03272550 0.04509633 0.069504861 0.071161202 0.084445355 0.006854829
#> [28,] 0.02820476 0.02854098 0.013259454 0.011535247 0.011121530 0.008623787
#>              [,7]         [,8]        [,9]       [,10]        [,11]
#>  [1,] 0.012635907 0.0345111366 0.076040440 0.110159742 0.0325710603
#>  [2,] 0.007552997 0.0235592374 0.057539137 0.105413277 0.0422239720
#>  [3,] 0.008545467 0.0381733605 0.084574047 0.118771536 0.0324497037
#>  [4,] 0.015190431 0.0246016818 0.091345984 0.136589695 0.0265719807
#>  [5,] 0.002163974 0.0308755774 0.100226627 0.150372515 0.0351852594
#>  [6,] 0.010629545 0.0094664492 0.057500734 0.105562756 0.0454659057
#>  [7,] 0.018758751 0.0088146221 0.060973183 0.119920601 0.0582538415
#>  [8,] 0.009707479 0.0070248754 0.009284702 0.006343006 0.0007590439
#>  [9,] 0.013152956 0.0199737870 0.015926428 0.021112447 0.0045872158
#> [10,] 0.012995884 0.0095565514 0.008615908 0.010552018 0.0046558621
#> [11,] 0.004736495 0.0063616130 0.006540508 0.008273813 0.0052234803
#> [12,] 0.013727905 0.0237629411 0.035964385 0.044094640 0.0106219700
#> [13,] 0.017485710 0.0289618004 0.040229347 0.048606665 0.0110262681
#> [14,] 0.015270182 0.0136434794 0.009343775 0.014808220 0.0072103373
#> [15,] 0.015447861 0.0052643668 0.006463093 0.007423185 0.0062578747
#> [16,] 0.004769091 0.0031395312 0.002539431 0.004803469 0.0073926899
#> [17,] 0.008054395 0.0152864520 0.023148034 0.037971823 0.0151285300
#> [18,] 0.008808176 0.0171088792 0.026238435 0.041346593 0.0156539829
#> [19,] 0.022516904 0.0136286036 0.011412669 0.013909778 0.0002166782
#> [20,] 0.004164008 0.0061287770 0.002216074 0.002699664 0.0002625998
#> [21,] 0.015660875 0.0253835007 0.022765288 0.022595283 0.0025400748
#> [22,] 0.016338439 0.0296181091 0.026026844 0.024768890 0.0032085867
#> [23,] 0.004946949 0.0009981894 0.003028267 0.001084037 0.0002915975
#> [24,] 0.016824835 0.0111546755 0.026618707 0.027447035 0.0017907530
#> [25,] 0.025559022 0.0150603987 0.030131129 0.032665880 0.0019734863
#> [26,] 0.006912700 0.0118231997 0.025919989 0.031108128 0.0024523076
#> [27,] 0.008316147 0.0163663580 0.031320234 0.036711683 0.0029309439
#> [28,] 0.010398057 0.0062078390 0.006016504 0.002358181 0.0009362008
#>              [,12] [,13]        [,14]        [,15]       [,16]       [,17]
#>  [1,] 1.040630e-02     0 1.164902e-02 4.204829e-02 0.134543244 0.095648309
#>  [2,] 1.276933e-02     0 1.089570e-02 3.910933e-02 0.122007318 0.089693065
#>  [3,] 6.267439e-03     0 1.678665e-02 5.603486e-02 0.121278648 0.075088547
#>  [4,] 7.006353e-03     0 1.125546e-02 2.877036e-02 0.124126947 0.101362316
#>  [5,] 1.093077e-02     0 1.200075e-02 6.427553e-02 0.130560278 0.073951720
#>  [6,] 2.032658e-02     0 1.877611e-02 7.476512e-02 0.094322149 0.044053262
#>  [7,] 2.190801e-02     0 1.982677e-02 8.307096e-02 0.107225101 0.052422580
#>  [8,] 1.623253e-04     0 4.557987e-04 1.463448e-03 0.007773402 0.006289221
#>  [9,] 3.362800e-03     0 6.533559e-04 2.881607e-03 0.020817558 0.017962139
#> [10,] 1.957335e-03     0 1.572967e-03 7.566119e-03 0.018121486 0.010187176
#> [11,] 1.383241e-03     0 2.215096e-03 3.495025e-03 0.029169640 0.025878108
#> [12,] 2.114031e-03     0 3.472683e-03 1.057402e-02 0.066073893 0.056075063
#> [13,] 3.254108e-03     0 3.538819e-03 1.091195e-02 0.066468155 0.055432680
#> [14,] 3.137903e-03     0 3.219365e-04 1.612032e-03 0.012596538 0.012370418
#> [15,] 2.377907e-03     0 8.200578e-04 5.435890e-03 0.009833125 0.006057058
#> [16,] 2.136754e-03     0 1.448610e-03 1.783064e-03 0.018598750 0.018648103
#> [17,] 2.873479e-03     0 2.586324e-03 8.376764e-03 0.052360163 0.047794559
#> [18,] 3.855192e-03     0 2.686450e-03 8.568660e-03 0.053571928 0.047993208
#> [19,] 2.045975e-03     0 5.051936e-04 2.105341e-03 0.002608169 0.004347026
#> [20,] 9.613066e-04     0 5.980695e-04 9.195487e-07 0.001499387 0.001222761
#> [21,] 1.130349e-03     0 2.172512e-03 4.312372e-03 0.027366544 0.022231718
#> [22,] 5.214653e-04     0 2.357629e-03 4.947864e-03 0.030091849 0.023513960
#> [23,] 2.024377e-04     0 8.139796e-05 2.278073e-03 0.001665226 0.005451690
#> [24,] 3.433409e-05     0 8.920290e-04 1.486432e-03 0.030988707 0.039076011
#> [25,] 2.852599e-04     0 9.069232e-04 1.174206e-03 0.031875678 0.038731295
#> [26,] 9.183757e-05     0 3.955442e-04 5.175883e-03 0.021081630 0.015812170
#> [27,] 2.845253e-04     0 5.291961e-04 6.232539e-03 0.023977626 0.017356998
#> [28,] 1.267331e-04     0 5.834451e-06 7.911580e-04 0.003982528 0.003640608
#>             [,18]        [,19]        [,20]       [,21]       [,22]       [,23]
#>  [1,] 0.042487843 0.0106463822 0.0049171837 0.223608521 0.192673633 0.100042059
#>  [2,] 0.038347943 0.0074420702 0.0031213380 0.221013974 0.164464609 0.080496336
#>  [3,] 0.042181108 0.0135991614 0.0060590213 0.237049354 0.176502923 0.119589142
#>  [4,] 0.049122862 0.0131807149 0.0046040366 0.264519469 0.199087503 0.118089411
#>  [5,] 0.029089442 0.0008059454 0.0070371123 0.205479242 0.143685400 0.050963673
#>  [6,] 0.010629566 0.0074674204 0.0249425848 0.143297831 0.080884399 0.027226127
#>  [7,] 0.015764281 0.0047428870 0.0243847732 0.176391820 0.099628581 0.031745531
#>  [8,] 0.008251357 0.0096907986 0.0077012507 0.010866329 0.013563288 0.011593977
#>  [9,] 0.007597600 0.0047053070 0.0036050973 0.024580727 0.027938245 0.012242492
#> [10,] 0.005336946 0.0044078932 0.0041772919 0.023145986 0.022593803 0.010828709
#> [11,] 0.012927402 0.0065187566 0.0043482660 0.053823714 0.054666538 0.038184489
#> [12,] 0.031573307 0.0157961102 0.0092958150 0.113079651 0.112847393 0.070583758
#> [13,] 0.029059236 0.0136730496 0.0078265796 0.109946905 0.111954549 0.069502717
#> [14,] 0.004377798 0.0058657778 0.0040183263 0.017991351 0.013758014 0.007134809
#> [15,] 0.003161499 0.0028137835 0.0031520756 0.014848275 0.008329188 0.004264062
#> [16,] 0.008723452 0.0037461567 0.0019333675 0.042459716 0.032061684 0.022740165
#> [17,] 0.026484981 0.0118861841 0.0054696232 0.101617982 0.084003351 0.051573564
#> [18,] 0.024655143 0.0109142310 0.0050677870 0.099074801 0.083634871 0.051537082
#> [19,] 0.002561615 0.0044773479 0.0033483441 0.006285256 0.006035002 0.003384348
#> [20,] 0.002674612 0.0049587453 0.0071346665 0.013073903 0.010624032 0.019157476
#> [21,] 0.023491309 0.0197374824 0.0134362025 0.074739245 0.065050674 0.064800681
#> [22,] 0.021975735 0.0203392339 0.0161521210 0.076337861 0.066970903 0.067684648
#> [23,] 0.003669885 0.0060804031 0.0083579692 0.010711292 0.010860714 0.015395395
#> [24,] 0.026903746 0.0223632534 0.0186015010 0.076165711 0.074754399 0.060983392
#> [25,] 0.023884413 0.0194777145 0.0163754149 0.072920624 0.072687844 0.061373953
#> [26,] 0.011260031 0.0045273984 0.0011704912 0.032311489 0.030611620 0.013912927
#> [27,] 0.010123069 0.0035442703 0.0009974872 0.031341715 0.031484138 0.015100289
#> [28,] 0.005745668 0.0073554835 0.0053779593 0.012360609 0.011709647 0.013174067
#>             [,24]       [,25]
#>  [1,] 0.046230251 0.020476941
#>  [2,] 0.034039561 0.016574904
#>  [3,] 0.063417681 0.046690905
#>  [4,] 0.053360895 0.034057050
#>  [5,] 0.008449433 0.001574758
#>  [6,] 0.009511592 0.018432498
#>  [7,] 0.011516544 0.019748169
#>  [8,] 0.011697306 0.010519283
#>  [9,] 0.011427620 0.006319430
#> [10,] 0.006878600 0.005632575
#> [11,] 0.026566803 0.018438287
#> [12,] 0.044682077 0.031917902
#> [13,] 0.042379734 0.028876100
#> [14,] 0.008525952 0.008171267
#> [15,] 0.003196416 0.005644399
#> [16,] 0.016278182 0.013737923
#> [17,] 0.031811625 0.026189122
#> [18,] 0.030723705 0.024288117
#> [19,] 0.005552698 0.003038684
#> [20,] 0.021835370 0.031279222
#> [21,] 0.054745546 0.064624548
#> [22,] 0.055360457 0.064022863
#> [23,] 0.016557831 0.024647105
#> [24,] 0.046004530 0.055025719
#> [25,] 0.044052337 0.051021511
#> [26,] 0.005440407 0.003660743
#> [27,] 0.005066223 0.003251848
#> [28,] 0.010746082 0.012487789
#> 
#> $measure
#> [1] "OR"
#> 
#> $threshold
#> [1] 0.28
#> 
#> $scenarios
#> [1] -1.098612289 -0.693147181 -0.000100005  0.693147181  1.098612289
#> 
#> attr(,"class")
#> [1] "robustness_index"
```
