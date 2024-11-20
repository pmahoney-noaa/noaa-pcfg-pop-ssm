# Projecting PCFG Gray Whale Abundance


<div class="cell-output-display">

| Model             | chain | num_divergent | num_max_treedepth | ebfmi |
|:------------------|:-----:|:-------------:|:-----------------:|:-----:|
| Base              |   1   |       0       |         0         | 0.47  |
| Base              |   2   |       0       |         0         | 0.51  |
| Base              |   3   |       0       |         0         | 0.33  |
| AR1v1             |   1   |       0       |         0         | 0.35  |
| AR1v1             |   2   |       0       |         0         | 0.33  |
| AR1v1             |   3   |       0       |         0         | 0.31  |
| AR1v2             |   1   |       0       |         0         | 0.54  |
| AR1v2             |   2   |       0       |         0         | 0.53  |
| AR1v2             |   3   |       0       |         0         | 0.51  |
| ENP Calves        |   1   |       0       |         0         | 0.65  |
| ENP Calves        |   2   |       0       |         0         | 0.59  |
| ENP Calves        |   3   |       0       |         0         | 0.68  |
| Calves/Strandings |   1   |       0       |         0         | 0.33  |
| Calves/Strandings |   2   |       0       |         0         | 0.38  |
| Calves/Strandings |   3   |       0       |         0         | 0.46  |
| Calves only       |   1   |       0       |         0         | 0.48  |
| Calves only       |   2   |       0       |         0         | 0.50  |
| Calves only       |   3   |       0       |         0         | 0.52  |
| Strandings only   |   1   |       1       |         0         | 0.35  |
| Strandings only   |   2   |       0       |         0         | 0.35  |
| Strandings only   |   3   |       0       |         0         | 0.38  |

</div>

## Leave-one-out Cross Validation (LOO)

<div class="cell-output-display">

| Model             | elpd_loo | p_loo | looic  | deltaLooic |
|:------------------|:--------:|:-----:|:------:|:----------:|
| AR1v1             |  25.55   | 2.37  | -51.11 |    0.00    |
| Calves only       |  24.78   | 2.84  | -49.56 |    1.54    |
| Strandings only   |  24.44   | 2.88  | -48.89 |    2.22    |
| Calves/Strandings |  24.32   | 3.09  | -48.65 |    2.46    |
| Base              |  24.22   | 3.11  | -48.44 |    2.66    |
| ENP Calves        |  23.86   | 3.37  | -47.71 |    3.39    |
| AR1v2             |  23.34   | 3.73  | -46.69 |    4.42    |

</div>

## Retrospective analysis

<div class="cell-output-display">

| Model             | Mean RSS | Mean Percentile (N) | Prop Below Threshold (N) | Prop Below Threshold (Nmin) | Number of Closures |
|:------------------|:--------:|:-------------------:|:------------------------:|:---------------------------:|:------------------:|
| AR1v1             | 1333782  |        0.519        |          0.123           |            0.017            |         0          |
| Strandings only   | 1377657  |        0.518        |          0.122           |            0.016            |         0          |
| Base              | 1413790  |        0.513        |          0.110           |            0.013            |         0          |
| Calves only       | 1435932  |        0.522        |          0.138           |            0.017            |         0          |
| Calves/Strandings | 1469416  |        0.521        |          0.145           |            0.021            |         0          |
| AR1v2             | 1608977  |        0.514        |          0.121           |            0.017            |         0          |
| ENP Calves        | 1610132  |        0.527        |          0.126           |            0.017            |         0          |

</div>

<div class="cell-output-display">

| Model             | Mean RSS | Mean Percentile (N) | Prop Below Threshold (N) | Prop Below Threshold (Nmin) | Number of Closures |
|:------------------|:--------:|:-------------------:|:------------------------:|:---------------------------:|:------------------:|
| Strandings only   | 2479358  |        0.509        |          0.154           |            0.036            |         0          |
| Base              | 2623142  |        0.504        |          0.147           |            0.030            |         0          |
| Calves only       | 2690216  |        0.509        |          0.185           |            0.042            |         0          |
| Calves/Strandings | 2761271  |        0.509        |          0.188           |            0.048            |         0          |
| AR1v2             | 2923377  |        0.509        |          0.163           |            0.039            |         0          |
| ENP Calves        | 3037646  |        0.512        |          0.165           |            0.040            |         0          |
| AR1v1             | 4039921  |        0.500        |          0.203           |            0.071            |         0          |

</div>

![](README_files/figure-commonmark/retroFig-1.png)

## Model-specific trends and projections

Note, the number of PCFG calves used in projected years (models
Calves/Strandings and Calves only) are not accurate and likely represent
underestimates of reality.

![](README_files/figure-commonmark/trendFig-1.png)

## Model-specific predictions for Y<sub>final</sub> + 2

![](README_files/figure-commonmark/fig-final-proj-1.png)
