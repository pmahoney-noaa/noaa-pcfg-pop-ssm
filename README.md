# Projecting PCFG Gray Whale Abundance


# Description of model set

The abundance of gray whales within the Pacific Coast Feeding Group
(PCFG) - a subset of the broader Eastern North Pacific (ENP)
population - is estimated within a Jolly-Seber modeling framework using
data derived from yearly mark-resight surveys (Harris et al. 2024). Due
to the time required to match individuals from many hundreds of
sightings each year to an extensive photographic catalog, the most
recent abundance estimates often lag behind the current calendar year by
one or more years. However, management decisions are contingent upon
abundance estimates for PCFG gray whales through the current calendar
year or into the next. To meet this need, we implemented a Bayesian
state-space model (SSM) to predict abundance beyond the last year for
which PCFG abundance was estimated. The SSM models abundance in
log-space with Gaussian errors and can be defined as the following:

*Observation process:*

![logN\_{est} \sim normal(logN, \sigma\_{obs})](https://latex.codecogs.com/svg.latex?logN_%7Best%7D%20%5Csim%20normal%28logN%2C%20%5Csigma_%7Bobs%7D%29 "logN_{est} \sim normal(logN, \sigma_{obs})")

*State process:*

![logN\_{est\_{t+1}} = logN\_{est\_{t}} + ln(\lambda\_{t})](https://latex.codecogs.com/svg.latex?logN_%7Best_%7Bt%2B1%7D%7D%20%3D%20logN_%7Best_%7Bt%7D%7D%20%2B%20ln%28%5Clambda_%7Bt%7D%29 "logN_{est_{t+1}} = logN_{est_{t}} + ln(\lambda_{t})")

![ln(\lambda\_{t}) \sim normal(\mu, \sigma)](https://latex.codecogs.com/svg.latex?ln%28%5Clambda_%7Bt%7D%29%20%5Csim%20normal%28%5Cmu%2C%20%5Csigma%29 "ln(\lambda_{t}) \sim normal(\mu, \sigma)")

Where ![logN](https://latex.codecogs.com/svg.latex?logN "logN") and
![\sigma\_{obs}](https://latex.codecogs.com/svg.latex?%5Csigma_%7Bobs%7D "\sigma_{obs}")
are the estimated (log) population abundance and standard error for each
year with observation data, respectively. Models vary by the
specification of
![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda") and
the factors hypothesized to be correlated with population growth, at
least given current knowledge and available data within this system. At
present, seven models were considered.

1.  Base:
    ![\mu = \mu\_{ln(\lambda)}; \sigma = \sigma\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cmu%20%3D%20%5Cmu_%7Bln%28%5Clambda%29%7D%3B%20%5Csigma%20%3D%20%5Csigma_%7Bln%28%5Clambda%29%7D "\mu = \mu_{ln(\lambda)}; \sigma = \sigma_{ln(\lambda)}")
2.  AR1v1:
    ![\mu = ln(\lambda\_{t-1}); \sigma = \sigma\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cmu%20%3D%20ln%28%5Clambda_%7Bt-1%7D%29%3B%20%5Csigma%20%3D%20%5Csigma_%7Bln%28%5Clambda%29%7D "\mu = ln(\lambda_{t-1}); \sigma = \sigma_{ln(\lambda)}")
3.  AR2v2:
    ![\mu = \mu\_{ln(\lambda)} + ln(\lambda\_{t-1}) \* \beta; \sigma = \sigma\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cmu%20%3D%20%5Cmu_%7Bln%28%5Clambda%29%7D%20%2B%20ln%28%5Clambda_%7Bt-1%7D%29%20%2A%20%5Cbeta%3B%20%5Csigma%20%3D%20%5Csigma_%7Bln%28%5Clambda%29%7D "\mu = \mu_{ln(\lambda)} + ln(\lambda_{t-1}) * \beta; \sigma = \sigma_{ln(\lambda)}")
4.  ENP calves:
    ![\mu = \mu\_{ln(\lambda)} + \beta\_{calves} \* C\_{t}; \sigma = \sigma\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cmu%20%3D%20%5Cmu_%7Bln%28%5Clambda%29%7D%20%2B%20%5Cbeta_%7Bcalves%7D%20%2A%20C_%7Bt%7D%3B%20%5Csigma%20%3D%20%5Csigma_%7Bln%28%5Clambda%29%7D "\mu = \mu_{ln(\lambda)} + \beta_{calves} * C_{t}; \sigma = \sigma_{ln(\lambda)}")
5.  (PCFG) Calves only:
    ![\mu = \mu\_{ln(\lambda)} + \beta\_{calves} \* C\_{t}; \sigma = \sigma\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cmu%20%3D%20%5Cmu_%7Bln%28%5Clambda%29%7D%20%2B%20%5Cbeta_%7Bcalves%7D%20%2A%20C_%7Bt%7D%3B%20%5Csigma%20%3D%20%5Csigma_%7Bln%28%5Clambda%29%7D "\mu = \mu_{ln(\lambda)} + \beta_{calves} * C_{t}; \sigma = \sigma_{ln(\lambda)}")
6.  (ENP) Strandings only:
    ![\mu = \mu\_{ln(\lambda)} + \beta\_{strand} \* S\_{t}; \sigma = \sigma\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cmu%20%3D%20%5Cmu_%7Bln%28%5Clambda%29%7D%20%2B%20%5Cbeta_%7Bstrand%7D%20%2A%20S_%7Bt%7D%3B%20%5Csigma%20%3D%20%5Csigma_%7Bln%28%5Clambda%29%7D "\mu = \mu_{ln(\lambda)} + \beta_{strand} * S_{t}; \sigma = \sigma_{ln(\lambda)}")
7.  (PCFG) Calves + Strandings:
    ![\mu = \mu\_{ln(\lambda)} + \beta\_{calves} \* C\_{t} + \beta\_{strand} \* S\_{t}; \sigma = \sigma\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cmu%20%3D%20%5Cmu_%7Bln%28%5Clambda%29%7D%20%2B%20%5Cbeta_%7Bcalves%7D%20%2A%20C_%7Bt%7D%20%2B%20%5Cbeta_%7Bstrand%7D%20%2A%20S_%7Bt%7D%3B%20%5Csigma%20%3D%20%5Csigma_%7Bln%28%5Clambda%29%7D "\mu = \mu_{ln(\lambda)} + \beta_{calves} * C_{t} + \beta_{strand} * S_{t}; \sigma = \sigma_{ln(\lambda)}")

<!--#
&#10;1. Base: $ln(\lambda_{t}) \sim normal(\mu_{ln(\lambda)}, \sigma_{ln(\lambda)})$ 
2. AR1v1: $ln(\lambda_{t}) \sim normal(ln(\lambda_{t-1}), \sigma_{ln(\lambda_{t})})$ 
3. AR2v2: $ln(\lambda_{t}) \sim normal(\mu_{ln(\lambda_{t})} + ln(\lambda_{t-1}) * \beta, \sigma_{ln(\lambda_{t})})$ 
4. ENP calves: $ln(\lambda_{t}) \sim normal(\mu_{ln(\lambda_{t})} + C_{t} * \beta_{calves}, \sigma_{ln(\lambda_{t})})$ 
5. (PCFG) Calves only: $ln(\lambda_{t}) \sim normal(\mu_{ln(\lambda_{t})} + C_{t} * \beta_{calves}, \sigma_{ln(\lambda_{t})})$ 
6. (ENP) Strandings only: $ln(\lambda_{t}) \sim normal(\mu_{ln\lambda} + S_{t} * \beta_{strand}, \sigma_{ln\lambda})$ 
7. (PCFG) Calves + Strandings: $ln(\lambda_{t}) \sim normal(\mu_{ln(\lambda_{t})} + C_{t} * \beta_{calves} + S_{t} * \beta_{strand}, \sigma_{ln(\lambda_{t})})$
&#10;-->

![C\_{t}](https://latex.codecogs.com/svg.latex?C_%7Bt%7D "C_{t}")
corresponds the number of calves in year *t* (either ENP or PCFG);
![S\_{t}](https://latex.codecogs.com/svg.latex?S_%7Bt%7D "S_{t}")
corresponds to the number of ENP strandings in year *t*.

### Additional prior specifications

Where appropriate, the following weakly informative hyper-priors were
implemented for both mean and scale parameters in the
![ln(\lambda)](https://latex.codecogs.com/svg.latex?ln%28%5Clambda%29 "ln(\lambda)")
priors.

![\mu\_{ln(\lambda)} \sim normal(0, 1)](https://latex.codecogs.com/svg.latex?%5Cmu_%7Bln%28%5Clambda%29%7D%20%5Csim%20normal%280%2C%201%29 "\mu_{ln(\lambda)} \sim normal(0, 1)")

![\sigma\_{ln(\lambda)} \sim lognormal(0, 1)](https://latex.codecogs.com/svg.latex?%5Csigma_%7Bln%28%5Clambda%29%7D%20%5Csim%20lognormal%280%2C%201%29 "\sigma_{ln(\lambda)} \sim lognormal(0, 1)")

Importantly,
![\sigma\_{ln(\lambda\_{t})}](https://latex.codecogs.com/svg.latex?%5Csigma_%7Bln%28%5Clambda_%7Bt%7D%29%7D "\sigma_{ln(\lambda_{t})}")
is constrained to be positive. For models with beta coefficients (models
4 - 7), additional normal priors are specified for each coefficient, *i*
(e.g.,
![ln(\lambda\_{t-1})](https://latex.codecogs.com/svg.latex?ln%28%5Clambda_%7Bt-1%7D%29 "ln(\lambda_{t-1})"),
ENP calves, PCFG calves, or number of strandings). At present, no
hyper-priors were used for coefficients.

![\beta_i \sim normal(0, 1)](https://latex.codecogs.com/svg.latex?%5Cbeta_i%20%5Csim%20normal%280%2C%201%29 "\beta_i \sim normal(0, 1)")

# Model fitting diagnostics

The table below provides various chain-specific diagnostics by model
specification. Chains with few to no divergences (num_divergent) and
that do not hit maximum tree depth (num_max_treedepth) indicate adequate
sampling of posteriors and overall good chain performance. Similarly, an
estimated Bayesian fraction of missing information (ebfmi) greater than
0.3 suggests reliable inference can be drawn from posterior estimates.

<div class="cell-output-display">

| Model                        | chain | num_divergent | num_max_treedepth | ebfmi |
|:-----------------------------|:-----:|:-------------:|:-----------------:|:-----:|
| Base                         |   1   |       0       |         0         | 0.47  |
| Base                         |   2   |       0       |         0         | 0.51  |
| Base                         |   3   |       0       |         0         | 0.33  |
| AR1v1                        |   1   |       0       |         0         | 0.35  |
| AR1v1                        |   2   |       0       |         0         | 0.33  |
| AR1v1                        |   3   |       0       |         0         | 0.31  |
| AR1v2                        |   1   |       0       |         0         | 0.54  |
| AR1v2                        |   2   |       0       |         0         | 0.53  |
| AR1v2                        |   3   |       0       |         0         | 0.51  |
| ENP Calves                   |   1   |       0       |         0         | 0.65  |
| ENP Calves                   |   2   |       0       |         0         | 0.59  |
| ENP Calves                   |   3   |       0       |         0         | 0.68  |
| PCFG Calves only             |   1   |       0       |         0         | 0.48  |
| PCFG Calves only             |   2   |       0       |         0         | 0.50  |
| PCFG Calves only             |   3   |       0       |         0         | 0.52  |
| ENP Strandings only          |   1   |       1       |         0         | 0.35  |
| ENP Strandings only          |   2   |       0       |         0         | 0.35  |
| ENP Strandings only          |   3   |       0       |         0         | 0.38  |
| PCFG Calves + ENP Strandings |   1   |       0       |         0         | 0.33  |
| PCFG Calves + ENP Strandings |   2   |       0       |         0         | 0.38  |
| PCFG Calves + ENP Strandings |   3   |       0       |         0         | 0.46  |

</div>

# Model results

## Coefficient estimates

![\beta](https://latex.codecogs.com/svg.latex?%5Cbeta "\beta")
coefficient estimates for models where coefficients on
![\lambda](https://latex.codecogs.com/svg.latex?%5Clambda "\lambda")
were estimated.
![lo\\ci](https://latex.codecogs.com/svg.latex?lo%5C_ci "lo\_ci") and
![hi\\ci](https://latex.codecogs.com/svg.latex?hi%5C_ci "hi\_ci")
correspond to the lower and upper 95% credible intervals.
![Rhat](https://latex.codecogs.com/svg.latex?Rhat "Rhat") values less
than 1.01 indicate effective mixing within and across chains, producing
reliable and consistent estimates for the parameters.
![ess](https://latex.codecogs.com/svg.latex?ess "ess") is the effective
sample size and is an indicator of chain sampling efficiency. Here,
![ess](https://latex.codecogs.com/svg.latex?ess "ess") values greater
than 300 indicate sufficient sampling has been achieved.

| model                        |                                                      variable                                                       |  mean  | median |  sd   | lo_ci  | hi_ci | rhat  | ess_bulk |
|:-----------------------------|:-------------------------------------------------------------------------------------------------------------------:|:------:|:------:|:-----:|:------:|:-----:|:-----:|:--------:|
| AR1v2                        | ![\beta\_{ln(\lambda)}](https://latex.codecogs.com/svg.latex?%5Cbeta_%7Bln%28%5Clambda%29%7D "\beta_{ln(\lambda)}") | -0.167 | -0.170 | 0.382 | -0.784 | 0.463 | 1.004 |  758.2   |
| ENP Calves                   |           ![\beta\_{Calves}](https://latex.codecogs.com/svg.latex?%5Cbeta_%7BCalves%7D "\beta_{Calves}")            | -0.004 | -0.004 | 0.036 | -0.064 | 0.055 | 1.000 |  1447.6  |
| PCFG Calves only             |           ![\beta\_{Calves}](https://latex.codecogs.com/svg.latex?%5Cbeta_%7BCalves%7D "\beta_{Calves}")            | 0.020  | 0.019  | 0.019 | -0.010 | 0.051 | 1.001 |  2793.6  |
| ENP Strandings only          |     ![\beta\_{Strandings}](https://latex.codecogs.com/svg.latex?%5Cbeta_%7BStrandings%7D "\beta_{Strandings}")      | -0.002 | -0.002 | 0.001 | -0.004 | 0.001 | 1.000 |  2198.9  |
| PCFG Calves + ENP Strandings |           ![\beta\_{Calves}](https://latex.codecogs.com/svg.latex?%5Cbeta_%7BCalves%7D "\beta_{Calves}")            | 0.017  | 0.016  | 0.020 | -0.016 | 0.050 | 1.004 |  2694.3  |
| PCFG Calves + ENP Strandings |     ![\beta\_{Strandings}](https://latex.codecogs.com/svg.latex?%5Cbeta_%7BStrandings%7D "\beta_{Strandings}")      | -0.001 | -0.001 | 0.002 | -0.004 | 0.001 | 1.002 |  1998.7  |

## Leave-one-out Cross Validation (LOO)

Derived estimates of model fit using an information criterion based on
leave-one-out cross validation (looic). Models with the lowest
![looic](https://latex.codecogs.com/svg.latex?looic "looic") value are
considered the best fitting model given the time series of abundance
data used to fit the model.
![\Delta looic](https://latex.codecogs.com/svg.latex?%5CDelta%20looic "\Delta looic")
represents the difference in
![looic](https://latex.codecogs.com/svg.latex?looic "looic") relative to
the model with the lowest
![looic](https://latex.codecogs.com/svg.latex?looic "looic").

<div class="cell-output-display">

| Model                        | elpd_loo | elpd_loo_se | p_loo | p_loo_se | looic  | looic_se | deltaLooic |
|:-----------------------------|:--------:|:-----------:|:-----:|:--------:|:------:|:--------:|:----------:|
| AR1v1                        |  25.55   |    0.85     | 2.37  |   0.61   | -51.11 |   1.70   |    0.00    |
| PCFG Calves only             |  24.78   |    0.65     | 2.84  |   0.54   | -49.56 |   1.31   |    1.54    |
| ENP Strandings only          |  24.44   |    0.79     | 2.88  |   0.45   | -48.89 |   1.57   |    2.22    |
| PCFG Calves + ENP Strandings |  24.32   |    0.74     | 3.09  |   0.53   | -48.65 |   1.49   |    2.46    |
| Base                         |  24.22   |    0.71     | 3.11  |   0.65   | -48.44 |   1.42   |    2.66    |
| ENP Calves                   |  23.86   |    0.74     | 3.37  |   0.54   | -47.71 |   1.48   |    3.39    |
| AR1v2                        |  23.34   |    0.68     | 3.73  |   0.66   | -46.69 |   1.36   |    4.42    |

</div>

## Predictive accuracy with retrospection

Model predictions when projecting one year forward for all data years
from 2002 through 2021. RSS is the residual sum of squares.

<div class="cell-output-display">

| Year | Model                        | Abundance Estimate | Mean N  | Median N | Lower 95% CI (N) | Upper 95% CI (N) |  Nmin   | Closure | Percentile for Abundance Estimate | Prop Below Threshold (N) | Prop Below Threshold (Nmin) |
|:-----|:-----------------------------|:------------------:|:-------:|:--------:|:----------------:|:----------------:|:-------:|:-------:|:---------------------------------:|:------------------------:|:---------------------------:|
| 2003 | Base                         |       208.7        | 200.434 | 200.357  |     168.471      |     234.323      | 187.406 |  FALSE  |               0.716               |          0.295           |            0.033            |
| 2003 | AR1v1                        |       208.7        | 205.196 | 203.962  |     171.201      |     247.051      | 191.053 |  FALSE  |               0.618               |          0.218           |            0.024            |
| 2003 | AR1v2                        |       208.7        | 199.524 | 198.434  |     165.133      |     239.238      | 184.499 |  FALSE  |               0.711               |          0.348           |            0.056            |
| 2003 | ENP Calves                   |       208.7        | 199.364 | 198.992  |     166.550      |     235.846      | 185.357 |  FALSE  |               0.723               |          0.329           |            0.042            |
| 2003 | PCFG Calves only             |       208.7        | 199.689 | 199.053  |     166.971      |     235.516      | 186.530 |  FALSE  |               0.724               |          0.325           |            0.037            |
| 2003 | ENP Strandings only          |       208.7        | 202.108 | 201.602  |     171.813      |     233.348      | 189.585 |  FALSE  |               0.677               |          0.250           |            0.020            |
| 2003 | PCFG Calves + ENP Strandings |       208.7        | 200.982 | 200.559  |     166.355      |     236.541      | 187.160 |  FALSE  |               0.693               |          0.291           |            0.040            |
| 2004 | Base                         |       214.9        | 206.328 | 205.552  |     173.120      |     244.794      | 192.109 |  FALSE  |               0.714               |          0.198           |            0.020            |
| 2004 | AR1v1                        |       214.9        | 209.020 | 207.079  |     175.212      |     251.902      | 194.803 |  FALSE  |               0.674               |          0.153           |            0.014            |
| 2004 | AR1v2                        |       214.9        | 210.755 | 209.385  |     173.956      |     254.552      | 194.252 |  FALSE  |               0.609               |          0.165           |            0.015            |
| 2004 | ENP Calves                   |       214.9        | 205.906 | 205.293  |     170.098      |     247.592      | 191.471 |  FALSE  |               0.716               |          0.211           |            0.028            |
| 2004 | PCFG Calves only             |       214.9        | 205.940 | 204.995  |     171.610      |     245.408      | 191.489 |  FALSE  |               0.711               |          0.207           |            0.023            |
| 2004 | ENP Strandings only          |       214.9        | 208.994 | 207.993  |     175.547      |     248.298      | 195.385 |  FALSE  |               0.671               |          0.138           |            0.014            |
| 2004 | PCFG Calves + ENP Strandings |       214.9        | 208.598 | 207.811  |     173.763      |     248.635      | 193.948 |  FALSE  |               0.666               |          0.172           |            0.018            |
| 2005 | Base                         |       220.0        | 210.088 | 208.484  |     177.052      |     249.344      | 195.903 |  FALSE  |               0.738               |          0.136           |            0.011            |
| 2005 | AR1v1                        |       220.0        | 209.253 | 207.997  |     172.831      |     253.563      | 194.643 |  FALSE  |               0.752               |          0.159           |            0.023            |
| 2005 | AR1v2                        |       220.0        | 210.969 | 210.328  |     172.587      |     253.275      | 195.209 |  FALSE  |               0.705               |          0.154           |            0.021            |
| 2005 | ENP Calves                   |       220.0        | 209.801 | 208.536  |     174.083      |     253.238      | 194.465 |  FALSE  |               0.737               |          0.162           |            0.018            |
| 2005 | PCFG Calves only             |       220.0        | 208.843 | 207.715  |     176.627      |     248.323      | 193.908 |  FALSE  |               0.756               |          0.166           |            0.014            |
| 2005 | ENP Strandings only          |       220.0        | 211.320 | 210.053  |     178.958      |     251.971      | 197.370 |  FALSE  |               0.725               |          0.124           |            0.011            |
| 2005 | PCFG Calves + ENP Strandings |       220.0        | 210.338 | 209.358  |     175.127      |     248.637      | 195.920 |  FALSE  |               0.728               |          0.132           |            0.014            |
| 2006 | Base                         |       199.2        | 209.609 | 208.795  |     174.333      |     251.336      | 194.743 |  FALSE  |               0.284               |          0.159           |            0.017            |
| 2006 | AR1v1                        |       199.2        | 207.736 | 206.882  |     171.124      |     250.977      | 192.671 |  FALSE  |               0.324               |          0.190           |            0.024            |
| 2006 | AR1v2                        |       199.2        | 210.908 | 210.403  |     170.587      |     254.448      | 194.348 |  FALSE  |               0.278               |          0.169           |            0.027            |
| 2006 | ENP Calves                   |       199.2        | 209.505 | 207.880  |     173.260      |     254.424      | 192.973 |  FALSE  |               0.318               |          0.187           |            0.018            |
| 2006 | PCFG Calves only             |       199.2        | 206.707 | 205.440  |     170.377      |     250.324      | 191.047 |  FALSE  |               0.353               |          0.211           |            0.028            |
| 2006 | ENP Strandings only          |       199.2        | 210.706 | 209.427  |     176.911      |     250.371      | 196.866 |  FALSE  |               0.248               |          0.130           |            0.012            |
| 2006 | PCFG Calves + ENP Strandings |       199.2        | 207.945 | 207.461  |     171.825      |     249.450      | 192.371 |  FALSE  |               0.321               |          0.191           |            0.024            |
| 2007 | Base                         |       197.7        | 207.294 | 207.061  |     169.376      |     247.231      | 192.952 |  FALSE  |               0.286               |          0.185           |            0.029            |
| 2007 | AR1v1                        |       197.7        | 207.291 | 206.788  |     171.536      |     247.495      | 192.214 |  FALSE  |               0.299               |          0.197           |            0.024            |
| 2007 | AR1v2                        |       197.7        | 208.388 | 208.036  |     168.666      |     251.035      | 191.912 |  FALSE  |               0.295               |          0.201           |            0.032            |
| 2007 | ENP Calves                   |       197.7        | 206.331 | 205.450  |     168.060      |     249.205      | 189.886 |  FALSE  |               0.331               |          0.228           |            0.035            |
| 2007 | PCFG Calves only             |       197.7        | 204.025 | 203.150  |     168.098      |     243.595      | 188.751 |  FALSE  |               0.372               |          0.255           |            0.036            |
| 2007 | ENP Strandings only          |       197.7        | 206.815 | 206.483  |     173.125      |     243.940      | 192.715 |  FALSE  |               0.292               |          0.189           |            0.020            |
| 2007 | PCFG Calves + ENP Strandings |       197.7        | 203.911 | 203.811  |     168.533      |     241.900      | 188.876 |  FALSE  |               0.364               |          0.247           |            0.034            |
| 2008 | Base                         |       211.8        | 208.051 | 207.802  |     169.935      |     248.740      | 193.011 |  FALSE  |               0.592               |          0.188           |            0.028            |
| 2008 | AR1v1                        |       211.8        | 207.829 | 207.236  |     171.502      |     247.874      | 193.167 |  FALSE  |               0.603               |          0.179           |            0.023            |
| 2008 | AR1v2                        |       211.8        | 209.231 | 208.222  |     170.397      |     250.617      | 192.880 |  FALSE  |               0.577               |          0.187           |            0.027            |
| 2008 | ENP Calves                   |       211.8        | 206.497 | 206.101  |     169.204      |     248.071      | 190.883 |  FALSE  |               0.624               |          0.216           |            0.031            |
| 2008 | PCFG Calves only             |       211.8        | 203.492 | 202.931  |     167.536      |     242.764      | 188.268 |  FALSE  |               0.689               |          0.268           |            0.042            |
| 2008 | ENP Strandings only          |       211.8        | 207.630 | 207.138  |     172.373      |     246.336      | 193.697 |  FALSE  |               0.610               |          0.175           |            0.022            |
| 2008 | PCFG Calves + ENP Strandings |       211.8        | 204.026 | 203.744  |     166.583      |     245.465      | 188.689 |  FALSE  |               0.677               |          0.249           |            0.039            |
| 2009 | Base                         |       213.8        | 210.511 | 210.118  |     174.370      |     251.949      | 195.516 |  FALSE  |               0.591               |          0.153           |            0.018            |
| 2009 | AR1v1                        |       213.8        | 209.714 | 209.408  |     173.074      |     248.187      | 195.272 |  FALSE  |               0.610               |          0.147           |            0.019            |
| 2009 | AR1v2                        |       213.8        | 210.977 | 210.427  |     170.917      |     255.120      | 194.250 |  FALSE  |               0.573               |          0.164           |            0.025            |
| 2009 | ENP Calves                   |       213.8        | 210.043 | 209.349  |     173.761      |     250.091      | 194.800 |  FALSE  |               0.603               |          0.157           |            0.021            |
| 2009 | PCFG Calves only             |       213.8        | 205.076 | 204.604  |     168.502      |     244.694      | 190.115 |  FALSE  |               0.701               |          0.237           |            0.032            |
| 2009 | ENP Strandings only          |       213.8        | 213.159 | 213.087  |     176.517      |     251.278      | 198.524 |  FALSE  |               0.519               |          0.118           |            0.016            |
| 2009 | PCFG Calves + ENP Strandings |       213.8        | 208.566 | 208.096  |     171.852      |     251.102      | 192.543 |  FALSE  |               0.623               |          0.193           |            0.021            |
| 2010 | Base                         |       205.3        | 211.895 | 211.057  |     175.200      |     252.207      | 196.591 |  FALSE  |               0.365               |          0.145           |            0.017            |
| 2010 | AR1v1                        |       205.3        | 211.698 | 212.044  |     172.035      |     251.151      | 196.503 |  FALSE  |               0.351               |          0.143           |            0.022            |
| 2010 | AR1v2                        |       205.3        | 212.490 | 211.587  |     175.790      |     253.274      | 196.017 |  FALSE  |               0.370               |          0.143           |            0.013            |
| 2010 | ENP Calves                   |       205.3        | 212.031 | 211.135  |     173.433      |     255.855      | 195.283 |  FALSE  |               0.373               |          0.152           |            0.020            |
| 2010 | PCFG Calves only             |       205.3        | 207.193 | 206.200  |     170.224      |     249.205      | 191.321 |  FALSE  |               0.479               |          0.209           |            0.028            |
| 2010 | ENP Strandings only          |       205.3        | 214.347 | 214.421  |     178.736      |     252.684      | 199.486 |  FALSE  |               0.300               |          0.105           |            0.010            |
| 2010 | PCFG Calves + ENP Strandings |       205.3        | 210.097 | 209.365  |     172.355      |     252.401      | 194.304 |  FALSE  |               0.410               |          0.165           |            0.022            |
| 2011 | Base                         |       209.9        | 213.210 | 213.413  |     176.704      |     251.150      | 198.470 |  FALSE  |               0.422               |          0.123           |            0.013            |
| 2011 | AR1v1                        |       209.9        | 217.820 | 217.834  |     180.712      |     255.829      | 203.572 |  FALSE  |               0.321               |          0.082           |            0.011            |
| 2011 | AR1v2                        |       209.9        | 214.658 | 214.540  |     174.913      |     259.078      | 198.208 |  FALSE  |               0.408               |          0.131           |            0.013            |
| 2011 | ENP Calves                   |       209.9        | 214.421 | 214.104  |     172.804      |     260.524      | 196.721 |  FALSE  |               0.417               |          0.144           |            0.022            |
| 2011 | PCFG Calves only             |       209.9        | 209.056 | 208.167  |     174.220      |     248.769      | 193.784 |  FALSE  |               0.538               |          0.171           |            0.018            |
| 2011 | ENP Strandings only          |       209.9        | 213.841 | 213.964  |     176.661      |     248.473      | 199.847 |  FALSE  |               0.401               |          0.108           |            0.015            |
| 2011 | PCFG Calves + ENP Strandings |       209.9        | 210.256 | 209.837  |     173.691      |     250.340      | 195.241 |  FALSE  |               0.503               |          0.152           |            0.018            |
| 2012 | Base                         |       223.5        | 217.462 | 217.271  |     179.918      |     255.988      | 203.490 |  FALSE  |               0.650               |          0.078           |            0.009            |
| 2012 | AR1v1                        |       223.5        | 228.133 | 227.300  |     190.881      |     271.906      | 213.322 |  FALSE  |               0.407               |          0.028           |            0.005            |
| 2012 | AR1v2                        |       223.5        | 218.527 | 217.822  |     178.648      |     261.657      | 202.088 |  FALSE  |               0.621               |          0.089           |            0.009            |
| 2012 | ENP Calves                   |       223.5        | 217.826 | 217.671  |     179.746      |     259.211      | 201.180 |  FALSE  |               0.625               |          0.095           |            0.010            |
| 2012 | PCFG Calves only             |       223.5        | 218.453 | 217.517  |     181.805      |     259.914      | 203.408 |  FALSE  |               0.636               |          0.068           |            0.008            |
| 2012 | ENP Strandings only          |       223.5        | 219.801 | 220.235  |     182.969      |     256.794      | 205.417 |  FALSE  |               0.583               |          0.063           |            0.008            |
| 2012 | PCFG Calves + ENP Strandings |       223.5        | 220.369 | 219.504  |     184.020      |     262.404      | 204.809 |  FALSE  |               0.588               |          0.061           |            0.006            |
| 2013 | Base                         |       243.0        | 226.232 | 225.705  |     190.458      |     265.025      | 211.452 |  FALSE  |               0.827               |          0.030           |            0.004            |
| 2013 | AR1v1                        |       243.0        | 238.894 | 237.415  |     199.802      |     286.503      | 223.813 |  FALSE  |               0.637               |          0.013           |            0.002            |
| 2013 | AR1v2                        |       243.0        | 227.570 | 226.756  |     187.901      |     274.075      | 210.749 |  FALSE  |               0.798               |          0.038           |            0.005            |
| 2013 | ENP Calves                   |       243.0        | 226.032 | 225.589  |     189.808      |     267.097      | 210.925 |  FALSE  |               0.836               |          0.033           |            0.003            |
| 2013 | PCFG Calves only             |       243.0        | 231.932 | 230.773  |     191.379      |     278.685      | 215.248 |  FALSE  |               0.729               |          0.027           |            0.002            |
| 2013 | ENP Strandings only          |       243.0        | 229.423 | 228.481  |     194.757      |     268.800      | 215.514 |  FALSE  |               0.806               |          0.020           |            0.001            |
| 2013 | PCFG Calves + ENP Strandings |       243.0        | 233.441 | 232.384  |     195.053      |     278.059      | 217.858 |  FALSE  |               0.717               |          0.019           |            0.003            |
| 2014 | Base                         |       250.7        | 237.563 | 236.253  |     202.687      |     281.043      | 221.685 |  FALSE  |               0.775               |          0.006           |            0.000            |
| 2014 | AR1v1                        |       250.7        | 247.082 | 244.678  |     209.622      |     293.458      | 230.635 |  FALSE  |               0.612               |          0.003           |            0.001            |
| 2014 | AR1v2                        |       250.7        | 238.031 | 237.200  |     198.760      |     283.083      | 221.138 |  FALSE  |               0.749               |          0.013           |            0.001            |
| 2014 | ENP Calves                   |       250.7        | 237.895 | 236.726  |     201.212      |     285.917      | 220.851 |  FALSE  |               0.762               |          0.009           |            0.001            |
| 2014 | PCFG Calves only             |       250.7        | 248.967 | 247.236  |     207.955      |     299.240      | 230.239 |  FALSE  |               0.566               |          0.005           |            0.000            |
| 2014 | ENP Strandings only          |       250.7        | 240.318 | 238.856  |     204.324      |     280.481      | 225.251 |  FALSE  |               0.741               |          0.005           |            0.000            |
| 2014 | PCFG Calves + ENP Strandings |       250.7        | 249.115 | 246.986  |     207.252      |     302.518      | 230.308 |  FALSE  |               0.571               |          0.005           |            0.001            |
| 2015 | Base                         |       253.2        | 243.759 | 241.745  |     206.670      |     290.244      | 227.000 |  FALSE  |               0.698               |          0.004           |            0.000            |
| 2015 | AR1v1                        |       253.2        | 248.364 | 246.402  |     210.453      |     297.609      | 231.833 |  FALSE  |               0.645               |          0.003           |            0.001            |
| 2015 | AR1v2                        |       253.2        | 242.397 | 240.528  |     201.550      |     293.205      | 224.310 |  FALSE  |               0.714               |          0.010           |            0.000            |
| 2015 | ENP Calves                   |       253.2        | 243.290 | 241.317  |     200.884      |     295.148      | 225.466 |  FALSE  |               0.713               |          0.009           |            0.002            |
| 2015 | PCFG Calves only             |       253.2        | 254.197 | 252.238  |     211.583      |     306.715      | 234.692 |  FALSE  |               0.518               |          0.004           |            0.000            |
| 2015 | ENP Strandings only          |       253.2        | 245.688 | 244.260  |     208.627      |     290.473      | 229.184 |  FALSE  |               0.680               |          0.003           |            0.000            |
| 2015 | PCFG Calves + ENP Strandings |       253.2        | 254.840 | 252.569  |     211.842      |     309.335      | 235.674 |  FALSE  |               0.508               |          0.002           |            0.000            |
| 2016 | Base                         |       254.9        | 244.516 | 242.334  |     209.210      |     292.166      | 227.140 |  FALSE  |               0.712               |          0.003           |            0.000            |
| 2016 | AR1v1                        |       254.9        | 243.759 | 242.084  |     204.913      |     293.835      | 226.712 |  FALSE  |               0.736               |          0.006           |            0.001            |
| 2016 | AR1v2                        |       254.9        | 242.712 | 240.914  |     200.786      |     295.549      | 225.084 |  FALSE  |               0.740               |          0.011           |            0.001            |
| 2016 | ENP Calves                   |       254.9        | 244.661 | 242.204  |     203.555      |     299.877      | 226.253 |  FALSE  |               0.711               |          0.005           |            0.000            |
| 2016 | PCFG Calves only             |       254.9        | 251.016 | 249.241  |     209.427      |     303.015      | 232.528 |  FALSE  |               0.606               |          0.003           |            0.000            |
| 2016 | ENP Strandings only          |       254.9        | 245.389 | 243.228  |     209.610      |     291.905      | 229.005 |  FALSE  |               0.712               |          0.004           |            0.001            |
| 2016 | PCFG Calves + ENP Strandings |       254.9        | 251.077 | 249.660  |     211.446      |     300.609      | 233.095 |  FALSE  |               0.604               |          0.001           |            0.000            |
| 2017 | Base                         |       228.3        | 240.594 | 237.931  |     201.948      |     291.142      | 223.329 |  FALSE  |               0.295               |          0.010           |            0.001            |
| 2017 | AR1v1                        |       228.3        | 233.570 | 232.904  |     194.100      |     277.459      | 217.915 |  FALSE  |               0.395               |          0.022           |            0.003            |
| 2017 | AR1v2                        |       228.3        | 238.394 | 236.673  |     195.839      |     291.298      | 219.974 |  FALSE  |               0.345               |          0.018           |            0.002            |
| 2017 | ENP Calves                   |       228.3        | 240.428 | 237.731  |     197.670      |     295.860      | 221.469 |  FALSE  |               0.316               |          0.013           |            0.000            |
| 2017 | PCFG Calves only             |       228.3        | 243.226 | 241.974  |     203.623      |     290.530      | 225.512 |  FALSE  |               0.243               |          0.009           |            0.000            |
| 2017 | ENP Strandings only          |       228.3        | 243.684 | 242.114  |     204.922      |     290.357      | 227.403 |  FALSE  |               0.216               |          0.004           |            0.001            |
| 2017 | PCFG Calves + ENP Strandings |       228.3        | 244.955 | 244.096  |     204.490      |     292.391      | 227.470 |  FALSE  |               0.215               |          0.007           |            0.000            |
| 2018 | Base                         |       214.9        | 230.268 | 229.174  |     191.466      |     276.456      | 213.966 |  FALSE  |               0.215               |          0.027           |            0.003            |
| 2018 | AR1v1                        |       214.9        | 223.418 | 223.085  |     184.796      |     265.790      | 208.178 |  FALSE  |               0.318               |          0.055           |            0.007            |
| 2018 | AR1v2                        |       214.9        | 228.355 | 227.245  |     184.906      |     274.927      | 211.244 |  FALSE  |               0.260               |          0.042           |            0.006            |
| 2018 | ENP Calves                   |       214.9        | 229.909 | 228.061  |     188.718      |     279.856      | 211.988 |  FALSE  |               0.245               |          0.038           |            0.003            |
| 2018 | PCFG Calves only             |       214.9        | 231.608 | 230.785  |     191.479      |     275.401      | 215.197 |  FALSE  |               0.194               |          0.026           |            0.003            |
| 2018 | ENP Strandings only          |       214.9        | 233.232 | 232.780  |     193.373      |     276.555      | 218.006 |  FALSE  |               0.160               |          0.022           |            0.001            |
| 2018 | PCFG Calves + ENP Strandings |       214.9        | 233.620 | 233.080  |     193.937      |     278.859      | 217.125 |  FALSE  |               0.171               |          0.020           |            0.001            |
| 2019 | Base                         |       212.0        | 222.427 | 222.318  |     185.206      |     261.471      | 207.304 |  FALSE  |               0.280               |          0.049           |            0.005            |
| 2019 | AR1v1                        |       212.0        | 215.168 | 214.762  |     177.768      |     254.938      | 200.699 |  FALSE  |               0.423               |          0.098           |            0.013            |
| 2019 | AR1v2                        |       212.0        | 220.029 | 219.316  |     178.464      |     265.616      | 203.054 |  FALSE  |               0.346               |          0.077           |            0.009            |
| 2019 | ENP Calves                   |       212.0        | 221.605 | 220.649  |     181.887      |     267.758      | 204.662 |  FALSE  |               0.322               |          0.068           |            0.008            |
| 2019 | PCFG Calves only             |       212.0        | 219.775 | 219.532  |     179.956      |     261.054      | 203.579 |  FALSE  |               0.348               |          0.078           |            0.010            |
| 2019 | ENP Strandings only          |       212.0        | 226.155 | 226.299  |     186.878      |     266.334      | 210.694 |  FALSE  |               0.219               |          0.043           |            0.003            |
| 2019 | PCFG Calves + ENP Strandings |       212.0        | 222.995 | 223.513  |     181.154      |     264.926      | 207.067 |  FALSE  |               0.282               |          0.066           |            0.012            |
| 2020 | Base                         |       206.9        | 216.152 | 215.827  |     179.796      |     255.144      | 200.957 |  FALSE  |               0.297               |          0.096           |            0.011            |
| 2020 | AR1v1                        |       206.9        | 209.842 | 209.097  |     174.208      |     249.666      | 195.072 |  FALSE  |               0.440               |          0.161           |            0.015            |
| 2020 | AR1v2                        |       206.9        | 215.200 | 214.514  |     175.823      |     256.764      | 199.026 |  FALSE  |               0.336               |          0.114           |            0.016            |
| 2020 | ENP Calves                   |       206.9        | 215.569 | 214.847  |     176.586      |     256.661      | 199.924 |  FALSE  |               0.332               |          0.109           |            0.013            |
| 2020 | PCFG Calves only             |       206.9        | 213.993 | 213.176  |     177.771      |     253.018      | 198.270 |  FALSE  |               0.355               |          0.116           |            0.013            |
| 2020 | ENP Strandings only          |       206.9        | 204.458 | 204.208  |     165.776      |     248.124      | 187.976 |  FALSE  |               0.556               |          0.262           |            0.046            |
| 2020 | PCFG Calves + ENP Strandings |       206.9        | 205.111 | 204.206  |     164.798      |     249.484      | 188.921 |  FALSE  |               0.557               |          0.249           |            0.047            |
| 2021 | Base                         |       210.1        | 211.887 | 211.868  |     176.738      |     251.386      | 195.878 |  FALSE  |               0.468               |          0.147           |            0.012            |
| 2021 | AR1v1                        |       210.1        | 204.365 | 203.839  |     168.548      |     245.298      | 189.418 |  FALSE  |               0.651               |          0.248           |            0.033            |
| 2021 | AR1v2                        |       210.1        | 211.665 | 211.028  |     172.444      |     254.382      | 195.077 |  FALSE  |               0.480               |          0.162           |            0.023            |
| 2021 | ENP Calves                   |       210.1        | 212.088 | 211.089  |     172.282      |     255.867      | 195.251 |  FALSE  |               0.482               |          0.155           |            0.022            |
| 2021 | PCFG Calves only             |       210.1        | 210.078 | 209.671  |     172.531      |     249.627      | 194.748 |  FALSE  |               0.509               |          0.164           |            0.019            |
| 2021 | ENP Strandings only          |       210.1        | 201.824 | 200.683  |     165.142      |     244.235      | 185.219 |  FALSE  |               0.683               |          0.314           |            0.047            |
| 2021 | PCFG Calves + ENP Strandings |       210.1        | 202.950 | 201.819  |     164.655      |     245.869      | 187.112 |  FALSE  |               0.655               |          0.290           |            0.049            |
| 2022 | Base                         |       202.0        | 210.319 | 209.917  |     170.389      |     252.092      | 193.944 |  FALSE  |               0.335               |          0.170           |            0.026            |
| 2022 | AR1v1                        |       202.0        | 200.704 | 199.157  |     159.748      |     248.448      | 182.796 |  FALSE  |               0.563               |          0.353           |            0.077            |
| 2022 | AR1v2                        |       202.0        | 210.176 | 209.075  |     169.389      |     256.883      | 193.025 |  FALSE  |               0.359               |          0.187           |            0.032            |
| 2022 | ENP Calves                   |       202.0        | 211.233 | 210.371  |     168.617      |     258.559      | 192.660 |  FALSE  |               0.353               |          0.194           |            0.033            |
| 2022 | PCFG Calves only             |       202.0        | 206.766 | 206.071  |     169.064      |     250.918      | 190.607 |  FALSE  |               0.416               |          0.221           |            0.031            |
| 2022 | ENP Strandings only          |       202.0        | 199.820 | 199.094  |     162.400      |     246.573      | 182.456 |  FALSE  |               0.556               |          0.362           |            0.069            |
| 2022 | PCFG Calves + ENP Strandings |       202.0        | 199.563 | 198.199  |     161.498      |     245.706      | 182.189 |  FALSE  |               0.571               |          0.379           |            0.078            |

</div>

Model predictions when projecting two years forward for all data years
from 2002 through 2021. RSS is the residual sum of squares.

<div class="cell-output-display">

| Year | Model                        | Abundance Estimate | Mean N  | Median N | Lower 95% CI (N) | Upper 95% CI (N) |  Nmin   | Closure | Percentile for Abundance Estimate | Prop Below Threshold (N) | Prop Below Threshold (Nmin) |
|:-----|:-----------------------------|:------------------:|:-------:|:--------:|:----------------:|:----------------:|:-------:|:-------:|:---------------------------------:|:------------------------:|:---------------------------:|
| 2003 | Base                         |       208.7        | 200.434 | 200.357  |     168.471      |     234.323      | 187.406 |  FALSE  |               0.716               |          0.295           |            0.033            |
| 2003 | AR1v1                        |       208.7        | 205.196 | 203.962  |     171.201      |     247.051      | 191.053 |  FALSE  |               0.618               |          0.218           |            0.024            |
| 2003 | AR1v2                        |       208.7        | 199.524 | 198.434  |     165.133      |     239.238      | 184.499 |  FALSE  |               0.711               |          0.348           |            0.056            |
| 2003 | ENP Calves                   |       208.7        | 199.364 | 198.992  |     166.550      |     235.846      | 185.357 |  FALSE  |               0.723               |          0.329           |            0.042            |
| 2003 | PCFG Calves only             |       208.7        | 199.689 | 199.053  |     166.971      |     235.516      | 186.530 |  FALSE  |               0.724               |          0.325           |            0.037            |
| 2003 | ENP Strandings only          |       208.7        | 202.108 | 201.602  |     171.813      |     233.348      | 189.585 |  FALSE  |               0.677               |          0.250           |            0.020            |
| 2003 | PCFG Calves + ENP Strandings |       208.7        | 200.982 | 200.559  |     166.355      |     236.541      | 187.160 |  FALSE  |               0.693               |          0.291           |            0.040            |
| 2004 | Base                         |       214.9        | 206.328 | 205.552  |     173.120      |     244.794      | 192.109 |  FALSE  |               0.714               |          0.198           |            0.020            |
| 2004 | AR1v1                        |       214.9        | 209.020 | 207.079  |     175.212      |     251.902      | 194.803 |  FALSE  |               0.674               |          0.153           |            0.014            |
| 2004 | AR1v2                        |       214.9        | 210.755 | 209.385  |     173.956      |     254.552      | 194.252 |  FALSE  |               0.609               |          0.165           |            0.015            |
| 2004 | ENP Calves                   |       214.9        | 205.906 | 205.293  |     170.098      |     247.592      | 191.471 |  FALSE  |               0.716               |          0.211           |            0.028            |
| 2004 | PCFG Calves only             |       214.9        | 205.940 | 204.995  |     171.610      |     245.408      | 191.489 |  FALSE  |               0.711               |          0.207           |            0.023            |
| 2004 | ENP Strandings only          |       214.9        | 208.994 | 207.993  |     175.547      |     248.298      | 195.385 |  FALSE  |               0.671               |          0.138           |            0.014            |
| 2004 | PCFG Calves + ENP Strandings |       214.9        | 208.598 | 207.811  |     173.763      |     248.635      | 193.948 |  FALSE  |               0.666               |          0.172           |            0.018            |
| 2005 | Base                         |       220.0        | 210.088 | 208.484  |     177.052      |     249.344      | 195.903 |  FALSE  |               0.738               |          0.136           |            0.011            |
| 2005 | AR1v1                        |       220.0        | 209.253 | 207.997  |     172.831      |     253.563      | 194.643 |  FALSE  |               0.752               |          0.159           |            0.023            |
| 2005 | AR1v2                        |       220.0        | 210.969 | 210.328  |     172.587      |     253.275      | 195.209 |  FALSE  |               0.705               |          0.154           |            0.021            |
| 2005 | ENP Calves                   |       220.0        | 209.801 | 208.536  |     174.083      |     253.238      | 194.465 |  FALSE  |               0.737               |          0.162           |            0.018            |
| 2005 | PCFG Calves only             |       220.0        | 208.843 | 207.715  |     176.627      |     248.323      | 193.908 |  FALSE  |               0.756               |          0.166           |            0.014            |
| 2005 | ENP Strandings only          |       220.0        | 211.320 | 210.053  |     178.958      |     251.971      | 197.370 |  FALSE  |               0.725               |          0.124           |            0.011            |
| 2005 | PCFG Calves + ENP Strandings |       220.0        | 210.338 | 209.358  |     175.127      |     248.637      | 195.920 |  FALSE  |               0.728               |          0.132           |            0.014            |
| 2006 | Base                         |       199.2        | 209.609 | 208.795  |     174.333      |     251.336      | 194.743 |  FALSE  |               0.284               |          0.159           |            0.017            |
| 2006 | AR1v1                        |       199.2        | 207.736 | 206.882  |     171.124      |     250.977      | 192.671 |  FALSE  |               0.324               |          0.190           |            0.024            |
| 2006 | AR1v2                        |       199.2        | 210.908 | 210.403  |     170.587      |     254.448      | 194.348 |  FALSE  |               0.278               |          0.169           |            0.027            |
| 2006 | ENP Calves                   |       199.2        | 209.505 | 207.880  |     173.260      |     254.424      | 192.973 |  FALSE  |               0.318               |          0.187           |            0.018            |
| 2006 | PCFG Calves only             |       199.2        | 206.707 | 205.440  |     170.377      |     250.324      | 191.047 |  FALSE  |               0.353               |          0.211           |            0.028            |
| 2006 | ENP Strandings only          |       199.2        | 210.706 | 209.427  |     176.911      |     250.371      | 196.866 |  FALSE  |               0.248               |          0.130           |            0.012            |
| 2006 | PCFG Calves + ENP Strandings |       199.2        | 207.945 | 207.461  |     171.825      |     249.450      | 192.371 |  FALSE  |               0.321               |          0.191           |            0.024            |
| 2007 | Base                         |       197.7        | 207.294 | 207.061  |     169.376      |     247.231      | 192.952 |  FALSE  |               0.286               |          0.185           |            0.029            |
| 2007 | AR1v1                        |       197.7        | 207.291 | 206.788  |     171.536      |     247.495      | 192.214 |  FALSE  |               0.299               |          0.197           |            0.024            |
| 2007 | AR1v2                        |       197.7        | 208.388 | 208.036  |     168.666      |     251.035      | 191.912 |  FALSE  |               0.295               |          0.201           |            0.032            |
| 2007 | ENP Calves                   |       197.7        | 206.331 | 205.450  |     168.060      |     249.205      | 189.886 |  FALSE  |               0.331               |          0.228           |            0.035            |
| 2007 | PCFG Calves only             |       197.7        | 204.025 | 203.150  |     168.098      |     243.595      | 188.751 |  FALSE  |               0.372               |          0.255           |            0.036            |
| 2007 | ENP Strandings only          |       197.7        | 206.815 | 206.483  |     173.125      |     243.940      | 192.715 |  FALSE  |               0.292               |          0.189           |            0.020            |
| 2007 | PCFG Calves + ENP Strandings |       197.7        | 203.911 | 203.811  |     168.533      |     241.900      | 188.876 |  FALSE  |               0.364               |          0.247           |            0.034            |
| 2008 | Base                         |       211.8        | 208.051 | 207.802  |     169.935      |     248.740      | 193.011 |  FALSE  |               0.592               |          0.188           |            0.028            |
| 2008 | AR1v1                        |       211.8        | 207.829 | 207.236  |     171.502      |     247.874      | 193.167 |  FALSE  |               0.603               |          0.179           |            0.023            |
| 2008 | AR1v2                        |       211.8        | 209.231 | 208.222  |     170.397      |     250.617      | 192.880 |  FALSE  |               0.577               |          0.187           |            0.027            |
| 2008 | ENP Calves                   |       211.8        | 206.497 | 206.101  |     169.204      |     248.071      | 190.883 |  FALSE  |               0.624               |          0.216           |            0.031            |
| 2008 | PCFG Calves only             |       211.8        | 203.492 | 202.931  |     167.536      |     242.764      | 188.268 |  FALSE  |               0.689               |          0.268           |            0.042            |
| 2008 | ENP Strandings only          |       211.8        | 207.630 | 207.138  |     172.373      |     246.336      | 193.697 |  FALSE  |               0.610               |          0.175           |            0.022            |
| 2008 | PCFG Calves + ENP Strandings |       211.8        | 204.026 | 203.744  |     166.583      |     245.465      | 188.689 |  FALSE  |               0.677               |          0.249           |            0.039            |
| 2009 | Base                         |       213.8        | 210.511 | 210.118  |     174.370      |     251.949      | 195.516 |  FALSE  |               0.591               |          0.153           |            0.018            |
| 2009 | AR1v1                        |       213.8        | 209.714 | 209.408  |     173.074      |     248.187      | 195.272 |  FALSE  |               0.610               |          0.147           |            0.019            |
| 2009 | AR1v2                        |       213.8        | 210.977 | 210.427  |     170.917      |     255.120      | 194.250 |  FALSE  |               0.573               |          0.164           |            0.025            |
| 2009 | ENP Calves                   |       213.8        | 210.043 | 209.349  |     173.761      |     250.091      | 194.800 |  FALSE  |               0.603               |          0.157           |            0.021            |
| 2009 | PCFG Calves only             |       213.8        | 205.076 | 204.604  |     168.502      |     244.694      | 190.115 |  FALSE  |               0.701               |          0.237           |            0.032            |
| 2009 | ENP Strandings only          |       213.8        | 213.159 | 213.087  |     176.517      |     251.278      | 198.524 |  FALSE  |               0.519               |          0.118           |            0.016            |
| 2009 | PCFG Calves + ENP Strandings |       213.8        | 208.566 | 208.096  |     171.852      |     251.102      | 192.543 |  FALSE  |               0.623               |          0.193           |            0.021            |
| 2010 | Base                         |       205.3        | 211.895 | 211.057  |     175.200      |     252.207      | 196.591 |  FALSE  |               0.365               |          0.145           |            0.017            |
| 2010 | AR1v1                        |       205.3        | 211.698 | 212.044  |     172.035      |     251.151      | 196.503 |  FALSE  |               0.351               |          0.143           |            0.022            |
| 2010 | AR1v2                        |       205.3        | 212.490 | 211.587  |     175.790      |     253.274      | 196.017 |  FALSE  |               0.370               |          0.143           |            0.013            |
| 2010 | ENP Calves                   |       205.3        | 212.031 | 211.135  |     173.433      |     255.855      | 195.283 |  FALSE  |               0.373               |          0.152           |            0.020            |
| 2010 | PCFG Calves only             |       205.3        | 207.193 | 206.200  |     170.224      |     249.205      | 191.321 |  FALSE  |               0.479               |          0.209           |            0.028            |
| 2010 | ENP Strandings only          |       205.3        | 214.347 | 214.421  |     178.736      |     252.684      | 199.486 |  FALSE  |               0.300               |          0.105           |            0.010            |
| 2010 | PCFG Calves + ENP Strandings |       205.3        | 210.097 | 209.365  |     172.355      |     252.401      | 194.304 |  FALSE  |               0.410               |          0.165           |            0.022            |
| 2011 | Base                         |       209.9        | 213.210 | 213.413  |     176.704      |     251.150      | 198.470 |  FALSE  |               0.422               |          0.123           |            0.013            |
| 2011 | AR1v1                        |       209.9        | 217.820 | 217.834  |     180.712      |     255.829      | 203.572 |  FALSE  |               0.321               |          0.082           |            0.011            |
| 2011 | AR1v2                        |       209.9        | 214.658 | 214.540  |     174.913      |     259.078      | 198.208 |  FALSE  |               0.408               |          0.131           |            0.013            |
| 2011 | ENP Calves                   |       209.9        | 214.421 | 214.104  |     172.804      |     260.524      | 196.721 |  FALSE  |               0.417               |          0.144           |            0.022            |
| 2011 | PCFG Calves only             |       209.9        | 209.056 | 208.167  |     174.220      |     248.769      | 193.784 |  FALSE  |               0.538               |          0.171           |            0.018            |
| 2011 | ENP Strandings only          |       209.9        | 213.841 | 213.964  |     176.661      |     248.473      | 199.847 |  FALSE  |               0.401               |          0.108           |            0.015            |
| 2011 | PCFG Calves + ENP Strandings |       209.9        | 210.256 | 209.837  |     173.691      |     250.340      | 195.241 |  FALSE  |               0.503               |          0.152           |            0.018            |
| 2012 | Base                         |       223.5        | 217.462 | 217.271  |     179.918      |     255.988      | 203.490 |  FALSE  |               0.650               |          0.078           |            0.009            |
| 2012 | AR1v1                        |       223.5        | 228.133 | 227.300  |     190.881      |     271.906      | 213.322 |  FALSE  |               0.407               |          0.028           |            0.005            |
| 2012 | AR1v2                        |       223.5        | 218.527 | 217.822  |     178.648      |     261.657      | 202.088 |  FALSE  |               0.621               |          0.089           |            0.009            |
| 2012 | ENP Calves                   |       223.5        | 217.826 | 217.671  |     179.746      |     259.211      | 201.180 |  FALSE  |               0.625               |          0.095           |            0.010            |
| 2012 | PCFG Calves only             |       223.5        | 218.453 | 217.517  |     181.805      |     259.914      | 203.408 |  FALSE  |               0.636               |          0.068           |            0.008            |
| 2012 | ENP Strandings only          |       223.5        | 219.801 | 220.235  |     182.969      |     256.794      | 205.417 |  FALSE  |               0.583               |          0.063           |            0.008            |
| 2012 | PCFG Calves + ENP Strandings |       223.5        | 220.369 | 219.504  |     184.020      |     262.404      | 204.809 |  FALSE  |               0.588               |          0.061           |            0.006            |
| 2013 | Base                         |       243.0        | 226.232 | 225.705  |     190.458      |     265.025      | 211.452 |  FALSE  |               0.827               |          0.030           |            0.004            |
| 2013 | AR1v1                        |       243.0        | 238.894 | 237.415  |     199.802      |     286.503      | 223.813 |  FALSE  |               0.637               |          0.013           |            0.002            |
| 2013 | AR1v2                        |       243.0        | 227.570 | 226.756  |     187.901      |     274.075      | 210.749 |  FALSE  |               0.798               |          0.038           |            0.005            |
| 2013 | ENP Calves                   |       243.0        | 226.032 | 225.589  |     189.808      |     267.097      | 210.925 |  FALSE  |               0.836               |          0.033           |            0.003            |
| 2013 | PCFG Calves only             |       243.0        | 231.932 | 230.773  |     191.379      |     278.685      | 215.248 |  FALSE  |               0.729               |          0.027           |            0.002            |
| 2013 | ENP Strandings only          |       243.0        | 229.423 | 228.481  |     194.757      |     268.800      | 215.514 |  FALSE  |               0.806               |          0.020           |            0.001            |
| 2013 | PCFG Calves + ENP Strandings |       243.0        | 233.441 | 232.384  |     195.053      |     278.059      | 217.858 |  FALSE  |               0.717               |          0.019           |            0.003            |
| 2014 | Base                         |       250.7        | 237.563 | 236.253  |     202.687      |     281.043      | 221.685 |  FALSE  |               0.775               |          0.006           |            0.000            |
| 2014 | AR1v1                        |       250.7        | 247.082 | 244.678  |     209.622      |     293.458      | 230.635 |  FALSE  |               0.612               |          0.003           |            0.001            |
| 2014 | AR1v2                        |       250.7        | 238.031 | 237.200  |     198.760      |     283.083      | 221.138 |  FALSE  |               0.749               |          0.013           |            0.001            |
| 2014 | ENP Calves                   |       250.7        | 237.895 | 236.726  |     201.212      |     285.917      | 220.851 |  FALSE  |               0.762               |          0.009           |            0.001            |
| 2014 | PCFG Calves only             |       250.7        | 248.967 | 247.236  |     207.955      |     299.240      | 230.239 |  FALSE  |               0.566               |          0.005           |            0.000            |
| 2014 | ENP Strandings only          |       250.7        | 240.318 | 238.856  |     204.324      |     280.481      | 225.251 |  FALSE  |               0.741               |          0.005           |            0.000            |
| 2014 | PCFG Calves + ENP Strandings |       250.7        | 249.115 | 246.986  |     207.252      |     302.518      | 230.308 |  FALSE  |               0.571               |          0.005           |            0.001            |
| 2015 | Base                         |       253.2        | 243.759 | 241.745  |     206.670      |     290.244      | 227.000 |  FALSE  |               0.698               |          0.004           |            0.000            |
| 2015 | AR1v1                        |       253.2        | 248.364 | 246.402  |     210.453      |     297.609      | 231.833 |  FALSE  |               0.645               |          0.003           |            0.001            |
| 2015 | AR1v2                        |       253.2        | 242.397 | 240.528  |     201.550      |     293.205      | 224.310 |  FALSE  |               0.714               |          0.010           |            0.000            |
| 2015 | ENP Calves                   |       253.2        | 243.290 | 241.317  |     200.884      |     295.148      | 225.466 |  FALSE  |               0.713               |          0.009           |            0.002            |
| 2015 | PCFG Calves only             |       253.2        | 254.197 | 252.238  |     211.583      |     306.715      | 234.692 |  FALSE  |               0.518               |          0.004           |            0.000            |
| 2015 | ENP Strandings only          |       253.2        | 245.688 | 244.260  |     208.627      |     290.473      | 229.184 |  FALSE  |               0.680               |          0.003           |            0.000            |
| 2015 | PCFG Calves + ENP Strandings |       253.2        | 254.840 | 252.569  |     211.842      |     309.335      | 235.674 |  FALSE  |               0.508               |          0.002           |            0.000            |
| 2016 | Base                         |       254.9        | 244.516 | 242.334  |     209.210      |     292.166      | 227.140 |  FALSE  |               0.712               |          0.003           |            0.000            |
| 2016 | AR1v1                        |       254.9        | 243.759 | 242.084  |     204.913      |     293.835      | 226.712 |  FALSE  |               0.736               |          0.006           |            0.001            |
| 2016 | AR1v2                        |       254.9        | 242.712 | 240.914  |     200.786      |     295.549      | 225.084 |  FALSE  |               0.740               |          0.011           |            0.001            |
| 2016 | ENP Calves                   |       254.9        | 244.661 | 242.204  |     203.555      |     299.877      | 226.253 |  FALSE  |               0.711               |          0.005           |            0.000            |
| 2016 | PCFG Calves only             |       254.9        | 251.016 | 249.241  |     209.427      |     303.015      | 232.528 |  FALSE  |               0.606               |          0.003           |            0.000            |
| 2016 | ENP Strandings only          |       254.9        | 245.389 | 243.228  |     209.610      |     291.905      | 229.005 |  FALSE  |               0.712               |          0.004           |            0.001            |
| 2016 | PCFG Calves + ENP Strandings |       254.9        | 251.077 | 249.660  |     211.446      |     300.609      | 233.095 |  FALSE  |               0.604               |          0.001           |            0.000            |
| 2017 | Base                         |       228.3        | 240.594 | 237.931  |     201.948      |     291.142      | 223.329 |  FALSE  |               0.295               |          0.010           |            0.001            |
| 2017 | AR1v1                        |       228.3        | 233.570 | 232.904  |     194.100      |     277.459      | 217.915 |  FALSE  |               0.395               |          0.022           |            0.003            |
| 2017 | AR1v2                        |       228.3        | 238.394 | 236.673  |     195.839      |     291.298      | 219.974 |  FALSE  |               0.345               |          0.018           |            0.002            |
| 2017 | ENP Calves                   |       228.3        | 240.428 | 237.731  |     197.670      |     295.860      | 221.469 |  FALSE  |               0.316               |          0.013           |            0.000            |
| 2017 | PCFG Calves only             |       228.3        | 243.226 | 241.974  |     203.623      |     290.530      | 225.512 |  FALSE  |               0.243               |          0.009           |            0.000            |
| 2017 | ENP Strandings only          |       228.3        | 243.684 | 242.114  |     204.922      |     290.357      | 227.403 |  FALSE  |               0.216               |          0.004           |            0.001            |
| 2017 | PCFG Calves + ENP Strandings |       228.3        | 244.955 | 244.096  |     204.490      |     292.391      | 227.470 |  FALSE  |               0.215               |          0.007           |            0.000            |
| 2018 | Base                         |       214.9        | 230.268 | 229.174  |     191.466      |     276.456      | 213.966 |  FALSE  |               0.215               |          0.027           |            0.003            |
| 2018 | AR1v1                        |       214.9        | 223.418 | 223.085  |     184.796      |     265.790      | 208.178 |  FALSE  |               0.318               |          0.055           |            0.007            |
| 2018 | AR1v2                        |       214.9        | 228.355 | 227.245  |     184.906      |     274.927      | 211.244 |  FALSE  |               0.260               |          0.042           |            0.006            |
| 2018 | ENP Calves                   |       214.9        | 229.909 | 228.061  |     188.718      |     279.856      | 211.988 |  FALSE  |               0.245               |          0.038           |            0.003            |
| 2018 | PCFG Calves only             |       214.9        | 231.608 | 230.785  |     191.479      |     275.401      | 215.197 |  FALSE  |               0.194               |          0.026           |            0.003            |
| 2018 | ENP Strandings only          |       214.9        | 233.232 | 232.780  |     193.373      |     276.555      | 218.006 |  FALSE  |               0.160               |          0.022           |            0.001            |
| 2018 | PCFG Calves + ENP Strandings |       214.9        | 233.620 | 233.080  |     193.937      |     278.859      | 217.125 |  FALSE  |               0.171               |          0.020           |            0.001            |
| 2019 | Base                         |       212.0        | 222.427 | 222.318  |     185.206      |     261.471      | 207.304 |  FALSE  |               0.280               |          0.049           |            0.005            |
| 2019 | AR1v1                        |       212.0        | 215.168 | 214.762  |     177.768      |     254.938      | 200.699 |  FALSE  |               0.423               |          0.098           |            0.013            |
| 2019 | AR1v2                        |       212.0        | 220.029 | 219.316  |     178.464      |     265.616      | 203.054 |  FALSE  |               0.346               |          0.077           |            0.009            |
| 2019 | ENP Calves                   |       212.0        | 221.605 | 220.649  |     181.887      |     267.758      | 204.662 |  FALSE  |               0.322               |          0.068           |            0.008            |
| 2019 | PCFG Calves only             |       212.0        | 219.775 | 219.532  |     179.956      |     261.054      | 203.579 |  FALSE  |               0.348               |          0.078           |            0.010            |
| 2019 | ENP Strandings only          |       212.0        | 226.155 | 226.299  |     186.878      |     266.334      | 210.694 |  FALSE  |               0.219               |          0.043           |            0.003            |
| 2019 | PCFG Calves + ENP Strandings |       212.0        | 222.995 | 223.513  |     181.154      |     264.926      | 207.067 |  FALSE  |               0.282               |          0.066           |            0.012            |
| 2020 | Base                         |       206.9        | 216.152 | 215.827  |     179.796      |     255.144      | 200.957 |  FALSE  |               0.297               |          0.096           |            0.011            |
| 2020 | AR1v1                        |       206.9        | 209.842 | 209.097  |     174.208      |     249.666      | 195.072 |  FALSE  |               0.440               |          0.161           |            0.015            |
| 2020 | AR1v2                        |       206.9        | 215.200 | 214.514  |     175.823      |     256.764      | 199.026 |  FALSE  |               0.336               |          0.114           |            0.016            |
| 2020 | ENP Calves                   |       206.9        | 215.569 | 214.847  |     176.586      |     256.661      | 199.924 |  FALSE  |               0.332               |          0.109           |            0.013            |
| 2020 | PCFG Calves only             |       206.9        | 213.993 | 213.176  |     177.771      |     253.018      | 198.270 |  FALSE  |               0.355               |          0.116           |            0.013            |
| 2020 | ENP Strandings only          |       206.9        | 204.458 | 204.208  |     165.776      |     248.124      | 187.976 |  FALSE  |               0.556               |          0.262           |            0.046            |
| 2020 | PCFG Calves + ENP Strandings |       206.9        | 205.111 | 204.206  |     164.798      |     249.484      | 188.921 |  FALSE  |               0.557               |          0.249           |            0.047            |
| 2021 | Base                         |       210.1        | 211.887 | 211.868  |     176.738      |     251.386      | 195.878 |  FALSE  |               0.468               |          0.147           |            0.012            |
| 2021 | AR1v1                        |       210.1        | 204.365 | 203.839  |     168.548      |     245.298      | 189.418 |  FALSE  |               0.651               |          0.248           |            0.033            |
| 2021 | AR1v2                        |       210.1        | 211.665 | 211.028  |     172.444      |     254.382      | 195.077 |  FALSE  |               0.480               |          0.162           |            0.023            |
| 2021 | ENP Calves                   |       210.1        | 212.088 | 211.089  |     172.282      |     255.867      | 195.251 |  FALSE  |               0.482               |          0.155           |            0.022            |
| 2021 | PCFG Calves only             |       210.1        | 210.078 | 209.671  |     172.531      |     249.627      | 194.748 |  FALSE  |               0.509               |          0.164           |            0.019            |
| 2021 | ENP Strandings only          |       210.1        | 201.824 | 200.683  |     165.142      |     244.235      | 185.219 |  FALSE  |               0.683               |          0.314           |            0.047            |
| 2021 | PCFG Calves + ENP Strandings |       210.1        | 202.950 | 201.819  |     164.655      |     245.869      | 187.112 |  FALSE  |               0.655               |          0.290           |            0.049            |
| 2022 | Base                         |       202.0        | 210.319 | 209.917  |     170.389      |     252.092      | 193.944 |  FALSE  |               0.335               |          0.170           |            0.026            |
| 2022 | AR1v1                        |       202.0        | 200.704 | 199.157  |     159.748      |     248.448      | 182.796 |  FALSE  |               0.563               |          0.353           |            0.077            |
| 2022 | AR1v2                        |       202.0        | 210.176 | 209.075  |     169.389      |     256.883      | 193.025 |  FALSE  |               0.359               |          0.187           |            0.032            |
| 2022 | ENP Calves                   |       202.0        | 211.233 | 210.371  |     168.617      |     258.559      | 192.660 |  FALSE  |               0.353               |          0.194           |            0.033            |
| 2022 | PCFG Calves only             |       202.0        | 206.766 | 206.071  |     169.064      |     250.918      | 190.607 |  FALSE  |               0.416               |          0.221           |            0.031            |
| 2022 | ENP Strandings only          |       202.0        | 199.820 | 199.094  |     162.400      |     246.573      | 182.456 |  FALSE  |               0.556               |          0.362           |            0.069            |
| 2022 | PCFG Calves + ENP Strandings |       202.0        | 199.563 | 198.199  |     161.498      |     245.706      | 182.189 |  FALSE  |               0.571               |          0.379           |            0.078            |

</div>

Model fit statistics for predicting one year forward for all data years
from 2002 through 2021. RSS is the residual sum of squares.

<div class="cell-output-display">

| Model                        |   RSS    | Mean Percentile (N) | Median Percentile (N) | Lower 95% CI for Percentile (N) | Upper 95% CI for Percentile (N) | Mean Prop Below Threshold (N) | Mean Prop Below Threshold (Nmin) | Number of Closures |
|:-----------------------------|:--------:|:-------------------:|:---------------------:|:-------------------------------:|:-------------------------------:|:-----------------------------:|:--------------------------------:|:------------------:|
| AR1v1                        | 815.825  |        0.519        |         0.583         |              0.308              |              0.744              |             0.123             |              0.017               |         0          |
| PCFG Calves only             | 1334.201 |        0.522        |         0.528         |              0.218              |              0.743              |             0.138             |              0.017               |         0          |
| PCFG Calves + ENP Strandings | 1348.650 |        0.521        |         0.571         |              0.192              |              0.722              |             0.145             |              0.021               |         0          |
| AR1v2                        | 1697.251 |        0.514        |         0.526         |              0.268              |              0.775              |             0.121             |              0.017               |         0          |
| ENP Strandings only          | 1789.330 |        0.518        |         0.569         |              0.186              |              0.775              |             0.122             |              0.016               |         0          |
| Base                         | 1860.969 |        0.513        |         0.530         |              0.246              |              0.803              |             0.110             |              0.013               |         0          |
| ENP Calves                   | 1876.833 |        0.527        |         0.542         |              0.279              |              0.801              |             0.126             |              0.017               |         0          |

</div>

Model fit statistics for predicting two years forward for all data years
from 2002 through 2021. RSS is the residual sum of squares.

<div class="cell-output-display">

| Model                        |   RSS    | Mean Percentile (N) | Median Percentile (N) | Lower 95% CI for Percentile (N) | Upper 95% CI for Percentile (N) | Mean Prop Below Threshold (N) | Mean Prop Below Threshold (Nmin) | Number of Closures |
|:-----------------------------|:--------:|:-------------------:|:---------------------:|:-------------------------------:|:-------------------------------:|:-----------------------------:|:--------------------------------:|:------------------:|
| AR1v1                        | 737.673  |        0.500        |         0.512         |              0.365              |              0.604              |             0.203             |              0.071               |         0          |
| PCFG Calves + ENP Strandings | 2669.565 |        0.509        |         0.512         |              0.151              |              0.725              |             0.188             |              0.048               |         0          |
| PCFG Calves only             | 2920.820 |        0.509        |         0.519         |              0.173              |              0.776              |             0.185             |              0.042               |         0          |
| ENP Strandings only          | 3352.119 |        0.509        |         0.576         |              0.153              |              0.817              |             0.154             |              0.036               |         0          |
| AR1v2                        | 3697.296 |        0.509        |         0.469         |              0.228              |              0.826              |             0.163             |              0.039               |         0          |
| Base                         | 4156.324 |        0.504        |         0.465         |              0.173              |              0.854              |             0.147             |              0.030               |         0          |
| ENP Calves                   | 4200.056 |        0.512        |         0.459         |              0.209              |              0.846              |             0.165             |              0.040               |         0          |

</div>

![](README_files/figure-commonmark/retroFig-1.png)

## Model-specific trends and projections

Note, the number of PCFG calves used in projected years (models
Calves/Strandings and Calves only) are not accurate and likely represent
underestimates of reality.

![](README_files/figure-commonmark/trendFig-1.png)

## Model-specific predictions for Y<sub>final</sub> + 2

![](README_files/figure-commonmark/fig-final-proj-1.png)
