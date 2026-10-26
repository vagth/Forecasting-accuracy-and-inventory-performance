# Forecast accuracy and inventory performance

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.ejor.2024.12.033-blue)](https://doi.org/10.1016/j.ejor.2024.12.033)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

R code accompanying the paper:

> Theodorou, E., Spiliotis, E., & Assimakopoulos, V. (2025). Forecast accuracy and inventory performance: Insights on their relationship from the M5 competition data. *European Journal of Operational Research*, 322(2), 414–426. https://doi.org/10.1016/j.ejor.2024.12.033

## Authors

- **Evangelos Theodorou** — Forecasting and Strategy Unit, School of Electrical and Computer Engineering, National Technical University of Athens, Greece
- **Evangelos Spiliotis** — Forecasting and Strategy Unit, School of Electrical and Computer Engineering, National Technical University of Athens, Greece
- **Vassilios Assimakopoulos** — Forecasting and Strategy Unit, School of Electrical and Computer Engineering, National Technical University of Athens, Greece

## Overview

The study examines how forecast accuracy relates to inventory performance (holding, ordering and lost-sales costs) using the M5 competition data. Forecasts from twelve methods (LightGBM, Naive, Seasonal Naive, Moving Average, SES, Croston, SBA, TSB, ETS, ARIMA, ADIDA, iMAPA) drive a periodic-review order-up-to inventory simulation across review periods, lead times and target service levels, and the resulting inventory metrics are compared with standard accuracy measures.

## Repository structure

| File | Description |
|---|---|
| `lgbm/lgbmForecasts.R` | Feature engineering and LightGBM training, producing rolling-origin forecasts |
| `Paper/Experiments.R` | Forecasting methods, order-up-to inventory simulation, and accuracy and inventory metrics |

## Requirements

R (≥ 4.1) with the following packages:

```r
install.packages(c("forecast", "smooth", "tsutils", "lightgbm", "data.table", "dplyr",
                   "plyr", "tidyr", "readr", "vroom", "stringr", "lubridate", "zoo",
                   "reshape2", "tsfeatures", "Matrix", "MLmetrics", "moments",
                   "tseries", "fracdiff", "rpart", "foreach", "doParallel", "ggplot2"))
```

## Data

No data are included in this repository. The scripts use the publicly available M5 competition data (`sales_train_evaluation.csv` and `calendar.csv`), available from the [M5 Forecasting – Accuracy](https://www.kaggle.com/competitions/m5-forecasting-accuracy/data) competition on Kaggle.

## Usage

Set the working directory paths at the top of each script, then run:

```r
source("lgbm/lgbmForecasts.R")   # LightGBM forecasts
source("Paper/Experiments.R")    # inventory simulation for all methods
```

## Citation

If you wish to cite the paper, please use:

```bibtex
@article{theodorou2025forecast,
  title   = {Forecast accuracy and inventory performance: Insights on their relationship from the {M5} competition data},
  author  = {Theodorou, Evangelos and Spiliotis, Evangelos and Assimakopoulos, Vassilios},
  journal = {European Journal of Operational Research},
  volume  = {322},
  number  = {2},
  pages   = {414--426},
  year    = {2025},
  doi     = {10.1016/j.ejor.2024.12.033}
}
```

## License

This code is released under the [MIT License](LICENSE).
