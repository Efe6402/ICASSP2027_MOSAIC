# MOSAIC S&P 500 experiment

This package learns a multilayer graph for 100 S&P 500 stocks observed through cumulative-return signals at 1, 20, 45, 70, and 120 trading-day horizons.

## Main experiment

Open MATLAB in this directory and run:

```matlab
RUN_SP500_SUBMISSION
```

The experiment uses five endpoint-offset configurations. Offsets 1 and 3 are estimation configurations, offsets 2 and 4 are their corresponding validation configurations, and offset 0 is excluded from model selection. Modal ranks 2--5 and 128 hyperparameter profiles are evaluated using validation-configuration smoothness. The selected rank and top-1 hyperparameter configuration are then fitted to the prepared offset-0 dataset.

Results are written to `results/run_<timestamp>`. The paper-reported top-1 configuration, fitted model, numerical tables, and principal figures are available in `selected_top1` before rerunning the experiment.

## Data preparation

The stock-price inputs, sector labels, and selected ticker list are under `data/stocks`. The exact offset-0 dataset used for final estimation is under `data/prepared_offset0`.

To regenerate the dated price CSV from the included RDS object, run:

```bash
Rscript convert_rds_to_csv_with_dates.R
```

To construct a new offset-0 multihorizon dataset from the packaged stock data, run:

```matlab
RUN_DATA_PREPARATION
```

The main experiment reconstructs all five offset configurations from the same packaged price data and verifies that its generated offset-0 signals agree numerically with the included prepared dataset before fitting the final model. Sector labels are used only for interpretation after graph learning.

Run `VERIFY_PACKAGE.m` for package and prepared-data checks.
