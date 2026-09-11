# MOSAIC synthetic experiments

This package contains the synthetic experiments used to evaluate MOSAIC against product graph learning (PGL) and multiplex graph learning (MXGL). Each experiment is self-contained: paths are resolved relative to its own directory, and generated data, fitted models, tables, and figures are written locally within that subpackage.

## Requirements

- MATLAB R2025b or a compatible release
- Statistics and Machine Learning Toolbox
- No external dataset download is required

## Experiments

1. `01_sample_size_and_model_selection` performs model selection across the graph-signal sample size \(p\), then evaluates the configurations selected at \(p=5000\) across the full sample-size grid.
2. `02_density_sweep` increases the density of the binary modal graphs through nested edge additions while holding \(p=5000\), the modal contribution profiles, and the copy coefficients fixed.
3. `03_view_profile_rho_sweep` varies the modal mixing parameter \(\rho\) while holding the modal graphs, copy coefficients, and \(p=5000\) fixed.
4. `04_noise_sweep` varies both the signal-to-noise ratio (SNR) and the degree of view-wise Gaussian noise heteroscedasticity at \(p=5000\).

## Running the experiments

From the package root, run:

```matlab
run('01_sample_size_and_model_selection/RUN_COMPLETE_SAMPLE_EXPERIMENT.m')
run('02_density_sweep/RUN_FULL_DENSITY_SWEEP.m')
run('03_view_profile_rho_sweep/RUN_FULL_RHO_SWEEP.m')
run('04_noise_sweep/RUN_FULL_NOISE_SWEEP.m')
```

The first experiment includes the full hyperparameter search and is the most computationally intensive. Each subpackage also provides `RUN_QUICK_TEST.m` for a small end-to-end installation check.

`CREATE_CONFERENCE_TOP1_SUPRA_F1.m` reads the latest completed runs from experiments 01--03 and saves their combined rank-one supra-adjacency F1 figure under `conference_figures`. Experiment 04 generates its noise-sweep figures within its own result directory.

## Selection and evaluation

Hyperparameters and graph-support thresholds are selected using validation realizations only. The reported candidate ranking is based on the mean validation F1 averaged over the within-view, cross-node/cross-view, and same-node copy components. The highest-ranked configurations and their component-specific validation thresholds are then evaluated on independent test realizations. AUC is computed from continuous edge scores, whereas F1 is computed from thresholded supports.

The density, modal-profile, and noise experiments use the included \(p=5000\) selections and thresholds without condition-specific retuning.
