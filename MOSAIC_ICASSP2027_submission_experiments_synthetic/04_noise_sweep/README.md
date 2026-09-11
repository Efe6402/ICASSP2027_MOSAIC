# View-wise Gaussian noise sweep

Run `RUN_FULL_NOISE_SWEEP.m`.

The experiment uses the \(p=5000\) synthetic setting and evaluates the clean signals together with prescribed SNR levels

```text
SNR (dB) = [clean, 20, 15, 12, 10, 8, 6, 4, 2, 0]
```

and heteroscedasticity levels

```text
kappa = [0, 0.2, 0.4, 0.6, 0.8, 1].
```

Zero-mean Gaussian noise is added directly to the graph signals. Its view-variance profile is

```text
v(kappa) = 1 + kappa * ([0.4, 0.8, 1.2, 1.6] - 1).
```

Thus, \(\kappa=0\) gives equal noise variance in all four views, whereas \(\kappa=1\) gives the full view-wise heteroscedastic profile. Every profile has unit mean variance, so SNR controls the total noise intensity while \(\kappa\) controls how that intensity is distributed across views. The profile is cyclically rotated across the ten test realizations, and common clean signals and Gaussian innovations provide paired comparisons across conditions.

MOSAIC, PGL, and MXGL are evaluated using the included \(p=5000\) component-F1-selected configurations and component-specific validation thresholds. No hyperparameter search, threshold recalibration, or reranking is performed during the noise sweep. The output reports supra-adjacency and component-level AUC and F1 across all \((\kappa,\mathrm{SNR})\) conditions.

`RUN_QUICK_TEST.m` provides a reduced end-to-end check.
