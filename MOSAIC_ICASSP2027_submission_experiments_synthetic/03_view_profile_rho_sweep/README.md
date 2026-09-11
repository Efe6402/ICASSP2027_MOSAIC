# Modal-profile mixing sweep

Run `RUN_FULL_RHO_SWEEP.m`.

The experiment uses \(n=16\) physical nodes, \(K=4\) layers, two modes, and \(p=5000\) graph signals. It evaluates

```text
rho = [5/102, 0.075, 0.125, 0.20, 0.35, 0.50].
```

At each value of \(\rho\), the two unnormalized layer-wise contribution profiles are

```text
[1, 1-rho,   rho, 0] and [0, rho, 1-rho, 1].
```

The corresponding column-simplex modal profiles are obtained by taking elementwise square roots and normalizing each profile to sum to one. The physical modal graphs and same-node copy coefficients remain fixed, while the modal profiles, supra-adjacency, and signal distribution vary with \(\rho\). Signal innovations are paired across the sweep within each split and realization.

MOSAIC, PGL, and MXGL are evaluated using the included \(p=5000\) component-F1-selected configurations and component-specific validation thresholds. Hyperparameters, thresholds, and candidate ranks are not recalibrated across \(\rho\). The run also saves the ground-truth modal profiles and structural supports for each condition.

`RUN_QUICK_TEST.m` provides a reduced end-to-end check.
