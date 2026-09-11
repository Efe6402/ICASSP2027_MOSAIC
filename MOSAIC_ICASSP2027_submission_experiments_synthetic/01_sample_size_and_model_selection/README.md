# Sample-size sweep and model selection

Run `RUN_COMPLETE_SAMPLE_EXPERIMENT.m` to execute the full experiment.

The experiment evaluates

```text
p = [20, 100, 300, 1000, 3000, 5000, 10000]
```

with ten validation and ten independent test realizations at each sample size. Within each realization, smaller signal matrices are prefixes of the draw generated at the largest \(p\), providing paired comparisons across the sample-size grid.

For each method and each \(p\), a coarse-to-fine hyperparameter search is performed on validation data. The final ranking uses the mean validation F1 averaged over the within-view, cross-node/cross-view, and same-node copy components. Separate relative support thresholds are selected on validation data for the supra-adjacency and its three components, after which the selected configurations are evaluated on the test realizations.

The complete entry point also evaluates the configurations and thresholds selected at \(p=5000\) across all sample sizes. To run only this fixed-configuration transfer stage using the included \(p=5000\) selections, use `RUN_P5000_FIXED_TRANSFER.m`.

`RUN_QUICK_TEST.m` provides a reduced end-to-end check.
