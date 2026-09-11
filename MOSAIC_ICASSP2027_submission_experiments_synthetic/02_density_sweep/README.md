# Modal-graph density sweep

Run `RUN_FULL_DENSITY_SWEEP.m`.

The experiment uses \(n=16\) physical nodes, \(K=4\) layers, two binary modal graphs, and \(p=5000\) graph signals. Starting from modal edge counts \((10,9)\), two physical edges are added to each mode at every level until the counts reach \((30,29)\). The supports are nested, and new edges are allocated between the designated eight-node community and the background according to a fixed within/background odds ratio.

The modal contribution profiles and same-node copy coefficients remain fixed. The changing modal graphs define a new supra-adjacency and signal distribution at each density level. Ground-truth quantities are used only for data generation and recovery evaluation; all model variables are estimated from the observed signals.

MOSAIC, PGL, and MXGL are evaluated using the included \(p=5000\) component-F1-selected configurations and component-specific validation thresholds. Hyperparameters, thresholds, and candidate ranks are not recalibrated across density levels. The output reports test AUC and F1 for the supra-adjacency and each structural component.

`RUN_QUICK_TEST.m` provides a reduced end-to-end check.
