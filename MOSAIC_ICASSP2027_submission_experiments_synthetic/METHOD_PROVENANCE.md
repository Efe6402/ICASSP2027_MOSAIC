# Method provenance

The package evaluates three topology-learning methods on a common signal bank and a common supra-adjacency evaluation protocol.

- **MOSAIC:** `solve_mosaic_crossview.m` implements the projected cyclic block updates for modal physical graphs, copy coefficients, and column-simplex view profiles.
- **PGL2021:** `Learn_PGL.m`, `PGL_solver.m`, and their linear-algebra helpers implement sparse product-graph factor learning from the physical-node and view-domain covariance summaries.
- **ZW2024:** `solve_zhang_wai_2024_paper.m` implements multiplex graph learning with inter-layer coupling. Its effective interaction matrix is used as the supra-adjacency score.

The SHA-256 inventory in `SOLVER_SHA256.csv` makes the solver snapshot auditable. Each subpackage contains a complete local solver tree, and setup functions add only files within that subpackage to the MATLAB path.

The embedded configuration tables contain the three candidates selected at `p=5000` by mean validation component F1. Their four component-specific thresholds were selected on validation realizations and are applied to held-out test realizations in the transfer experiments.
