# COIL-20 six-view MOSAIC experiment

This package contains the six-view COIL-20 experiment used for the paper. It performs auxiliary cross-configuration selection over modal ranks `2:6` and a 25-profile focused hyperparameter bank, then fits the selected rank/profile once on the phase-0 target data.

## Run

Open MATLAB in this folder and run:

```matlab
RUN_COIL20_SUBMISSION
```

The run uses the packaged 32-by-32 COIL-20 signals. Phase 5, 25, and 45 datasets are fitted during selection and evaluated on phases 10, 30, and 50, respectively. Phase 0 is used only for the final fit. Every fit uses the complete 1024-dimensional image vector.

The selection score is the mean held-out smoothness of the mass-normalized learned supra-graph. For each rank, the lowest-scoring profile is retained, and the modal rank is chosen by the geometric endpoint-chord elbow. The final target fit uses the resulting top-1 profile; target data do not rerank candidates.

The selected paper configuration is recorded in [`selected_top1/SELECTED_TOP1.csv`](selected_top1/SELECTED_TOP1.csv), and its fitted model is included in the same folder. `GENERATE_PAPER_FIGURES.m` regenerates the four separate paper figures without rerunning model selection.

## Main outputs

- `results/selection/profile_rank_validation.csv`
- `results/selection/rank_selection.csv`
- `results/selection/SELECTED_TOP1.csv`
- `results/final/FINAL_TOP1_FIT.csv`
- `results/final/top01_r3_profile2508_G08_beta_c_010.mat`
- `paper_figures/01_modal_view_profiles.png`
- `paper_figures/02_within_view_graph.png`
- `paper_figures/03_crossnode_crossview_pairs.png`
- `paper_figures/04_induced_view_graph.png`

The `reference/coil.png` file is the paper montage assembled from these four panels.
