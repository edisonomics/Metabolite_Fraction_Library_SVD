# Archive — Original Manuscript Scripts & Manual Downstream Analyses

This directory preserves the original manuscript scripts, the saved decomposition state (`archive/result/output.mat`) from the manuscript run, and the downstream **integration** and **component subtraction** workflows that depend on manual metabolite-to-component-index assignments tied to `archive/result/output.mat`.

The reusable core helper functions (`src/integrateComponentsSVD.m` and `src/subtractComponentsSVD.m`) remain in `src/` so users can apply them to any new dataset or decomposition output.

## Contents

```text
archive/
├── README.md                                  # This file
├── code/
│   ├── SVD_mFL_Analysis.m                     # Original monolithic SVD decomposition script
│   ├── GISSMO_Matching.m                      # Original monolithic GISSMO matching script
│   ├── svd_subtraction_example_workflow.m     # Manuscript full optimal_k basis subtraction workflow (window_mask)
│   └── svd_integration_example_workflow.m     # Manuscript 23-metabolite qNMR integration & manual comparison
└── result/
    └── output.mat                             # Archived manuscript best_state (1x1 struct) and ppm (1x32768)
```

## Why Integration and Subtraction Are Archived

The manuscript qNMR integration script references specific component row indices (`cfg.component_idx`) matched to the cluster numbering in `archive/result/output.mat`, and the archived subtraction workflow reproduces the manuscript residual library by subtracting the full `optimal_k` ranked basis set (`best_state.rank_idx(1:best_state.optimal_k)`) across all fractions using `'window_mask'` (`+/- 0.03 ppm`) directly against `archive/result/output.mat`. Running both archived workflows against `archive/result/output.mat` and `test_data/crassa.mat` reproduces the exact manuscript quantification (including the manual integration benchmark and $R^2$) and residual library figures.

## Running the Archived Workflows

Open MATLAB in the **repository root directory** (`Metabolite_Fraction_Library_SVD`) and run either script directly (each automatically invokes `setup` to add `src/` to the path):

```matlab
% Reproduce manuscript full optimal_k basis subtraction (window_mask, all fractions)
run(fullfile('archive', 'code', 'svd_subtraction_example_workflow.m'));

% Reproduce manuscript 23-metabolite qNMR integration (SVD vs. manual + R^2)
run(fullfile('archive', 'code', 'svd_integration_example_workflow.m'));
```

## Exported Outputs (`archive/result/`)

By default, both archived workflows write outputs under `cfg.results_dir = fullfile('archive', 'result')` (configurable to any custom directory path):

| Analysis | Workflow Script | Exported File(s) |
| --- | --- | --- |
| Library Subtraction | `archive/code/svd_subtraction_example_workflow.m` | `archive/result/figures/subtraction/residual_library_3d.png`, `archive/result/figures/subtraction/original_vs_residual_2d.png` (when `cfg.show_overlay_2d = true`) |
| qNMR Integration | `archive/code/svd_integration_example_workflow.m` | `archive/result/figures/integration/metabolite_concentrations.png`, `archive/result/tables/metabolite_concentrations.csv` |
