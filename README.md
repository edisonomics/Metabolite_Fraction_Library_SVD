# Metabolite_Fraction_Library_SVD

An automated MATLAB pipeline for extracting purified component spectra and elution profiles from a 2D NMR metabolite fraction library (mFL) using peak tracing, Singular Value Decomposition with Alternating Least Squares (SVD-ALS), and hierarchical clustering.

## Table of Contents

- [Overview](#overview)
- [Availability & Reproducibility](#availability--reproducibility)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start (Test Data)](#quick-start-test-data)
- [Input Data Format & Outputs](#input-data-format--outputs)
- [Figure Map](#figure-map)
- [Archive (Original Manuscript Analysis)](#archive-original-manuscript-analysis)
- [Configuration & Optional Parameters](#configuration--optional-parameters)
- [Repository Structure](#repository-structure)
- [Reproducing Manuscript Results](#reproducing-manuscript-results)
- [Citation](#citation)
- [License](#license)
- [Support & Contact](#support--contact)

## Overview

`Metabolite_Fraction_Library_SVD` decomposes multi-fraction 1D NMR spectral libraries into a compact, non-redundant basis set of purified metabolite spectra and their corresponding chromatographic elution profiles across modular stages:

1. **Downsampling (Optional):** Linearly interpolates spectra onto a uniform chemical shift (`ppm`) grid (`downsampleData`).
2. **Preprocessing:** Replaces solvent and reference peak regions (e.g., DSS, residual water, formate) with synthetic baseline noise matched to a signal-free region and clips negative values (`preprocessFractions`).
3. **Peak Detection, Tracing & Deduplication:** Detects peaks per fraction (`findLibraryPeaks`), links peaks across consecutive fractions within a chemical shift tolerance (`tracePeaks`), and removes highly correlated duplicate traces sharing the same fractional span (`removeRedundantTraces`).
4. **SVD-ALS Component Extraction:** Computes per-node economy SVD, constructs gated covariance targets (`computeCovarianceTarget`), refines elution and spectral weights via alternating least squares (`refineWeightsALS`), and applies a correlation-modulated intensity filter (`filterExtractedSpectrum`) via `extractComponentsSVD`.
5. **Quality Control Filtering:** Removes extracted components whose Pearson correlation against their target covariance spectrum falls below a threshold (`filterComponentsQC`).
6. **Hierarchical Clustering:** Groups purified spectra by pairwise Pearson correlation distance (`clusterExtractedSpectra`).
7. **Correlation Cutoff Optimization:** Sweeps clustering correlation cutoffs, optimizes cluster elution assignments, ranks components by a weighted purity score, and selects the cutoff and component count (`optimal_k`) maximizing cumulative variance explained (%) (`optimizeClusterCutoff`).
8. **Visualization & Export:** Generates cumulative variance-explained scree plots (`plotVarianceExplained`) and stepped 3D library reconstructions (`plotSteppedReconstructions`, `displaySpectra3D`), and exports 300-DPI `.png` figures, CSV tables, consensus CSVs, and a `.mat` checkpoint.
9. **GISSMO Reference Library Matching (Downstream Workflow):** Loads and L2-normalizes a GISSMO NMR reference library (`gissmo_load_library`), matches extracted consensus spectra via global or split-half shift scanning (`gissmo_match_spectra`), exports the sorted match table, and optionally plots overlaid spectral matches (`gissmo_plot_matches`).
10. **Reusable Downstream Utilities (`src/`):** Includes `subtractComponentsSVD` (spectral library component subtraction) and `integrateComponentsSVD` (automated qNMR concentration quantification against an internal DSS reference). The manual component-indexed manuscript workflows using these helpers are archived in [`archive/`](archive/README.md).

## Availability & Reproducibility

- **Source Code:** Freely available to non-commercial users at [https://github.com/edisonomics/Metabolite_Fraction_Library_SVD](https://github.com/edisonomics/Metabolite_Fraction_Library_SVD).
- **Archived Release & Test Data (DOI):** [https://doi.org/10.5281/zenodo.24537564](https://doi.org/10.5281/zenodo.24537564)
- **Code Ocean Compute Capsule (DOI):** [https://doi.org/10.24433/CO.6916450.v1](https://doi.org/10.24433/CO.6916450.v1)

## Requirements

- **MATLAB:** R2019b or later (developed and tested with MATLAB R2023b on NMRbox running Ubuntu 24).
- **Required MathWorks Toolboxes:**
  - Optimization Toolbox (`optimproblem`, `fmincon`)
  - Signal Processing Toolbox (`findpeaks`)
  - Image Processing Toolbox (`imdilate`)
  - Statistics and Machine Learning Toolbox (`zscore`, `corr`, `pdist`, `linkage`, `dendrogram`, `cluster`)
- **Optional MathWorks Toolboxes:**
  - Parallel Computing Toolbox (only required if `use_parallel = true` in `extractComponentsSVD` or `optimizeClusterCutoff`)
- **Third-Party Dependencies:** None; the repository is self-contained.

## Installation

1. **Acquire the repository** by cloning with Git (or downloading and extracting the release archive):
   ```bash
   cd Metabolite_Fraction_Library_SVD
   ```
2. **Open MATLAB in the repository root directory** (the folder containing `setup.m` and `src/`). You can run `setup` manually, or run any workflow script directly—each workflow automatically invokes `setup` if `src/` is not yet on the MATLAB path.

## Quick Start (Test Data)

Default test datasets (`test_data/crassa.mat` and `test_data/GISSMO_Lib.mat`) are loaded automatically from `test_data/` unless you specify a different path in Section 2 (`USER CONFIGURATION`) of a workflow script.

1. Open MATLAB with your current working directory set to the repository root (`Metabolite_Fraction_Library_SVD`).
2. Run the primary SVD decomposition workflow (automatically runs `setup`, loads `test_data/crassa.mat`, and writes outputs to `results/`):
   ```matlab
   run(fullfile('workflows', 'mFL_svd_example_workflow.m'));
   ```
3. Run the GISSMO reference library matching workflow (uses workspace variables if present, or automatically loads `results/data/svd_results.mat` in a fresh MATLAB session):
   ```matlab
   run(fullfile('workflows', 'gissmo_matching_example_workflow.m'));
   ```

To write outputs to a custom folder anywhere on your computer instead of `<repo_root>/results/`, set `cfg.results_dir` in Section 2 of each workflow script (for example, `cfg.results_dir = '/Users/yourname/Desktop/mfl_results';` or `'C:\mfl_results'`).

## Input Data Format & Outputs

### Required Inputs

By default, `workflows/mFL_svd_example_workflow.m` loads `test_data/crassa.mat`. To analyze your own dataset, point `cfg.data_file` to a `.mat` file containing `X` and `ppm`, or set `cfg.data_file = ''` to use `X` and `ppm` already loaded in the workspace:

| Variable | Dimensions | Type | Description |
| --- | --- | --- | --- |
| `X` | `M x P` | `double` | Spectral intensity matrix (`M` fractions as rows, `P` chemical shift points as columns). |
| `ppm` | `1 x P` | `double` | Chemical shift vector in parts per million (ppm) corresponding to the columns of `X`. |

### Exported Output Directory (`results/`)

All exported files from the main workflows are written under `cfg.results_dir` (default: `<repo_root>/results/`). Figures are exported as **300-DPI `.png` files only** to keep runs fast and accessible without MATLAB.

```text
results/
├── figures/
│   ├── svd/
│   │   ├── preprocessed_library_3d.png
│   │   ├── detected_peaks_3d.png
│   │   ├── deduplicated_traces_3d.png
│   │   ├── qc_correlation_filter.png          # when cfg.show_qc_plot = true
│   │   ├── hierarchical_clustering.png
│   │   ├── variance_explained.png
│   │   ├── reconstruction_k05_3d.png          # per k in cfg.recon_k_steps
│   │   ├── reconstruction_k10_3d.png
│   │   ├── reconstruction_k25_3d.png
│   │   ├── reconstruction_full_3d.png
│   │   └── reconstruction_raw_3d.png
│   └── gissmo/
│       └── target_<ID>_<MatchName>.png        # when cfg.show_match_plots = true
├── tables/
│   ├── variance_explained.csv                 # from mFL_svd_example_workflow.m
│   └── gissmo_matches_sorted.csv              # from gissmo_matching_example_workflow.m
└── data/
    ├── consensus_spectra.csv                  # ppm column + Component_1..Component_C
    ├── consensus_elutions.csv                 # Fraction column + Component_1..Component_C
    └── svd_results.mat                        # Checkpoint for downstream workflows
```

### Primary Workspace Outputs

| Variable | Type | Workflow | Description |
| --- | --- | --- | --- |
| `best_state` | `struct` | `mFL_svd_example_workflow.m` | Optimal basis set (`cutoff`, `optimal_k`, `consensus_spectra`, `consensus_elutions`, `basis_spectra`, `basis_elutions`, `optimal_reconstructed_lib`, `purity_score`, `cumulative_var_explained`, `cluster_idx`). |
| `sweep_results` | `struct` | `mFL_svd_example_workflow.m` | Full cutoff sweep diagnostics (`cutoff_range`, `max_var_explained`, `optimal_ks`, `num_clusters`, `global_best_var_explained`, `best_cutoff`). |
| `variance_explained_pct` | `double` vector | `mFL_svd_example_workflow.m` | Cumulative variance explained (%) by intensity-sorted components up to `best_state.optimal_k`. |
| `sorted_gissmo_results` | `table` | `gissmo_matching_example_workflow.m` | GISSMO spectral match results ranked by descending cosine score. |

## Figure Map

The table below maps the GISSMO matching, subtraction, and integration figures to their generating scripts and exported PNG paths:

| Analysis | Workflow Script | Exported PNG File(s) | Description |
| --- | --- | --- | --- |
| GISSMO Matching | `workflows/gissmo_matching_example_workflow.m` | `results/figures/gissmo/target_<ID>_<MatchName>.png` | Overlaid consensus spectrum vs. top GISSMO reference match per selected component (enabled when `cfg.show_match_plots = true`). |
| Library Subtraction *(Archived)* | `archive/code/svd_subtraction_example_workflow.m` | `archive/result/figures/subtraction/residual_library_3d.png` | 3D stacked plot of the residual spectral library after subtracting the full `optimal_k` basis set across all fractions using `window_mask` (`cfg.show_residual_3d = true`). |
| Library Subtraction *(Archived)* | `archive/code/svd_subtraction_example_workflow.m` | `archive/result/figures/subtraction/original_vs_residual_2d.png` | 2D overlay of original (black) vs. residual (red) spectral library (enabled when `cfg.show_overlay_2d = true`). |
| qNMR Integration *(Archived)* | `archive/code/svd_integration_example_workflow.m` | `archive/result/figures/integration/metabolite_concentrations.png` | Publication-style log-scale grouped bar chart comparing manual vs. SVD metabolite concentrations (mM) with $R^2$ (`cfg.show_conc_plot = true`). |

## Archive (Original Manuscript Analysis)

Original monolithic scripts, the saved manuscript decomposition (`archive/result/output.mat`), and the manual component-indexed **subtraction** and **qNMR integration** workflows are preserved in [`archive/`](archive/README.md). These archived workflows load `archive/result/output.mat` and `test_data/crassa.mat` directly and export their figures and `metabolite_concentrations.csv` under `archive/result/`. See [`archive/README.md`](archive/README.md) for full instructions.

## Configuration & Optional Parameters

All user-configurable parameters in `workflows/mFL_svd_example_workflow.m` (and their corresponding function defaults in `src/`) are documented below.

### Stage 1–2: Input, Output, Downsampling & Preprocessing

| Parameter (`cfg.*`) | Default | Function | Purpose |
| --- | --- | --- | --- |
| `data_file` | `fullfile('test_data', 'crassa.mat')` | Workflow | Path to a `.mat` file containing `X` and `ppm`. Set to `''` to use variables already loaded in the workspace. |
| `results_dir` | `'results'` | All Workflows | Output root directory (relative or absolute path anywhere on your computer). |
| `save_outputs` | `true` | All Workflows | Export 300-DPI PNG figures, CSV tables, and data outputs to `cfg.results_dir`. |
| `do_downsample` | `false` | Workflow | Set `true` to linearly interpolate spectra to `num_points` via `downsampleData`. |
| `num_points` | `8192` | `downsampleData` | Target number of `ppm` points when `do_downsample = true` (must be strictly less than `length(ppm)`). |
| `noise_region` | `[0.2, 0.6]` | `preprocessFractions` | `[min_ppm, max_ppm]` signal-free baseline region used to estimate synthetic noise mean and standard deviation. |
| `remove_regions` | `[-0.244, 0.331; 4.73, 4.82; 8.40644, 8.48747]` | `preprocessFractions` | `K x 2` matrix of `[min_ppm, max_ppm]` regions to replace with synthetic baseline noise (e.g., DSS, residual water, formate). |
| `rng_seed` | `42` | Workflow | Fixed random seed (`rng(cfg.rng_seed, 'twister')`) for reproducible synthetic baseline noise in `preprocessFractions`. Leave empty (`[]`) to skip seeding. |
| `show_preprocess_plot` | `true` | `preprocessFractions` | Display 3D stacked plot of the preprocessed library via `displaySpectra3D`. |

### Stage 3: Peak Detection, Tracing & Redundancy Removal

| Parameter (`cfg.*`) | Default | Function | Purpose |
| --- | --- | --- | --- |
| `noise_threshold` | `0.05` | `findLibraryPeaks` | Minimum peak height (`MinPeakHeight`) passed to `findpeaks`. Choose ~2.5x baseline noise to avoid noise streaks. |
| `ppm_tol` | `0.003` | `tracePeaks` | Maximum chemical shift (`ppm`) difference allowed when linking a peak across adjacent fractions. |
| `min_trace_length` | `3` | `tracePeaks` | Minimum number of points required to retain a trace (including zero-padded boundary anchors). |
| `trace_corr_threshold` | `0.99` | `removeRedundantTraces` | Pearson correlation cutoff above which lower-intensity traces sharing the exact same fractional span are removed. |
| `show_peaks_plot` | `true` | `plotPeaksAndTraces` | Display 3D library overlay with detected peak markers. |
| `show_traces_plot` | `true` | `plotPeaksAndTraces` | Display 3D library overlay with deduplicated peak traces. |

### Stage 4–6: SVD-ALS Extraction, Quality Control & Clustering

| Parameter (`cfg.*`) | Default | Function | Purpose |
| --- | --- | --- | --- |
| `svd_opts.use_parallel` | `false` | `extractComponentsSVD` | Set `true` to run per-node SVD-ALS across workers with `parfor` (requires Parallel Computing Toolbox). |
| `svd_opts.verbose` | `true` | `extractComponentsSVD` | Print start and completion status messages to the Command Window. |
| `qc_corr_threshold` | `0.2` | `filterComponentsQC` | Minimum Pearson correlation between an extracted `clean_spectrum` and its target covariance spectrum required to retain a component. |
| `qc_target_line` | `0.8` | `filterComponentsQC` | Horizontal reference line drawn on the QC correlation scatter plot. |
| `show_qc_plot` | `false` *(fn default: `true`)* | `filterComponentsQC` | Display QC scatter and histogram diagnostic figure. |
| `linkage_method` | `'complete'` | `clusterExtractedSpectra` | Hierarchical linkage method passed to MATLAB `linkage` (e.g., `'complete'`, `'average'`, `'single'`). |
| `show_cluster_plot` | `true` | `clusterExtractedSpectra` | Display paired dendrogram and sorted pairwise correlation heatmap. |

### Stage 7: Correlation Cutoff Optimization

| Parameter (`cfg.sweep_opts.*`) | Default | Function | Purpose |
| --- | --- | --- | --- |
| `cutoff_range` | `0.5:0.02:0.95` | `optimizeClusterCutoff` | Vector of Pearson correlation cutoffs evaluated during hierarchical clustering optimization. |
| `weight_covar` | `0.8` | `optimizeClusterCutoff` | Purity score weight assigned to target covariance match (`spec_correlations`). |
| `weight_sim` | `0.0` | `optimizeClusterCutoff` | Purity score weight assigned to intra-cluster spectral similarity. |
| `weight_intensity` | `0.2` | `optimizeClusterCutoff` | Purity score weight assigned to initial peak intensity. |
| `max_iters` | `20` | `optimizeClusterCutoff` | Maximum iterations for alternating cluster elution member selection. |
| `snap_spectra_to_winner` | `false` | `optimizeClusterCutoff` | If `true`, replaces each cluster's median consensus spectrum with the winning member's spectrum. |
| `use_parallel` | `false` | `optimizeClusterCutoff` | Set `true` to evaluate cutoffs in parallel with `parfor`. |
| `verbose` | `true` | `optimizeClusterCutoff` | Print per-cutoff progress and best-cutoff summary. |

### Stage 8: Final Visualizations

| Parameter (`cfg.*`) | Default | Function | Purpose |
| --- | --- | --- | --- |
| `show_variance_plot` | `true` | `plotVarianceExplained` | Render scree plot of cumulative variance explained (%) up to `optimal_k`. |
| `variance_k_markers` | `[5, 10, 25]` | `plotVarianceExplained` | Component counts at which to draw vertical dotted reference lines on the scree plot. |
| `show_stepped_recon` | `true` | `plotSteppedReconstructions` | Render stepped 3D reconstructions at increasing component counts. |
| `recon_k_steps` | `[5, 10, 25]` | `plotSteppedReconstructions` | Component counts plotted as stepped 3D reconstructions. |
| `recon_include_full` | `true` | `plotSteppedReconstructions` | Also plot the full `optimal_k` intensity-sorted reconstruction. |
| `recon_include_raw` | `true` | `plotSteppedReconstructions` | Also plot the raw preprocessed library for side-by-side comparison. |

### Standalone Visualization Optional Arguments (`displaySpectra3D`)

`displaySpectra3D(X, ppm, isbackground)` accepts one optional argument: `isbackground` (`logical`, default `false`), which colors the last spectrum black when set to `true`. The function is completely self-contained with no external third-party function dependencies.

### GISSMO Reference Library Matching (`gissmo_*`)

An optional GISSMO NMR reference library (`test_data/GISSMO_Lib.mat`, containing `lib_combined`, `ppm_lib`, and `combined_names`) can be loaded with `gissmo_load_library(mat_path)` and matched against extracted consensus spectra using `gissmo_match_spectra(spec_matrix, ppm_target, gissmo_source, ...)` and `gissmo_plot_matches(match_out, ...)`. Their optional Name-Value parameters are:

| Option | Default | Function | Purpose |
| --- | --- | --- | --- |
| `'MinFrac'` | `0.0` | `gissmo_match_spectra` | Peak removal threshold as a fraction of the maximum peak intensity (`0` to `1`). |
| `'ShiftMax'` | `50` | `gissmo_match_spectra` | Maximum chemical shift point offset tested in either direction. |
| `'ShiftStep'` | `5` | `gissmo_match_spectra` | Point step size between tested shift offsets. |
| `'SearchMode'` | `"global"` | `gissmo_match_spectra` | Shift alignment mode: `"global"` (uniform shift) or `"split"` (independent left/right halves). |
| `'PlotMatches'` | `false` | `gissmo_match_spectra` | Automatically invoke `gissmo_plot_matches` after matching. |
| `'Verbose'` | `true` | `gissmo_match_spectra` | Print formatted match summary tables to the Command Window. |
| `'TargetIDs'` | `[]` (all) | `gissmo_plot_matches` | Vector of specific `Target_ID` values to plot. |
| `'MinScore'` | `-Inf` | `gissmo_plot_matches` | Plot only matches with cosine score `>= MinScore`. |
| `'SortByScore'` | `false` | `gissmo_plot_matches` | Plot matches in descending score order instead of `Target_ID` order. |
| `'SaveDir'` | `""` *(workflow defaults to `<results_dir>/figures/gissmo`)* | `gissmo_plot_matches` | Output directory to save 300-DPI PNG overlays. |

### Archived Subtraction & Integration Workflows

- **`archive/code/svd_subtraction_example_workflow.m`:** Loads `archive/result/output.mat` and `test_data/crassa.mat`, subtracts the full `optimal_k` ranked basis set across all fractions using `'window_mask'` (`+/- 0.03 ppm`) via `subtractComponentsSVD`, and saves 300-DPI PNG figures under `archive/result/figures/subtraction/`.
- **`archive/code/svd_integration_example_workflow.m`:** Loads `archive/result/output.mat` and `test_data/crassa.mat`, runs automated SVD integration via `integrateComponentsSVD` alongside the manual peak-bound benchmark ($R^2$), and exports `archive/result/tables/metabolite_concentrations.csv` and `archive/result/figures/integration/metabolite_concentrations.png`.

## Repository Structure

```text
Metabolite_Fraction_Library_SVD/
├── setup.m                                    # Adds src/ and workflows/ to the MATLAB path (excludes archive/)
├── README.md                                  # Project documentation
├── LICENSE.txt                                # License file (PolyForm Noncommercial 1.0.0)
├── NOTICE.txt                                 # Required copyright and attribution notice
├── test_data/
│   ├── crassa.mat                             # Example NMR metabolite fraction library dataset
│   └── GISSMO_Lib.mat                         # GISSMO NMR reference spectral library
├── workflows/
│   ├── mFL_svd_example_workflow.m             # End-to-end configurable SVD decomposition pipeline
│   └── gissmo_matching_example_workflow.m     # GISSMO reference library matching workflow
├── archive/
│   ├── README.md                              # Documentation for archived manuscript scripts & outputs
│   ├── code/
│   │   ├── SVD_mFL_Analysis.m                 # Original monolithic SVD decomposition script
│   │   ├── GISSMO_Matching.m                  # Original monolithic GISSMO matching script
│   │   ├── svd_subtraction_example_workflow.m # Archived manuscript subtraction workflow
│   │   └── svd_integration_example_workflow.m # Archived manuscript qNMR integration workflow
│   └── result/
│       └── output.mat                         # Archived manuscript best_state and ppm
└── src/
    ├── downsampleData.m                       # Linear interpolation onto a uniform ppm grid
    ├── preprocessFractions.m                  # Solvent/reference replacement with noise
    ├── findLibraryPeaks.m                     # Per-fraction peak detection
    ├── tracePeaks.m                           # Cross-fraction peak tracing
    ├── removeRedundantTraces.m                # Same-span correlated trace deduplication
    ├── extractComponentsSVD.m                 # Per-node SVD-ALS extraction driver
    ├── computeCovarianceTarget.m              # Compute target covariance spectrum
    ├── refineWeightsALS.m                     # Alternating least squares weight refinement
    ├── filterExtractedSpectrum.m              # Correlation-modulated thresholding & smoothing
    ├── filterComponentsQC.m                   # Target-correlation quality control filtering
    ├── clusterExtractedSpectra.m              # Hierarchical clustering of purified spectra
    ├── optimizeClusterCutoff.m                # Correlation cutoff sweep & basis selection
    ├── subtractComponentsSVD.m                # Subtract selected consensus components from library
    ├── integrateComponentsSVD.m               # Automated qNMR integration against DSS reference
    ├── gissmo_load_library.m                  # Load and L2-normalize GISSMO reference library
    ├── gissmo_match_spectra.m                 # Match consensus spectra against GISSMO library
    ├── gissmo_plot_matches.m                  # Plot consensus vs. matched GISSMO spectra
    ├── plotPeaksAndTraces.m                   # 3D peak/trace overlay visualization
    ├── plotVarianceExplained.m                # Cumulative variance explained scree plot
    ├── plotSteppedReconstructions.m           # Stepped 3D reconstruction comparisons
    ├── displaySpectra3D.m                     # 3D stacked spectral visualization utility
    └── private/
        ├── optimizeElution.m                  # Time-domain elution weight optimization
        ├── optimizeSpectra.m                  # Frequency-domain spectral weight optimization
        ├── optimizeClusterElutions.m          # Cluster elution member selection
        └── rankAndEvaluateClusters.m          # Purity ranking & cumulative variance evaluation
```

## Reproducing Manuscript Results

To reproduce the automated decomposition, GISSMO matching, and diagnostic figures on the included test dataset (`test_data/crassa.mat`):

1. Open MATLAB in the repository root directory (`Metabolite_Fraction_Library_SVD`).
2. Run `workflows/mFL_svd_example_workflow.m` with the default `cfg` parameters (`cfg.rng_seed = 42`, `cfg.svd_opts.use_parallel = false`, `cfg.sweep_opts.use_parallel = false`).
3. Run `workflows/gissmo_matching_example_workflow.m` (either in the same session or in a fresh session using `results/data/svd_results.mat`).
4. Inspect the exported 300-DPI `.png` figures in `results/figures/`, CSV tables in `results/tables/`, and consensus spectra/elutions CSVs in `results/data/`.
5. To reproduce the original manuscript component subtraction and 23-metabolite qNMR integration comparison against `archive/result/output.mat`, follow [`archive/README.md`](archive/README.md).

- **Expected Runtime:** Approximately **1 hour** in serial mode (`use_parallel = false`, default). Runtime can be reduced by running in parallel mode (`cfg.svd_opts.use_parallel = true` and `cfg.sweep_opts.use_parallel = true`; requires the Parallel Computing Toolbox).
- **Browser-Based Reproducibility:** Available via the [Code Ocean Compute Capsule](https://doi.org/10.24433/CO.6916450.v1).

## Citation

If you use this software in your research, please cite our manuscript and the archived software release:

- **Manuscript (bioRxiv):** [https://doi.org/10.64898/2026.10.30.763183](https://doi.org/10.64898/2026.10.30.763183)
- **Archived Software Release (Zenodo):** [https://doi.org/10.5281/zenodo.24537564](https://doi.org/10.5281/zenodo.24537564)


## License

This software is licensed under the **PolyForm Noncommercial License 1.0.0** (`PolyForm-Noncommercial-1.0.0`). See [LICENSE.txt](LICENSE.txt) and [NOTICE.txt](NOTICE.txt) for full terms. Non-commercial academic and research use is freely permitted without request. For commercial licensing inquiries, please contact the Edison Lab via the repository maintainers.


## Support & Contact

For bug reports, questions, or feature requests, please open an issue on the [GitHub Issues page](https://github.com/edisonomics/Metabolite_Fraction_Library_SVD/issues).
