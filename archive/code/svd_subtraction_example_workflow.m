%% Archived Spectral Library Component Subtraction Workflow
% =========================================================================
% Subtracts the full optimized-k ranked SVD basis set from the preprocessed
% spectral library using the archived manuscript decomposition
% (archive/result/output.mat), calculates Frobenius-norm variance
% removed/remaining metrics across all fractions, and exports 300-DPI PNG
% residual visualizations under cfg.results_dir (default: archive/result).
%
% Pipeline Stages:
%   1. Load archived best_state & ppm (archive/result/output.mat) and
%      reconstruct preprocessed library (lib) from test_data/crassa.mat
%      if not already present in the workspace
%   2. Subtract the full optimal_k ranked basis set via subtractComponentsSVD
%   3. Report Frobenius-norm variance removed and remaining metrics
%   4. Optionally display and export 300-DPI PNG 3D residual and 2D overlay
%      plots under archive/result/figures/subtraction/
%
% External Functions Required:
%   subtractComponentsSVD, preprocessFractions,
%   displaySpectra3D (when cfg.show_residual_3d is true)
%
% Required MATLAB Toolboxes:
%   None (compatible with MATLAB R2019b or later)
% =========================================================================

%% 1. VERIFY PATH SETUP
% Run from the repository root directory (the folder containing setup.m and src/).
if isempty(which('subtractComponentsSVD'))
    if isfile('setup.m')
        setup;
    else
        error('svd_subtraction_example_workflow:PathNotSet', ...
            'Run this workflow from the repository root directory (the folder containing setup.m and src/).');
    end
else
    disp('Repository is on MATLAB path and archived SVD subtraction workflow is ready.');
end

%% 2. USER CONFIGURATION
% Edit the parameters below to match your dataset and subtraction settings.

% --- Output Directory & Exports ---
% Defaults to archive/result; set to any relative or absolute path if desired.
cfg.results_dir  = fullfile('archive', 'result');
cfg.save_outputs = true;      % Export 300-DPI PNG subtraction figures

% --- Input Data ---
% Archived manuscript best_state & ppm, paired with raw crassa.mat.
cfg.svd_results_file = fullfile('archive', 'result', 'output.mat');
cfg.raw_data_file    = fullfile('test_data', 'crassa.mat');

% Preprocessing parameters used to build 'lib' from crassa.mat when 'lib'
% is not present in output.mat or the workspace.
cfg.noise_region   = [0.2, 0.6];
cfg.remove_regions = [ ...
    -0.244,   0.331;    % DSS / reference peak
     4.73,    4.82;     % Residual water
     8.40644, 8.48747   % Formate
];
cfg.rng_seed = 42;

% --- Components to Subtract ---
% Leave empty ([]) to subtract the full optimized-k ranked basis set:
% best_state.rank_idx(1:best_state.optimal_k). Supply a custom vector of
% row indices in best_state.consensus_spectra to override.
cfg.component_idx = [];

% --- Subtraction Mode ---
% Options: 'zero_mask', 'window_mask', 'linear', 'linear_clipped'
cfg.method     = 'window_mask';
cfg.ppm_window = 0.03;   % Half-width in ppm (used when method = 'window_mask')

% --- Variance Metric Fraction Window ---
% Empty ([]) evaluates variance metrics across all fractions in lib.
cfg.frac_idx = [];

% --- Visualization & Console Output ---
cfg.show_residual_3d = true;   % Display 3D residual library via displaySpectra3D
cfg.show_overlay_2d  = false;  % Display 2D overlay of original lib (black) vs. residual (red)
cfg.overlay_ylim     = [-2, 80];
cfg.verbose          = true;   % Print variance metrics to Command Window

%% 3. LOAD & VALIDATE INPUTS
total_timer = tic;

if ~isfile(cfg.svd_results_file)
    error('svd_subtraction_example_workflow:FileNotFound', ...
        'Archived SVD results file not found: %s', cfg.svd_results_file);
end
fprintf('Loading archived SVD results from: %s\n', cfg.svd_results_file);
loaded_svd = load(cfg.svd_results_file);

if ~isfield(loaded_svd, 'best_state')
    error('svd_subtraction_example_workflow:MissingBestState', ...
        'MAT file must contain variable ''best_state''.');
end
best_state = loaded_svd.best_state;

if isfield(loaded_svd, 'ppm')
    ppm_target = loaded_svd.ppm;
elseif isfield(loaded_svd, 'ppm_work')
    ppm_target = loaded_svd.ppm_work;
else
    error('svd_subtraction_example_workflow:MissingPPM', ...
        'MAT file must contain variable ''ppm'' or ''ppm_work''.');
end

if isfield(loaded_svd, 'lib')
    lib = loaded_svd.lib;
else
    if ~isfile(cfg.raw_data_file)
        error('svd_subtraction_example_workflow:MissingRawData', ...
            'Preprocessed ''lib'' not in %s and raw file not found: %s', ...
            cfg.svd_results_file, cfg.raw_data_file);
    end
    fprintf('Reconstructing preprocessed library ''lib'' from: %s\n', cfg.raw_data_file);
    loaded_raw = load(cfg.raw_data_file);
    if ~isfield(loaded_raw, 'X')
        error('svd_subtraction_example_workflow:MissingX', ...
            'Raw MAT file must contain variable ''X''.');
    end
    if ~isempty(cfg.rng_seed)
        rng(cfg.rng_seed, 'twister');
    end
    lib = preprocessFractions(loaded_raw.X, ppm_target, ...
        cfg.noise_region, cfg.remove_regions, false);
end

if isempty(cfg.component_idx)
    if ~isfield(best_state, 'rank_idx') || ~isfield(best_state, 'optimal_k')
        error('svd_subtraction_example_workflow:MissingBasisFields', ...
            'best_state must contain ''rank_idx'' and ''optimal_k'' to subtract the full optimized-k basis set.');
    end
    comp_to_subtract = best_state.rank_idx(1:best_state.optimal_k);
else
    comp_to_subtract = cfg.component_idx;
end

if cfg.save_outputs && (cfg.show_residual_3d || cfg.show_overlay_2d)
    fig_dir = fullfile(cfg.results_dir, 'figures', 'subtraction');
    if ~isfolder(fig_dir), mkdir(fig_dir); end
end

fprintf('Subtracting %d components (full optimal_k basis set) from library (%d fractions x %d ppm points) using ''%s''.\n', ...
    numel(comp_to_subtract), size(lib, 1), numel(ppm_target), cfg.method);

%% 4. SUBTRACT COMPONENTS & COMPUTE VARIANCE METRICS
fprintf('\n--- Running Spectral Library Subtraction ---\n');

[lib_residual, lib_delete, subtraction_metrics, fig_residual_3d] = subtractComponentsSVD( ...
    lib, ppm_target, best_state, comp_to_subtract, ...
    cfg.method, cfg.ppm_window, cfg.frac_idx, ...
    cfg.show_residual_3d, cfg.verbose);

if cfg.show_residual_3d && cfg.save_outputs && ~isempty(fig_residual_3d) && isgraphics(fig_residual_3d)
    print(fig_residual_3d, fullfile(fig_dir, 'residual_library_3d.png'), '-dpng', '-r300');
end

%% 5. OPTIONAL 2D OVERLAY PLOT (ORIGINAL VS. RESIDUAL)
if cfg.show_overlay_2d
    fig_overlay_2d = figure('Name', 'Original vs. Residual Library Overlay', ...
                            'Color', 'w', 'Position', [10, 10, 1000, 1200]);
    plot(ppm_target, lib, 'k', ppm_target, lib_residual, 'r');
    set(gca, 'XDir', 'reverse', 'TickDir', 'out');
    ylim(cfg.overlay_ylim);
    xlabel('Chemical Shift (ppm)');
    ylabel('Signal Intensity');
    title(sprintf('Original (Black) vs. Residual (Red) — %s', cfg.method), ...
          'Interpreter', 'none');

    if cfg.save_outputs && isgraphics(fig_overlay_2d)
        print(fig_overlay_2d, fullfile(fig_dir, 'original_vs_residual_2d.png'), '-dpng', '-r300');
    end
else
    fig_overlay_2d = [];
end

if cfg.save_outputs && (cfg.show_residual_3d || cfg.show_overlay_2d)
    fprintf('Saved subtraction PNG figure(s) to: %s\n', fig_dir);
end

fprintf('\nTotal archived SVD subtraction workflow runtime: %.1f s.\n', toc(total_timer));
