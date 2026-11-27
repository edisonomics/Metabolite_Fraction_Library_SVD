%% GISSMO Reference Library Spectral Matching Workflow
% =========================================================================
% Matches purified component spectra produced by mFL_svd_example_workflow.m
% against a GISSMO NMR reference spectral library and optionally plots the
% top spectral overlays.
%
% Pipeline Stages:
%   1. Load or validate SVD decomposition outputs (best_state and ppm_work)
%   2. Load and L2-normalize the GISSMO reference spectral library
%   3. Align and match component spectra against GISSMO via shift scanning
%   4. Export sorted GISSMO matches table and optional 300-DPI PNG overlays
%
% Required Input Data (already in workspace from mFL_svd_example_workflow.m,
% loaded automatically from <cfg.results_dir>/data/svd_results.mat, or
% loaded from cfg.svd_results_file):
%   best_state - Struct returned by optimizeClusterCutoff containing
%                'basis_spectra' (optimal_k x P) and 'consensus_spectra'.
%   ppm_work   - 1xP numeric chemical shift vector (ppm). Falls back to
%                'ppm' if 'ppm_work' is not present in the workspace.
%
% Primary Outputs (left in workspace and exported to cfg.results_dir):
%   gissmo_lib            - Struct returned by gissmo_load_library.
%   gissmo_match          - Struct returned by gissmo_match_spectra containing
%                           .results, .sorted_results, .ppm_lib,
%                           .shifted_spectra, .matched_lib_spectra, and .options.
%   gissmo_results        - Table of matches ordered by Target_ID.
%   sorted_gissmo_results - Table of matches ranked by Score (descending).
%   figs_gissmo           - Array of MATLAB Figure handles (when
%                           cfg.show_match_plots is true).
%
% External Functions Required:
%   gissmo_load_library, gissmo_match_spectra, gissmo_plot_matches
%
% Required MATLAB Toolboxes:
%   None (compatible with MATLAB R2019b or later)
% =========================================================================

%% 1. VERIFY PATH SETUP
% Run from the repository root directory (the folder containing setup.m and src/).
if isempty(which('gissmo_match_spectra'))
    if isfile('setup.m')
        setup;
    else
        error('gissmo_matching_example_workflow:PathNotSet', ...
            'Run this workflow from the repository root directory (the folder containing setup.m and src/).');
    end
else
    disp('Repository is on MATLAB path and GISSMO matching workflow is ready.');
end

%% 2. USER CONFIGURATION
% Edit the parameters below to match your dataset and matching preferences.

% --- Output Directory & Exports ---
% Set to any relative or absolute directory path on your computer.
cfg.results_dir  = 'results';
cfg.save_outputs = true;      % Export sorted GISSMO matches CSV and PNG plots

% --- Input Data ---
% Leave cfg.svd_results_file empty ('') to use 'best_state' and 'ppm_work'
% already in the workspace, or automatically fall back to
% <cfg.results_dir>/data/svd_results.mat.
cfg.svd_results_file  = '';
cfg.use_basis_spectra = false; % false = best_state.consensus_spectra (all clusters),
                               % true  = best_state.basis_spectra (optimal_k only)

% --- GISSMO Reference Library ---
cfg.gissmo_lib_file = fullfile('test_data', 'GISSMO_Lib.mat');

% --- Matching Options (gissmo_match_spectra) ---
cfg.min_frac    = 0.0;      % Peak removal threshold as fraction of max peak (0..1)
cfg.shift_max   = 50;       % Max point shift tested in either direction
cfg.shift_step  = 5;        % Point step size between tested shifts
cfg.search_mode = "global"; % "global" (uniform shift) or "split" (L/R halves)
cfg.verbose     = true;     % Print formatted match summary tables to console

% --- Plotting Options (gissmo_plot_matches) ---
cfg.show_match_plots   = false; % Set true to generate overlaid match figures
cfg.plot_target_ids    = [];    % Specific Target_IDs to plot ([] = all)
cfg.plot_min_score     = -Inf;  % Minimum cosine score to plot (e.g., 0.7)
cfg.plot_sort_by_score = false; % Plot matches in Target_ID order (true = score descending)
cfg.plot_save_dir      = "";    % Optional override ("" defaults to <cfg.results_dir>/figures/gissmo)

%% 3. LOAD & VALIDATE INPUTS
total_timer = tic;

load_path = cfg.svd_results_file;
if isempty(load_path) && (~exist('best_state', 'var') || (~exist('ppm_work', 'var') && ~exist('ppm', 'var')))
    default_ckpt = fullfile(cfg.results_dir, 'data', 'svd_results.mat');
    if isfile(default_ckpt)
        load_path = default_ckpt;
    end
end

if ~isempty(load_path)
    if ~isfile(load_path)
        error('gissmo_matching_example_workflow:FileNotFound', ...
            'SVD results file not found: %s', load_path);
    end
    fprintf('Loading SVD results from: %s\n', load_path);
    loaded = load(load_path);
    if ~isfield(loaded, 'best_state')
        error('gissmo_matching_example_workflow:MissingBestState', ...
            'MAT file must contain variable ''best_state''.');
    end
    best_state = loaded.best_state;
    if isfield(loaded, 'ppm_work')
        ppm_target = loaded.ppm_work;
    elseif isfield(loaded, 'ppm')
        ppm_target = loaded.ppm;
    else
        error('gissmo_matching_example_workflow:MissingPPM', ...
            'MAT file must contain variable ''ppm_work'' or ''ppm''.');
    end
else
    if ~exist('best_state', 'var')
        error('gissmo_matching_example_workflow:MissingBestState', ...
            ['Variable ''best_state'' not found in workspace and no checkpoint was found at %s. ' ...
             'Run workflows/mFL_svd_example_workflow.m first or set cfg.svd_results_file in Section 2.'], ...
             fullfile(cfg.results_dir, 'data', 'svd_results.mat'));
    end
    if exist('ppm_work', 'var')
        ppm_target = ppm_work;
    elseif exist('ppm', 'var')
        ppm_target = ppm;
    else
        error('gissmo_matching_example_workflow:MissingPPM', ...
            ['Neither ''ppm_work'' nor ''ppm'' was found in workspace. Run ' ...
             'workflows/mFL_svd_example_workflow.m first or set cfg.svd_results_file in Section 2.']);
    end
end

if cfg.use_basis_spectra
    if ~isfield(best_state, 'basis_spectra') || isempty(best_state.basis_spectra)
        error('gissmo_matching_example_workflow:EmptySpectra', ...
            'best_state.basis_spectra is missing or empty.');
    end
    spec_matrix = best_state.basis_spectra;
    spec_source_label = 'best_state.basis_spectra';
else
    if ~isfield(best_state, 'consensus_spectra') || isempty(best_state.consensus_spectra)
        error('gissmo_matching_example_workflow:EmptySpectra', ...
            'best_state.consensus_spectra is missing or empty.');
    end
    spec_matrix = best_state.consensus_spectra;
    spec_source_label = 'best_state.consensus_spectra';
end

if ~isfile(cfg.gissmo_lib_file)
    error('gissmo_matching_example_workflow:MissingGISSMO', ...
        'GISSMO reference library not found: %s', cfg.gissmo_lib_file);
end

fprintf('Prepared %d target spectra from %s (%d ppm points).\n', ...
    size(spec_matrix, 1), spec_source_label, numel(ppm_target));

%% 4. LOAD GISSMO REFERENCE LIBRARY
fprintf('\n--- Loading GISSMO Reference Library ---\n');
gissmo_lib = gissmo_load_library(cfg.gissmo_lib_file);
fprintf('Loaded %d GISSMO reference spectra (%d ppm points).\n', ...
    size(gissmo_lib.spectra_norm, 1), numel(gissmo_lib.ppm));

%% 5. MATCH SPECTRA AGAINST GISSMO LIBRARY
fprintf('\n--- Matching Spectra Against GISSMO Library ---\n');
t_match = tic;

gissmo_match = gissmo_match_spectra( ...
    spec_matrix, ppm_target, gissmo_lib, ...
    'MinFrac',     cfg.min_frac, ...
    'ShiftMax',    cfg.shift_max, ...
    'ShiftStep',   cfg.shift_step, ...
    'SearchMode',  cfg.search_mode, ...
    'PlotMatches', false, ...
    'Verbose',     cfg.verbose);

gissmo_results        = gissmo_match.results;
sorted_gissmo_results = gissmo_match.sorted_results;

fprintf('GISSMO matching finished in %.1f s.\n', toc(t_match));

if cfg.save_outputs
    table_dir = fullfile(cfg.results_dir, 'tables');
    if ~isfolder(table_dir), mkdir(table_dir); end
    writetable(sorted_gissmo_results, fullfile(table_dir, 'gissmo_matches_sorted.csv'));
    fprintf('Saved sorted GISSMO matches table to: %s\n', ...
        fullfile(table_dir, 'gissmo_matches_sorted.csv'));
end

%% 6. PLOT MATCHED SPECTRA
if cfg.show_match_plots
    fprintf('\n--- Plotting GISSMO Matches ---\n');
    if strlength(cfg.plot_save_dir) > 0
        effective_save_dir = cfg.plot_save_dir;
    elseif cfg.save_outputs
        effective_save_dir = string(fullfile(cfg.results_dir, 'figures', 'gissmo'));
    else
        effective_save_dir = "";
    end

    figs_gissmo = gissmo_plot_matches( ...
        gissmo_match, ...
        'TargetIDs',   cfg.plot_target_ids, ...
        'MinScore',    cfg.plot_min_score, ...
        'SortByScore', cfg.plot_sort_by_score, ...
        'SaveDir',     effective_save_dir);
end

fprintf('\nTotal GISSMO matching workflow runtime: %.1f s.\n', toc(total_timer));
