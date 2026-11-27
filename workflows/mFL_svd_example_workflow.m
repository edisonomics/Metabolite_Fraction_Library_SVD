%% NMR / Spectral Fraction Library Decomposition Workflow
% =========================================================================
% End-to-end pipeline for extracting purified component spectra and elution
% profiles from a 2D NMR/spectral fraction library.
%
% Pipeline Stages:
%   1. Downsample / interpolate spectra onto a uniform ppm grid (optional)
%   2. Replace solvent/reference regions with synthetic noise and clip negatives
%   3. Detect peaks per fraction, trace across consecutive fractions, and
%      remove redundant traces
%   4. Extract component spectra & elutions via per-node SVD-ALS
%   5. Quality-control filter against target covariance spectra
%   6. Hierarchically cluster purified spectra by Pearson correlation
%   7. Sweep correlation cutoffs to select the optimal spectral basis set
%   8. Visualize cumulative variance explained and stepped 3D reconstructions
%   9. Export PNG figures, variance explained table, consensus CSVs, and
%      the SVD checkpoint for downstream workflows
%
% Required Input Data (loaded automatically from cfg.data_file or already
% present in workspace when cfg.data_file is empty):
%   X   - MxP numeric matrix of spectral intensities (rows = fractions,
%         columns = chemical shift points).
%   ppm - 1xP numeric vector of chemical shift values (ppm).
%
% Primary Outputs (left in workspace and exported to cfg.results_dir):
%   best_state             - Struct containing the optimal basis set
%                            (best_state.consensus_spectra,
%                             best_state.basis_spectra,
%                             best_state.basis_elutions,
%                             best_state.optimal_reconstructed_lib, etc.).
%   sweep_results          - Struct of error and component count metrics
%                            across the correlation cutoff sweep.
%   variance_explained_pct - Cumulative variance explained (%) up to
%                            best_state.optimal_k (intensity-sorted).
%
% Required MATLAB Toolboxes:
%   - Signal Processing Toolbox (findpeaks)
%   - Statistics and Machine Learning Toolbox (zscore, corr, pdist,
%     linkage, dendrogram, cluster)
%   - Optimization Toolbox (optimproblem, fmincon)
%   - Image Processing Toolbox (imdilate)
%   - Parallel Computing Toolbox (optional, if use_parallel = true)
% =========================================================================

%% 1. VERIFY PATH SETUP
% Run from the repository root directory (the folder containing setup.m and src/).
if isempty(which('downsampleData'))
    if isfile('setup.m')
        setup;
    else
        error('mFL_svd_example_workflow:PathNotSet', ...
            'Run this workflow from the repository root directory (the folder containing setup.m and src/).');
    end
else
    disp('Repository is on MATLAB path and workflow is ready.');
end

%% 2. USER CONFIGURATION
% Edit the parameters below to match your dataset and analysis preferences.

% --- Input Data ---
% Default loads test_data/crassa.mat automatically. Provide a different .mat
% file path containing 'X' and 'ppm', OR set to '' to use 'X' and 'ppm'
% already loaded in your MATLAB workspace.
cfg.data_file = fullfile('test_data', 'crassa.mat');

% --- Output Directory & Exports ---
% Set to any relative or absolute directory path on your computer.
cfg.results_dir  = 'results';
cfg.save_outputs = true;      % Export 300-DPI PNG figures, variance table, and data files

% --- Step 1: Downsampling / Interpolation ---
cfg.do_downsample = false;       % Set false to keep original spectral resolution
cfg.num_points    = 8192;        % Target number of ppm points if do_downsample = true

% --- Step 2: Preprocessing (Noise Replacement & Negative Clipping) ---
% Define [min_ppm, max_ppm] for a signal-free baseline noise region and a
% Kx2 matrix of regions to remove (e.g., DSS, water, formate).
% Adjust these chemical shift bounds to match your spectra!
cfg.noise_region        = [0.2, 0.6];
cfg.remove_regions      = [ ...
    -0.244,   0.331;    % Example: DSS / reference peak
     4.73,    4.82;     % Example: residual water
     8.40644, 8.48747   % Example: formate
];
cfg.rng_seed             = 42;
cfg.show_preprocess_plot = true; % Calls displaySpectra3D on preprocessed library

% --- Step 3: Peak Detection, Tracing & Redundancy Removal ---
cfg.noise_threshold      = 0.05;  % Minimum peak height for findpeaks
cfg.ppm_tol              = 0.003; % Max ppm shift allowed when linking across fractions
cfg.min_trace_length     = 3;     % Minimum trace length (including 0-padded anchors)
cfg.trace_corr_threshold = 0.99;  % Pearson r cutoff to drop redundant same-span traces
cfg.show_peaks_plot      = true;  % 3D overlay of raw detected peaks
cfg.show_traces_plot     = true;  % 3D overlay of deduplicated traces

% --- Step 4: SVD-ALS Component Extraction ---
cfg.svd_opts.use_parallel = false; % Set true to run parfor across nodes
cfg.svd_opts.verbose      = true;  % Print start/end extraction summary

% --- Step 5: Quality Control ---
cfg.qc_corr_threshold = 0.2;   % Minimum Pearson r vs. target covariance spectrum
cfg.qc_target_line    = 0.8;   % Reference dashed line drawn on QC scatter plot
cfg.show_qc_plot      = false; % Display QC scatter + histogram figure

% --- Step 6: Hierarchical Clustering ---
cfg.linkage_method    = 'complete'; % Linkage criterion for hierarchical tree
cfg.show_cluster_plot = true;       % Display dendrogram + sorted correlation heatmap

% --- Step 7: Correlation Cutoff Optimization (Deep Sweep) ---
cfg.sweep_opts.cutoff_range           = 0.5:0.02:0.95;
cfg.sweep_opts.weight_covar           = 0.8;   % Purity weight: target covariance match
cfg.sweep_opts.weight_sim             = 0.0;   % Purity weight: intra-cluster similarity
cfg.sweep_opts.weight_intensity       = 0.2;   % Purity weight: initial peak intensity
cfg.sweep_opts.max_iters              = 20;    % Max iterations for elution selection
cfg.sweep_opts.snap_spectra_to_winner = false; % true = replace median spectrum with winner
cfg.sweep_opts.use_parallel           = false; % Set true to run parfor across cutoffs
cfg.sweep_opts.verbose                = true;  % Print per-cutoff progress and summary

% --- Step 8: Final Visualizations ---
cfg.show_variance_plot = true;         % Scree plot of cumulative variance explained
cfg.variance_k_markers = [5, 10, 25];  % Vertical dotted reference lines on scree plot
cfg.show_stepped_recon = true;         % 3D reconstructions at stepped component counts
cfg.recon_k_steps      = [5, 10, 25];  % Component counts for stepped 3D reconstructions
cfg.recon_include_full = true;         % Also display full optimal_k reconstruction
cfg.recon_include_raw  = true;         % Also display raw library for comparison


%% 3. LOAD & VALIDATE INPUT DATA & OUTPUT FOLDERS
total_timer = tic;

if ~isempty(cfg.data_file)
    if ~isfile(cfg.data_file)
        error('mFL_svd_example_workflow:FileNotFound', ...
            'Input data file not found: %s', cfg.data_file);
    end
    fprintf('Loading input data from: %s\n', cfg.data_file);
    loaded = load(cfg.data_file);
    if ~isfield(loaded, 'X') || ~isfield(loaded, 'ppm')
        error('mFL_svd_example_workflow:MissingFields', ...
            'MAT file must contain variables ''X'' and ''ppm''.');
    end
    X   = loaded.X;
    ppm = loaded.ppm;
elseif ~exist('X', 'var') || ~exist('ppm', 'var')
    error('mFL_svd_example_workflow:MissingInput', ...
        ['Variables ''X'' and ''ppm'' not found in workspace. ' ...
         'Either load them first or set cfg.data_file in Section 2 (USER CONFIGURATION).']);
end

if ~ismatrix(X) || ~isvector(ppm) || size(X, 2) ~= numel(ppm)
    error('mFL_svd_example_workflow:DimensionMismatch', ...
        'Input dimensions mismatch: X is %dx%d while ppm has %d elements.', ...
        size(X, 1), size(X, 2), numel(ppm));
end

if cfg.save_outputs
    fig_dir   = fullfile(cfg.results_dir, 'figures', 'svd');
    table_dir = fullfile(cfg.results_dir, 'tables');
    data_dir  = fullfile(cfg.results_dir, 'data');
    if ~isfolder(fig_dir),   mkdir(fig_dir);   end
    if ~isfolder(table_dir), mkdir(table_dir); end
    if ~isfolder(data_dir),  mkdir(data_dir);  end
end

fprintf('Loaded spectral library: %d fractions x %d ppm points.\n\n', size(X, 1), size(X, 2));

%% 4. DOWNSAMPLE & PREPROCESS SPECTRAL LIBRARY
fprintf('--- Stage 1-2: Downsampling & Preprocessing ---\n');

if cfg.do_downsample
    [X_work, ppm_work] = downsampleData(X, ppm, cfg.num_points);
    fprintf('Interpolated spectra to %d ppm points.\n', length(ppm_work));
else
    X_work   = X;
    ppm_work = ppm;
end

if ~isempty(cfg.rng_seed)
    rng(cfg.rng_seed, 'twister');
end

lib = preprocessFractions(X_work, ppm_work, cfg.noise_region, ...
                          cfg.remove_regions, cfg.show_preprocess_plot);

if cfg.show_preprocess_plot && cfg.save_outputs
    fig_preprocess = gcf;
    print(fig_preprocess, fullfile(fig_dir, 'preprocessed_library_3d.png'), '-dpng', '-r300');
end

%% 5. PEAK DETECTION, TRACING & REDUNDANCY REMOVAL
fprintf('\n--- Stage 3: Peak Detection, Tracing & Deduplication ---\n');

peaks_struct = findLibraryPeaks(lib, ppm_work, cfg.noise_threshold);
if cfg.show_peaks_plot
    fig_peaks = plotPeaksAndTraces(lib, ppm_work, peaks_struct);
    if cfg.save_outputs && ~isempty(fig_peaks) && isgraphics(fig_peaks)
        print(fig_peaks, fullfile(fig_dir, 'detected_peaks_3d.png'), '-dpng', '-r300');
    end
end

nodes = tracePeaks(peaks_struct, cfg.ppm_tol, cfg.min_trace_length);
nodes_filtered = removeRedundantTraces(nodes, cfg.trace_corr_threshold);

fprintf('Detected traces: %d raw -> %d after redundancy removal.\n', ...
    length(nodes), length(nodes_filtered));

if isempty(nodes_filtered)
    error('mFL_svd_example_workflow:NoTraces', ...
        'No peak traces survived filtering. Consider lowering cfg.noise_threshold or cfg.min_trace_length.');
end

if cfg.show_traces_plot
    fig_traces = plotPeaksAndTraces(lib, ppm_work, nodes_filtered);
    if cfg.save_outputs && ~isempty(fig_traces) && isgraphics(fig_traces)
        print(fig_traces, fullfile(fig_dir, 'deduplicated_traces_3d.png'), '-dpng', '-r300');
    end
end

%% 6. OPTIONAL: REORDER TRACES BY INTENSITY
index_max = zeros(length(nodes_filtered), 1);
for i = 1:length(index_max)
    index_max(i) = max(nodes_filtered(i).peaks);
end
[~, idx] = sort(index_max, 'descend');
nodes_filtered = nodes_filtered(idx);
clear i index_max idx

%% 7. SVD-ALS COMPONENT EXTRACTION
fprintf('\n--- Stage 4: SVD-ALS Component Extraction ---\n');
t_svd = tic;

nodes_svd = extractComponentsSVD(lib, nodes_filtered, cfg.svd_opts);

fprintf('SVD-ALS extraction finished in %.1f s.\n', toc(t_svd));

%% 8. QUALITY CONTROL FILTERING
fprintf('\n--- Stage 5: Quality Control Filtering ---\n');

[nodes_qc, spec_correlations_filt, spec_correlations_all, fig_qc] = ...
    filterComponentsQC(lib, nodes_svd, cfg.qc_corr_threshold, ...
                       cfg.show_qc_plot, cfg.qc_target_line);

if cfg.show_qc_plot && cfg.save_outputs && ~isempty(fig_qc) && isgraphics(fig_qc)
    print(fig_qc, fullfile(fig_dir, 'qc_correlation_filter.png'), '-dpng', '-r300');
end

fprintf('QC retained %d / %d components (r >= %.2f).\n', ...
    length(nodes_qc), length(nodes_svd), cfg.qc_corr_threshold);

if isempty(nodes_qc)
    error('mFL_svd_example_workflow:NoNodesAfterQC', ...
        'No components survived QC filtering. Consider lowering cfg.qc_corr_threshold.');
end

%% 9. HIERARCHICAL CLUSTERING OF EXTRACTED SPECTRA
fprintf('\n--- Stage 6: Hierarchical Clustering ---\n');

[cluster_info, fig_cluster] = clusterExtractedSpectra( ...
    nodes_qc, cfg.linkage_method, cfg.show_cluster_plot);

if cfg.show_cluster_plot && cfg.save_outputs && ~isempty(fig_cluster) && isgraphics(fig_cluster)
    print(fig_cluster, fullfile(fig_dir, 'hierarchical_clustering.png'), '-dpng', '-r300');
end

%% 10. CORRELATION CUTOFF OPTIMIZATION (DEEP SWEEP)
fprintf('\n--- Stage 7: Correlation Cutoff Optimization ---\n');
t_sweep = tic;

[best_state, sweep_results] = optimizeClusterCutoff( ...
    lib, nodes_qc, spec_correlations_filt, cluster_info, cfg.sweep_opts);

fprintf('Cutoff sweep finished in %.1f s.\n', toc(t_sweep));

%% 11. FINAL VISUALIZATIONS (VARIANCE EXPLAINED & STEPPED RECONSTRUCTIONS)
fprintf('\n--- Stage 8: Final Visualizations ---\n');

if cfg.show_variance_plot
    [variance_explained_pct, intensity_sort_idx, fig_variance] = ...
        plotVarianceExplained(lib, best_state, cfg.variance_k_markers, true);
    if cfg.save_outputs && ~isempty(fig_variance) && isgraphics(fig_variance)
        print(fig_variance, fullfile(fig_dir, 'variance_explained.png'), '-dpng', '-r300');
    end
end

if cfg.show_stepped_recon
    figs_stepped = plotSteppedReconstructions( ...
        lib, ppm_work, best_state, cfg.recon_k_steps, ...
        cfg.recon_include_full, cfg.recon_include_raw);
    if cfg.save_outputs
        for i_step = 1:numel(figs_stepped.steps)
            f_step = figs_stepped.steps(i_step);
            if isgraphics(f_step)
                print(f_step, fullfile(fig_dir, sprintf('reconstruction_k%02d_3d.png', figs_stepped.step_ks(i_step))), '-dpng', '-r300');
            end
        end
        if ~isempty(figs_stepped.full_reconstruct) && isgraphics(figs_stepped.full_reconstruct)
            print(figs_stepped.full_reconstruct, fullfile(fig_dir, 'reconstruction_full_3d.png'), '-dpng', '-r300');
        end
        if ~isempty(figs_stepped.raw) && isgraphics(figs_stepped.raw)
            print(figs_stepped.raw, fullfile(fig_dir, 'reconstruction_raw_3d.png'), '-dpng', '-r300');
        end
    end
end

%% 12. EXPORT TABLE, CONSENSUS CSVS & CHECKPOINT
if cfg.save_outputs
    fprintf('\n--- Stage 9: Exporting Results ---\n');

    % 1. Variance explained table
    if exist('variance_explained_pct', 'var') && ~isempty(variance_explained_pct)
        n_var = numel(variance_explained_pct);
        variance_table = table( ...
            (1:n_var)', ...
            intensity_sort_idx(1:n_var), ...
            variance_explained_pct(:), ...
            'VariableNames', {'Rank', 'Component_ID', 'Cumulative_Variance_Explained_Pct'});
    else
        n_var = best_state.optimal_k;
        variance_table = table( ...
            (1:n_var)', ...
            best_state.rank_idx(1:n_var), ...
            best_state.cumulative_var_explained(1:n_var), ...
            'VariableNames', {'Rank', 'Component_ID', 'Cumulative_Variance_Explained_Pct'});
    end
    writetable(variance_table, fullfile(table_dir, 'variance_explained.csv'));

    % 2. Consensus spectra & elutions as separate CSVs
    num_comp = size(best_state.consensus_spectra, 1);
    comp_names = arrayfun(@(k) sprintf('Component_%d', k), 1:num_comp, 'UniformOutput', false);

    spec_tbl = array2table([ppm_work(:), best_state.consensus_spectra'], ...
        'VariableNames', [{'ppm'}, comp_names]);
    writetable(spec_tbl, fullfile(data_dir, 'consensus_spectra.csv'));

    num_fracs = size(best_state.consensus_elutions, 2);
    elut_tbl = array2table([(1:num_fracs)', best_state.consensus_elutions'], ...
        'VariableNames', [{'Fraction'}, comp_names]);
    writetable(elut_tbl, fullfile(data_dir, 'consensus_elutions.csv'));

    % 3. Checkpoint MAT file for downstream workflows
    checkpoint_file = fullfile(data_dir, 'svd_results.mat');
    if exist('variance_explained_pct', 'var') && exist('intensity_sort_idx', 'var')
        save(checkpoint_file, 'best_state', 'sweep_results', 'lib', ...
            'X_work', 'ppm_work', 'variance_explained_pct', 'intensity_sort_idx');
    else
        save(checkpoint_file, 'best_state', 'sweep_results', 'lib', ...
            'X_work', 'ppm_work');
    end

    fprintf('Saved PNG figures, variance table, consensus CSVs, and checkpoint to: %s\n', cfg.results_dir);
end

fprintf('\nTotal workflow runtime: %.1f s.\n', toc(total_timer));
