function [best_state, sweep_results] = optimizeClusterCutoff(lib, nodes, spec_correlations, cluster_info, opts)
% OPTIMIZECLUSTERCUTOFF Sweep correlation cutoffs to select the optimal spectral basis set.
%   Evaluates a grid of hierarchical clustering correlation cutoffs. For
%   each cutoff, optimizes cluster elution assignments via alternating
%   member selection, ranks clusters by a weighted purity confidence score,
%   and finds the number of components maximizing cumulative variance
%   explained (%). Returns the winning basis state and full sweep diagnostics.
%
%   External Functions Required:
%       private/optimizeClusterElutions
%       private/rankAndEvaluateClusters
%
%   Toolboxes Required:
%       Statistics and Machine Learning Toolbox (cluster)
%       Parallel Computing Toolbox (optional, when opts.use_parallel is true)
%
%   Inputs:
%       lib               - MxP spectral library matrix (rows = fractions, columns = ppm).
%       nodes             - 1xN struct array (e.g., from filterComponentsQC) with
%                           fields 'fracs', 'peaks', 'optimized_elution', and
%                           'singular_value'.
%       spec_correlations - Nx1 vector of target spectrum correlations matching nodes.
%       cluster_info      - Struct returned by clusterExtractedSpectra containing
%                           'spectra_matrix', 'corr_matrix', and 'clustering_tree'.
%       opts              - (Optional) Struct of sweep and scoring controls:
%                           'cutoff_range'           (default: 0.5:0.02:0.95)
%                           'weight_covar'           (default: 0.8)
%                           'weight_sim'             (default: 0.0)
%                           'weight_intensity'       (default: 0.2)
%                           'max_iters'              (default: 20)
%                           'snap_spectra_to_winner' (default: false)
%                           'use_parallel'           (default: false)
%                           'verbose'                (default: true)
%
%   Outputs:
%       best_state    - Struct for the maximum-variance-explained cutoff containing
%                       'cutoff', 'cluster_idx', 'num_total_clusters',
%                       'consensus_spectra', 'consensus_elutions',
%                       'best_member_selected', 'cluster_sizes', 'rank_idx',
%                       'optimal_k', 'purity_score', 'score_covar', 'score_sim',
%                       'score_intensity', 'optimal_reconstructed_lib',
%                       'cumulative_var_explained', 'basis_spectra'
%                       (optimal_k x P ranked basis), and 'basis_elutions'
%                       (optimal_k x M ranked elutions).
%       sweep_results - Struct summarizing the full sweep with fields
%                       'cutoff_range', 'max_var_explained', 'optimal_ks',
%                       'num_clusters', 'global_best_var_explained', and
%                       'best_cutoff'.

    % Parse optional opts struct
    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    cutoff_range           = getOpt(opts, 'cutoff_range', 0.5:0.02:0.95);
    weight_covar           = getOpt(opts, 'weight_covar', 0.8);
    weight_sim             = getOpt(opts, 'weight_sim', 0.0);
    weight_intensity       = getOpt(opts, 'weight_intensity', 0.2);
    max_iters              = getOpt(opts, 'max_iters', 20);
    snap_spectra_to_winner = getOpt(opts, 'snap_spectra_to_winner', false);
    use_parallel           = getOpt(opts, 'use_parallel', false);
    verbose                = getOpt(opts, 'verbose', true);

    % Input validation
    if isempty(lib) || ~ismatrix(lib) || isscalar(lib)
        error('optimizeClusterCutoff:InvalidLibrary', ...
            'Input ''lib'' must be a non-empty 2D numeric matrix.');
    end

    num_nodes = length(nodes);
    spec_correlations = spec_correlations(:);
    if length(spec_correlations) ~= num_nodes
        error('optimizeClusterCutoff:LengthMismatch', ...
            ['Length of spec_correlations (%d) does not match length of nodes (%d). ' ...
             'Ensure filterComponentsQC uses consistent >= / < boundary logic.'], ...
            length(spec_correlations), num_nodes);
    end

    % Edge case: no nodes surviving QC
    if num_nodes == 0
        best_state = struct();
        sweep_results = struct( ...
            'cutoff_range', cutoff_range(:), ...
            'max_var_explained', [], ...
            'optimal_ks', [], ...
            'num_clusters', [], ...
            'global_best_var_explained', NaN, ...
            'best_cutoff', NaN);
        return;
    end

    if verbose
        fprintf('Starting Deep Sweep for Optimal Correlation Cutoff...\n');
    end

    % Precompute sweep-invariant quantities once
    num_elution_points = size(lib, 1);
    lib_norm_sq = norm(lib, 'fro')^2;
    spectra_matrix = cluster_info.spectra_matrix;
    corr_matrix = cluster_info.corr_matrix;
    clustering_tree = cluster_info.clustering_tree;

    all_elutions = zeros(num_nodes, num_elution_points);
    node_max_peaks = zeros(num_nodes, 1);
    for n = 1:num_nodes
        all_elutions(n, nodes(n).fracs) = ...
            nodes(n).optimized_elution(:)' * nodes(n).singular_value;
        node_max_peaks(n) = max(nodes(n).peaks);
    end
    elution_norms_sq = sum(all_elutions.^2, 2);

    % Setup loop trackers
    num_sweeps              = length(cutoff_range);
    sweep_max_var_explained = zeros(num_sweeps, 1);
    sweep_optimal_ks        = zeros(num_sweeps, 1);
    sweep_num_clusters      = zeros(num_sweeps, 1);

    if use_parallel
        sweep_states = cell(num_sweeps, 1);

        parfor s = 1:num_sweeps
            current_cutoff = cutoff_range(s);

            if num_nodes == 1
                cluster_idx = 1;
            else
                cluster_idx = cluster(clustering_tree, ...
                    'cutoff', 1 - current_cutoff, 'criterion', 'distance');
            end
            num_total_clusters = max(cluster_idx);
            sweep_num_clusters(s) = num_total_clusters;

            [consensus_spectra, consensus_elutions, best_member_selected, cluster_sizes] = ...
                optimizeClusterElutions(lib, spectra_matrix, all_elutions, elution_norms_sq, ...
                    cluster_idx, max_iters, snap_spectra_to_winner);

            [purity_score, score_covar, score_sim, score_intensity, rank_idx, ...
             optimal_k, max_var_explained, cumulative_var_explained] = ...
                rankAndEvaluateClusters(lib, lib_norm_sq, node_max_peaks, spec_correlations, ...
                    corr_matrix, cluster_idx, consensus_spectra, consensus_elutions, ...
                    weight_covar, weight_sim, weight_intensity);

            sweep_max_var_explained(s) = max_var_explained;
            sweep_optimal_ks(s) = optimal_k;

            st = struct();
            st.cutoff = current_cutoff;
            st.cluster_idx = cluster_idx;
            st.num_total_clusters = num_total_clusters;
            st.consensus_spectra = consensus_spectra;
            st.consensus_elutions = consensus_elutions;
            st.best_member_selected = best_member_selected;
            st.cluster_sizes = cluster_sizes;
            st.rank_idx = rank_idx;
            st.optimal_k = optimal_k;
            st.purity_score = purity_score;
            st.score_covar = score_covar;
            st.score_sim = score_sim;
            st.score_intensity = score_intensity;
            st.cumulative_var_explained = cumulative_var_explained;
            sweep_states{s} = st;
        end

        if verbose
            for s = 1:num_sweeps
                fprintf('Cutoff %.2f | Clusters: %3d | Opt Components: %3d | Var Explained: %.2f%%\n', ...
                    cutoff_range(s), sweep_num_clusters(s), sweep_optimal_ks(s), sweep_max_var_explained(s));
            end
        end

        [global_best_var_explained, best_s] = max(sweep_max_var_explained);
        best_state = sweep_states{best_s};
    else
        global_best_var_explained = -inf;
        best_state = struct();

        for s = 1:num_sweeps
            current_cutoff = cutoff_range(s);

            if num_nodes == 1
                cluster_idx = 1;
            else
                cluster_idx = cluster(clustering_tree, ...
                    'cutoff', 1 - current_cutoff, 'criterion', 'distance');
            end
            num_total_clusters = max(cluster_idx);
            sweep_num_clusters(s) = num_total_clusters;

            [consensus_spectra, consensus_elutions, best_member_selected, cluster_sizes] = ...
                optimizeClusterElutions(lib, spectra_matrix, all_elutions, elution_norms_sq, ...
                    cluster_idx, max_iters, snap_spectra_to_winner);

            [purity_score, score_covar, score_sim, score_intensity, rank_idx, ...
             optimal_k, max_var_explained, cumulative_var_explained] = ...
                rankAndEvaluateClusters(lib, lib_norm_sq, node_max_peaks, spec_correlations, ...
                    corr_matrix, cluster_idx, consensus_spectra, consensus_elutions, ...
                    weight_covar, weight_sim, weight_intensity);

            sweep_max_var_explained(s) = max_var_explained;
            sweep_optimal_ks(s) = optimal_k;

            if verbose
                fprintf('Cutoff %.2f | Clusters: %3d | Opt Components: %3d | Var Explained: %.2f%%\n', ...
                    current_cutoff, num_total_clusters, optimal_k, max_var_explained);
            end

            if max_var_explained > global_best_var_explained
                global_best_var_explained = max_var_explained;
                best_state.cutoff = current_cutoff;
                best_state.cluster_idx = cluster_idx;
                best_state.num_total_clusters = num_total_clusters;
                best_state.consensus_spectra = consensus_spectra;
                best_state.consensus_elutions = consensus_elutions;
                best_state.best_member_selected = best_member_selected;
                best_state.cluster_sizes = cluster_sizes;
                best_state.rank_idx = rank_idx;
                best_state.optimal_k = optimal_k;
                best_state.purity_score = purity_score;
                best_state.score_covar = score_covar;
                best_state.score_sim = score_sim;
                best_state.score_intensity = score_intensity;
                best_state.cumulative_var_explained = cumulative_var_explained;
            end
        end
    end

    % Compute final ranked basis set and reconstruction once for winning state
    k_opt = best_state.optimal_k;
    top_idx = best_state.rank_idx(1:k_opt);
    best_state.basis_spectra = best_state.consensus_spectra(top_idx, :);
    best_state.basis_elutions = best_state.consensus_elutions(top_idx, :);
    best_state.optimal_reconstructed_lib = best_state.basis_elutions' * best_state.basis_spectra;

    sweep_results = struct( ...
        'cutoff_range', cutoff_range(:), ...
        'max_var_explained', sweep_max_var_explained, ...
        'optimal_ks', sweep_optimal_ks, ...
        'num_clusters', sweep_num_clusters, ...
        'global_best_var_explained', global_best_var_explained, ...
        'best_cutoff', best_state.cutoff);

    if verbose
        fprintf('\n>>> DEEP SWEEP COMPLETE <<<\n');
        fprintf('Winning Cutoff: r >= %.2f (Var Explained: %.2f%% using %d components)\n\n', ...
            best_state.cutoff, global_best_var_explained, best_state.optimal_k);
    end
end

function val = getOpt(opts, field, default_val)
% GETOPT Safely retrieve an option field or fall back to its default.
    if isfield(opts, field) && ~isempty(opts.(field))
        val = opts.(field);
    else
        val = default_val;
    end
end

