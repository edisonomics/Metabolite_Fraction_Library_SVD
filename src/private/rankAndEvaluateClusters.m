function [purity_score, score_covar, score_sim, score_intensity, rank_idx, optimal_k, max_var_explained, cumulative_var_explained] = rankAndEvaluateClusters(lib, lib_norm_sq, node_max_peaks, spec_correlations, corr_matrix, cluster_idx, consensus_spectra, consensus_elutions, weight_covar, weight_sim, weight_intensity)
% RANKANDEVALUATECLUSTERS Score clusters by purity and evaluate incremental reconstruction.
%   Computes per-cluster covariance correlation, intra-cluster similarity,
%   and normalized peak intensity scores, combines them into a weighted
%   purity confidence score to rank clusters descending, and incrementally
%   reconstructs lib to identify the optimal component count k maximizing
%   cumulative variance explained (%).
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       lib                - MxP spectral library matrix.
%       lib_norm_sq        - Scalar squared Frobenius norm of lib (norm(lib, 'fro')^2).
%       node_max_peaks     - Nx1 vector of maximum peak intensity per node.
%       spec_correlations  - Nx1 vector of target covariance correlations per node.
%       corr_matrix        - NxN pairwise correlation matrix across nodes.
%       cluster_idx        - Nx1 cluster assignment vector (1..C).
%       consensus_spectra  - CxP matrix of cluster consensus spectra.
%       consensus_elutions - CxM matrix of cluster elution profiles.
%       weight_covar       - Purity weight for target covariance correlation.
%       weight_sim         - Purity weight for intra-cluster similarity.
%       weight_intensity   - Purity weight for initial library peak intensity.
%
%   Outputs:
%       purity_score             - Cx1 weighted purity confidence score per cluster.
%       score_covar              - Cx1 median target correlation score per cluster.
%       score_sim                - Cx1 mean intra-cluster correlation score per cluster.
%       score_intensity          - Cx1 normalized median peak intensity per cluster.
%       rank_idx                 - Cx1 cluster indices sorted by purity_score descending.
%       optimal_k                - Component count maximizing cumulative_var_explained.
%       max_var_explained        - Maximum cumulative variance explained (%).
%       cumulative_var_explained - Cx1 vector of cumulative variance explained (%) for k = 1..C.

    num_total_clusters = max(cluster_idx);
    score_covar     = zeros(num_total_clusters, 1);
    score_sim       = zeros(num_total_clusters, 1);
    score_intensity = zeros(num_total_clusters, 1);

    % 1. Purity confidence scoring
    for c = 1:num_total_clusters
        members = find(cluster_idx == c);
        num_members = length(members);

        score_covar(c) = median(spec_correlations(members));
        if num_members > 1
            sub_corr = corr_matrix(members, members);
            pairwise_vals = sub_corr(triu(true(num_members), 1));
            score_sim(c) = mean(pairwise_vals);
        else
            score_sim(c) = 1.0;
        end
        score_intensity(c) = median(node_max_peaks(members));
    end

    score_covar(score_covar < 0) = 0;
    score_sim(score_sim < 0) = 0;

    max_intensity_overall = max(score_intensity);
    if max_intensity_overall > 0
        score_intensity = score_intensity ./ max_intensity_overall;
    end

    total_weight = weight_covar + weight_sim + weight_intensity;
    if total_weight > 0
        purity_score = ((score_covar .* weight_covar) + ...
                        (score_sim .* weight_sim) + ...
                        (score_intensity .* weight_intensity)) ./ total_weight;
    else
        purity_score = zeros(num_total_clusters, 1);
    end

    [~, rank_idx] = sort(purity_score, 'descend');

    % 2. Incremental reconstruction and cumulative variance explained per k
    reordered_spectra  = consensus_spectra(rank_idx, :);
    reordered_elutions = consensus_elutions(rank_idx, :);

    cumulative_var_explained = zeros(num_total_clusters, 1);
    for k = 1:num_total_clusters
        partial_recon = reordered_elutions(1:k, :)' * reordered_spectra(1:k, :);
        if lib_norm_sq > 0
            residual_norm_sq = norm(lib - partial_recon, 'fro')^2;
            cumulative_var_explained(k) = (1 - (residual_norm_sq / lib_norm_sq)) * 100;
        else
            cumulative_var_explained(k) = 100;
        end
    end

    [max_var_explained, optimal_k] = max(cumulative_var_explained);

end

