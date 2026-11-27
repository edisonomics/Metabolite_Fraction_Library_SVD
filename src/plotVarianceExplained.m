function [variance_explained_pct, sort_idx, fig] = plotVarianceExplained(lib, best_state, k_markers, verbose)
% PLOTVARIANCEEXPLAINED Compute and plot cumulative variance explained by intensity-sorted components.
%   Sorts the winning consensus components by descending peak elution
%   intensity, incrementally reconstructs the library up to optimal_k to
%   calculate cumulative variance explained (%), and generates a clean
%   publication-style scree plot with vertical reference markers.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       lib        - MxP spectral library matrix (rows = fractions, columns = ppm).
%       best_state - Struct returned by optimizeClusterCutoff containing
%                    'optimal_k', 'consensus_spectra', and 'consensus_elutions'.
%       k_markers  - (Optional) Row vector of component counts at which to draw
%                    dotted vertical reference lines (default: [5, 10, 25]).
%       verbose    - (Optional) Logical flag to print final variance explained
%                    to the command window (default: true).
%
%   Outputs:
%       variance_explained_pct - max_k x 1 vector of cumulative variance explained (%).
%       sort_idx               - Cx1 permutation sorting components by descending
%                                maximum elution intensity.
%       fig                    - Figure handle for the variance explained plot.

    if nargin < 3 || isempty(k_markers)
        k_markers = [5, 10, 25];
    end
    if nargin < 4 || isempty(verbose)
        verbose = true;
    end

    if isempty(fieldnames(best_state)) || isempty(best_state.consensus_spectra)
        variance_explained_pct = [];
        sort_idx = [];
        fig = [];
        warning('plotVarianceExplained:EmptyState', ...
            'best_state is empty; skipping variance explained calculation.');
        return;
    end

    % 1. Sort components by highest elution intensity (descending)
    opt_k    = best_state.optimal_k;
    spectra  = best_state.consensus_spectra;
    elutions = best_state.consensus_elutions;

    elution_max_ints = max(elutions, [], 2);
    [~, sort_idx] = sort(elution_max_ints, 'descend');

    spectra  = spectra(sort_idx, :);
    elutions = elutions(sort_idx, :);

    max_k = min(size(spectra, 1), opt_k);

    % 2. Incrementally calculate cumulative variance explained (rank-1 updates)
    lib_norm_sq = norm(lib, 'fro')^2;
    variance_explained_pct = zeros(max_k, 1);
    residual = lib;

    for k = 1:max_k
        residual = residual - elutions(k, :)' * spectra(k, :);
        if lib_norm_sq > 0
            residual_norm_sq = norm(residual, 'fro')^2;
            variance_explained_pct(k) = (1 - (residual_norm_sq / lib_norm_sq)) * 100;
        else
            variance_explained_pct(k) = 100;
        end
    end

    % 3. Publication-quality plot
    fig = figure('Name', 'Cumulative Variance Explained (Intensity Sorted)', ...
                 'Position', [150, 150, 600, 600], 'Color', 'w');
    hold on;

    min_y = max(0, min(variance_explained_pct) - 5);
    max_y = min(100, max(variance_explained_pct) + 5);

    % Vertical reference lines in the background
    for m = 1:length(k_markers)
        plot([k_markers(m), k_markers(m)], [min_y, max_y], ':', ...
             'Color', 'k', 'LineWidth', 1.5);
    end

    % Scree curve on top so white markers mask lines behind them
    plot(1:max_k, variance_explained_pct, '-o', 'LineWidth', 2.5, 'MarkerSize', 7, ...
         'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'k');

    set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on', ...
             'TickDir', 'out', 'Layer', 'top');
    xlim([0, max_k + 1]);
    ylim([min_y, max_y]);
    hold off;

    if verbose
        fprintf('Variance explained at optimal cutoff (k = %d): %.2f%%\n', ...
            max_k, variance_explained_pct(max_k));
    end
end
