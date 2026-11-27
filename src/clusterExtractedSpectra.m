function [cluster_info, fig] = clusterExtractedSpectra(nodes, linkage_method, show_plot)
% CLUSTEREXTRACTEDSPECTRA Hierarchically cluster extracted spectra by correlation.
%   Extracts purified spectra from a node struct array into a 2D matrix,
%   computes a NaN-sanitized pairwise Pearson correlation matrix (replacing
%   NaNs with 0) and pairwise correlation distances via pdist, builds a
%   hierarchical linkage tree, and optionally renders a paired dendrogram
%   and sorted 2D correlation heatmap.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       Statistics and Machine Learning Toolbox (corr, pdist, linkage, dendrogram)
%
%   Inputs:
%       nodes          - 1xN struct array (e.g., from filterComponentsQC)
%                        containing field 'clean_spectrum'.
%       linkage_method - (Optional) Linkage method passed to linkage()
%                        (default: 'complete').
%       show_plot      - (Optional) Logical flag to display the dendrogram
%                        and sorted correlation heatmap (default: true).
%
%   Outputs:
%       cluster_info - Struct containing clustering results with fields:
%                      'spectra_matrix'  (NxP matrix of clean spectra),
%                      'corr_matrix'     (NxN symmetric Pearson correlation matrix),
%                      'clustering_tree' ((N-1)x3 linkage tree, or [] if N < 2),
%                      'cluster_order'   (1xN dendrogram leaf permutation),
%                      'linkage_method'  (char/string method used).
%       fig          - Figure handle if show_plot is true and N >= 2; otherwise [].

    % Handle optional defaults
    if nargin < 2 || isempty(linkage_method)
        linkage_method = 'complete';
    end
    if nargin < 3 || isempty(show_plot)
        show_plot = true;
    end

    fig = [];
    num_nodes = length(nodes);

    % Edge case: empty input
    if num_nodes == 0
        cluster_info = struct( ...
            'spectra_matrix', [], ...
            'corr_matrix', [], ...
            'clustering_tree', [], ...
            'cluster_order', [], ...
            'linkage_method', linkage_method);
        return;
    end

    % 1. Extract all purified spectra from the struct into a single 2D matrix
    num_spectral_points = length(nodes(1).clean_spectrum);
    spectra_matrix = zeros(num_nodes, num_spectral_points);
    for i = 1:num_nodes
        spectra_matrix(i, :) = nodes(i).clean_spectrum(:)';
    end

    % Edge case: single node (cannot build a linkage tree)
    if num_nodes == 1
        if show_plot
            warning('clusterExtractedSpectra:SingleNode', ...
                'Fewer than 2 nodes provided; skipping clustering plot.');
        end
        cluster_info = struct( ...
            'spectra_matrix', spectra_matrix, ...
            'corr_matrix', 1, ...
            'clustering_tree', [], ...
            'cluster_order', 1, ...
            'linkage_method', linkage_method);
        return;
    end

    % 2. Calculate Pairwise Correlation Matrix
    corr_matrix = corr(spectra_matrix');
    corr_matrix(isnan(corr_matrix)) = 0;

    % 3. Calculate Hierarchical Linkage
    dist_metric = pdist(spectra_matrix, 'correlation');
    clustering_tree = linkage(dist_metric, linkage_method);


    % 4. Visualization (or silent leaf-order extraction)
    if show_plot
        fig = figure('Name', 'Extracted Spectra Hierarchical Correlation Clustering', ...
                     'Position', [150, 150, 1100, 800], 'Color', 'w');

        % Left Panel: Dendrogram
        ax1 = subplot(1, 5, 1);
        [H, ~, cluster_order] = dendrogram(clustering_tree, 0, 'Orientation', 'left');
        set(H, 'LineWidth', 1.5, 'Color', 'k');
        set(ax1, 'YDir', 'reverse', 'YTick', []);
        title(ax1, 'Dendrogram', 'FontSize', 12, 'FontWeight', 'bold');
        box(ax1, 'off');
        axis(ax1, 'off');

        % Right Panel: Sorted 2D Correlation Matrix
        sorted_corr_matrix = corr_matrix(cluster_order, cluster_order);

        ax2 = subplot(1, 5, 2:5);
        imagesc(ax2, sorted_corr_matrix);
        colormap(ax2, jet(256));
        cbar = colorbar(ax2);
        cbar.Label.String = 'Pearson Correlation (r)';
        cbar.Label.FontWeight = 'bold';
        cbar.Label.FontSize = 12;
        set(ax2, 'CLim', [0 1]);

        title(ax2, '2D Correlation Matrix (Hierarchically Clustered)', ...
              'FontSize', 14, 'FontWeight', 'bold');
        xlabel(ax2, 'Clustered Node Index', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(ax2, 'Clustered Node Index', 'FontSize', 12, 'FontWeight', 'bold');

        set(ax2, 'XTick', 1:num_nodes, 'XTickLabel', cluster_order, 'XTickLabelRotation', 90);
        set(ax2, 'YTick', 1:num_nodes, 'YTickLabel', cluster_order);
        axis(ax2, 'square');
        grid(ax2, 'on');
    else
        tmp_fig = figure('Visible', 'off');
        cleanup_obj = onCleanup(@() delete(tmp_fig)); 
        [~, ~, cluster_order] = dendrogram(clustering_tree, 0);
    end

    % 5. Bundle outputs
    cluster_info = struct( ...
        'spectra_matrix', spectra_matrix, ...
        'corr_matrix', corr_matrix, ...
        'clustering_tree', clustering_tree, ...
        'cluster_order', cluster_order, ...
        'linkage_method', linkage_method);
end
