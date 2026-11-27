function nodes_filtered = removeRedundantTraces(nodes, corr_threshold)
% REMOVEREDUNDANTTRACES Removes highly correlated duplicate traces sharing the same fraction span.
%
% Groups traces by their exact fractional span, sorts each group by maximum
% absolute intensity (descending), and removes lower-intensity traces whose
% Pearson correlation with a stronger trace meets or exceeds corr_threshold.
%
% External Functions Required:
%   None
%
% Toolboxes Required:
%   Statistics and Machine Learning Toolbox (corr)
%
% Inputs:
%   nodes          - Struct array of traces from tracePeaks (fields: fracs, ppms, locs, peaks)
%   corr_threshold - (Optional) Pearson correlation threshold for removal (default: 0.99)
%
% Outputs:
%   nodes_filtered - Filtered struct array of non-redundant traces

    if nargin < 2
        corr_threshold = 0.99;
    end

    if isempty(nodes)
        nodes_filtered = nodes;
        return;
    end

    % Create a string representation of the exact fractions each trace spans
    frac_strs = cell(1, length(nodes));
    for i = 1:length(nodes)
        frac_strs{i} = mat2str(nodes(i).fracs);
    end

    % Group traces by unique fractional spans
    [unique_spans, ~, ic] = unique(frac_strs);
    indices_to_keep = [];

    for i = 1:length(unique_spans)
        idx = find(ic == i);

        if isscalar(idx)
            indices_to_keep = [indices_to_keep, idx]; %#ok<AGROW>
            continue;
        end

        % Extract peak intensities into a matrix for correlation
        L = length(nodes(idx(1)).fracs);
        P = zeros(length(idx), L);
        for j = 1:length(idx)
            P(j, :) = nodes(idx(j)).peaks;
        end

        % Sort traces by maximum absolute intensity in descending order
        max_ints = max(abs(P), [], 2);
        [~, sort_order] = sort(max_ints, 'descend');
        sorted_idx = idx(sort_order);
        P_sorted = P(sort_order, :);

        % Keep the highest-intensity trace and remove highly correlated weaker ones
        keep_mask = true(length(sorted_idx), 1);
        for m = 1:length(sorted_idx)
            if keep_mask(m) && m < length(sorted_idx)
                v1 = P_sorted(m, :)';
                V_others = P_sorted(m+1:end, :)';

                C = corr(v1, V_others);
                C(isnan(C)) = 0; % Handle flat/zero-variance arrays safely

                corr_idx = find(C >= corr_threshold) + m;
                keep_mask(corr_idx) = false;
            end
        end

        indices_to_keep = [indices_to_keep, reshape(sorted_idx(keep_mask), 1, [])]; %#ok<AGROW>
    end

    nodes_filtered = nodes(indices_to_keep);
end
