function [nodes_qc, spec_correlations_filt, spec_correlations, fig] = filterComponentsQC(lib, nodes, corr_threshold, show_plot, target_line)
% FILTERCOMPONENTSQC Remove extracted components with poor target spectrum correlation.
%   Recalculates each node's target covariance spectrum using
%   computeCovarianceTarget, computes the Pearson correlation against its
%   extracted clean_spectrum, optionally plots correlation diagnostics, and
%   removes nodes falling below corr_threshold or yielding NaN correlation.
%
%   External Functions Required:
%       computeCovarianceTarget
%
%   Toolboxes Required:
%       Statistics and Machine Learning Toolbox (zscore, via computeCovarianceTarget)
%
%   Inputs:
%       lib            - MxP spectral data matrix (rows = fractions, columns = ppm).
%       nodes          - 1xN struct array from extractComponentsSVD containing
%                        fields 'fracs', 'peaks', and 'clean_spectrum'.
%       corr_threshold - (Optional) Minimum Pearson correlation required to
%                        retain a component (default: 0.2).
%       show_plot      - (Optional) Logical flag to display scatter and
%                        histogram QC plots (default: true).
%       target_line    - (Optional) Reference correlation threshold line
%                        plotted on the scatter panel (default: 0.8).
%
%   Outputs:
%       nodes_qc               - Filtered struct array with poor/NaN nodes removed.
%       spec_correlations_filt - Column vector of correlations for retained nodes.
%       spec_correlations      - Nx1 vector of pre-filter correlations for all
%                                input nodes (NaN for zero-variance edge cases).
%       fig                    - Figure handle if show_plot is true; otherwise [].

    % Handle optional defaults
    if nargin < 3 || isempty(corr_threshold)
        corr_threshold = 0.2;
    end
    if nargin < 4 || isempty(show_plot)
        show_plot = true;
    end
    if nargin < 5 || isempty(target_line)
        target_line = 0.8;
    end

    % Preallocate arrays
    num_nodes = length(nodes);
    spec_correlations = zeros(num_nodes, 1);

    for i = 1:num_nodes
        % 1. Recalculate the Target Covariance Spectrum (specIn)
        matrix_loop = lib(nodes(i).fracs(1):nodes(i).fracs(end), :);

        raw_peak = nodes(i).peaks(:);
        Scl = sum(raw_peak.^2);
        peakIn = raw_peak ./ sqrt(Scl);

        [specIn, ~] = computeCovarianceTarget(matrix_loop, peakIn);

        % 2. Extract the cleaned spectrum
        clean_spec = nodes(i).clean_spectrum(:);
        specIn = specIn(:);

        % 3. Calculate Correlation between Target and Cleaned Spectrum
        R = corrcoef(specIn, clean_spec);

        if size(R, 1) > 1 && ~isnan(R(1, 2))
            spec_correlations(i) = R(1, 2);
        else
            spec_correlations(i) = NaN; % Handle flatline/zero-variance edge cases
        end
    end

    % 4. Visualization
    fig = [];
    if show_plot
        fig = figure('Name', 'Spectral Correlation Analysis', 'Position', [100, 100, 800, 600]);

        % Panel 1: Scatter plot of correlation by node
        subplot(2, 1, 1);
        scatter(1:num_nodes, spec_correlations, 35, 'filled', 'MarkerFaceColor', '#7E2F8E');
        yline(target_line, 'r--', sprintf('Target Threshold (%g)', target_line), ...
            'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);
        title('Similarity: Clean Spectrum vs. Target Covariance Spectrum');
        xlabel('Node Index');
        ylabel('Correlation Coefficient (r)');
        ylim([min([0, min(spec_correlations) - 0.1]) 1.1]);
        grid on;

        % Panel 2: Histogram of correlation distribution
        subplot(2, 1, 2);
        histogram(spec_correlations, 30, 'FaceColor', '#EDB120', 'EdgeColor', 'w');
        title('Distribution of Spectral Correlations');
        xlabel('Correlation Coefficient (r)');
        ylabel('Frequency (Number of Nodes)');
        grid on;
    end

    % 5. Remove any peaks that have a correlation below threshold or NaN
    keep_mask = spec_correlations >= corr_threshold & ~isnan(spec_correlations);
    nodes_qc = nodes(keep_mask);
    spec_correlations_filt = spec_correlations(keep_mask);

end
