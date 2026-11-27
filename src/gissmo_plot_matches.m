function figs = gissmo_plot_matches(match_out, options)
% GISSMO_PLOT_MATCHES Plot consensus spectra overlaid with their GISSMO matches.
%
%   figs = gissmo_plot_matches(match_out)
%   figs = gissmo_plot_matches(match_out, 'Name', Value)
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       None
%
%   Input:
%       match_out - Struct returned by gissmo_match_spectra.
%
%   Name-Value Options:
%       'TargetIDs'   - Vector of specific Target_IDs to plot (default: all)
%       'MinScore'    - Only plot matches with Score >= MinScore (default: -Inf)
%       'SortByScore' - Plot in descending score order (default: false)
%       'SaveDir'     - Directory to save 300-DPI PNGs (default: "" = do not save)
%
%   Output:
%       figs - Array of generated MATLAB Figure handles.

    arguments
        match_out (1,1) struct
        options.TargetIDs   double {mustBeInteger, mustBePositive} = []
        options.MinScore    (1,1) double  {mustBeReal} = -Inf
        options.SortByScore (1,1) logical = false
        options.SaveDir     (1,1) string  = ""
    end

    if ~all(isfield(match_out, {'results', 'sorted_results', 'ppm_lib', ...
                                'shifted_spectra', 'matched_lib_spectra'}))
        error('gissmo_plot_matches:InvalidInput', ...
            'Input must be a struct returned by gissmo_match_spectra.');
    end

    if options.SortByScore
        tbl = match_out.sorted_results;
    else
        tbl = match_out.results;
    end

    % Filter by TargetIDs and MinScore
    mask = tbl.Score >= options.MinScore;
    if ~isempty(options.TargetIDs)
        if ~isvector(options.TargetIDs)
            error('gissmo_plot_matches:InvalidTargetIDs', 'TargetIDs must be a vector.');
        end
        mask = mask & ismember(tbl.Target_ID, options.TargetIDs(:));
    end
    tbl = tbl(mask, :);

    n_plots = height(tbl);
    if n_plots == 0
        disp('No GISSMO matches met the plotting criteria.');
        figs = gobjects(0);
        return;
    end

    save_figs = strlength(options.SaveDir) > 0;
    if save_figs && ~isfolder(options.SaveDir)
        mkdir(options.SaveDir);
    end

    fprintf('Generating %d GISSMO match figure(s)...\n', n_plots);
    figs = gobjects(n_plots, 1);

    for k = 1:n_plots
        id          = tbl.Target_ID(k);
        match_name  = tbl.Match_Name(k);
        match_idx   = tbl.Match_Index(k);
        score_val   = tbl.Score(k);

        fig_name = sprintf('GISSMO Match: Target %d (%s)', id, match_name);
        fig = figure('Name', fig_name, 'Color', 'w', ...
                     'Units', 'inches', 'Position', [2, 2, 7.5, 5]);

        p1 = plot(match_out.ppm_lib, match_out.matched_lib_spectra(:, id), ...
            'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5);
        hold on;
        p2 = plot(match_out.ppm_lib, match_out.shifted_spectra(:, id), ...
            'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
        hold off;

        ax = gca;
        ax.XDir      = 'reverse';
        ax.FontName  = 'Arial';
        ax.FontSize  = 14;
        ax.LineWidth = 1.5;
        ax.Box       = 'off';
        ax.TickDir   = 'out';
        ax.XColor    = 'k';
        ax.YColor    = 'k';

        xlabel('Chemical Shift (ppm)', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Intensity (a.u.)',     'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

        title_str = sprintf('Target-%d vs. %s (#%d) | Score: %.3f', ...
            id, match_name, match_idx, score_val);
        title(title_str, 'Interpreter', 'none', 'FontName', 'Arial', ...
            'FontSize', 16, 'FontWeight', 'bold');

        legend_str = {sprintf('GISSMO: %s', match_name), 'Consensus Spectrum'};
        lgd = legend([p1, p2], legend_str, 'Location', 'northwest', 'Interpreter', 'none');
        lgd.FontName = 'Arial';
        lgd.FontSize = 14;
        legend('boxoff');

        drawnow;
        figs(k) = fig;

        if save_figs
            safe_name = regexprep(char(match_name), '[^a-zA-Z0-9_-]', '_');
            file_path = fullfile(options.SaveDir, ...
                sprintf('target_%02d_%s.png', id, safe_name));
            print(fig, file_path, '-dpng', '-r300');
        end
    end
end
