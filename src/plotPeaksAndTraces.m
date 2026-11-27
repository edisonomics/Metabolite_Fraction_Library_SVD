function fig = plotPeaksAndTraces(lib, ppm, overlay_struct)
% PLOTPEAKSANDTRACES Plots a 3D spectral library with optional peak markers or trace overlays.
%
% Automatically detects whether overlay_struct contains detected peaks
% (from findLibraryPeaks) or traced paths (from tracePeaks or
% removeRedundantTraces) and renders the matching overlay.
%
% External Functions Required:
%   None
%
% Toolboxes Required:
%   None
%
% Inputs:
%   lib            - Library data matrix (rows = fractions, columns = ppm points)
%   ppm            - Chemical shift vector (ppm)
%   overlay_struct - (Optional) peaks_struct (magenta circles) or nodes struct (trace lines)
%
% Outputs:
%   fig            - Handle to the generated figure

    num_fracs = size(lib, 1);
    num_pts   = size(lib, 2);

    % Build 3D coordinate grids matching the library dimensions
    x = repmat(reshape(ppm, 1, []), num_fracs, 1);
    y = lib;
    z = repmat((1:num_fracs)', 1, num_pts);

    colors = winter(num_fracs);
    has_overlay = (nargin >= 3 && ~isempty(overlay_struct));

    % Use a 900x750 figure window for 3D trace overlays and default sizing for peak markers
    if has_overlay && isfield(overlay_struct, 'fracs')
        fig = figure('Position', [10 10 900 750]);
    else
        fig = figure;
    end
    ax = axes(fig);
    hold(ax, 'on');

    if has_overlay && isfield(overlay_struct, 'peak_locs')
        % Plot spectra with magenta peak markers
        for i = 1:num_fracs
            plot3(ax, x(i, :), z(i, :), y(i, :), 'Color', colors(i, :), 'LineWidth', 0.5);
            locs = overlay_struct(i).peak_locs;
            plot3(ax, x(i, locs), z(i, locs), y(i, locs), "o", "Color", "m");
        end
    else
        % Plot base spectra
        for i = 1:num_fracs
            plot3(ax, x(i, :), z(i, :), y(i, :), 'Color', colors(i, :), 'LineWidth', 0.5);
        end

        % Overlay 3D traces if a nodes struct was provided
        if has_overlay && isfield(overlay_struct, 'fracs')
            for i = 1:length(overlay_struct)
                plot3(ax, ppm(overlay_struct(i).locs), ...
                          overlay_struct(i).fracs, ...
                          overlay_struct(i).peaks, 'LineWidth', 1.5);
            end
        end
    end

    hold(ax, 'off');
    set(ax, 'XDir', 'rev');
    view(ax, [0 70]);
    ax.Clipping = "off";
end
