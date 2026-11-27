function displaySpectra3D(X, ppm, isbackground)
% DISPLAYSPECTRA3D Create an interactive 3D stacked plot of NMR spectra.
%   Visualizes a spectral fraction library (or any set of 1D spectra) as a
%   min-max normalized 3D stacked plot with interactive vertical zoom and
%   clipping toggle controls. Use xlim([min_ppm, max_ppm]) after calling
%   this function to zoom into a specific chemical shift region.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       X            - MxP numeric spectral intensity matrix (rows = fractions,
%                      columns = ppm points).
%       ppm          - 1xP numeric chemical shift vector (ppm).
%       isbackground - (Optional) Logical flag; if true, colors the last
%                      spectrum black (default: false).
%
%   Outputs:
%       None (opens an interactive 3D figure window).

    arguments
        X double
        ppm (1,:) double
        isbackground logical = false
    end

    % Normalize spectral intensities to [0, 1] for 3D display
    X_shrunk = X;
    ppm_shrunk = ppm;
    X_shrunk = (X_shrunk - min(X_shrunk(:))) / (max(X_shrunk(:)) - min(X_shrunk(:)));

    rows = size(X_shrunk, 1);

    % Build 3D coordinate matrices for plot3
    for i = 1:rows
        x(i, :) = ppm_shrunk;
    end
    y = X_shrunk;
    for i = 1:rows
        z(i, :) = zeros(1, length(ppm_shrunk)) + i;
    end

    colors = winter(rows);
    fig = figure();
    ax1 = axes(fig);
    hold on

    % Plot each fraction
    for i = 1:rows - 1
        plot3(x(i, :), z(i, :), y(i, :), 'Color', colors(i, :), 'LineWidth', 0.5)
    end
    if isbackground
        plot3(x(rows, :), z(rows, :), y(rows, :), 'Color', "k", 'LineWidth', 0.5)
    else
        plot3(x(rows, :), z(rows, :), y(rows, :), 'Color', colors(rows, :), 'LineWidth', 0.5)
    end

    hold off
    set(gca, 'XDir', 'rev');
    view(0, 45)
    originalZLim = zlim(ax1);
    ax1.Clipping = "off";
    clipping_variable = false;

    % Create a slider control for vertical zoom
    slider = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 50, 'Value', 1, ...
        'Position', [100, 50, 220, 20], ...
        'Callback', @(src, event) adjustVerticalZoom(ax1, src, originalZLim));

    % Create a text label for the slider
    sliderLabel = uicontrol('Style', 'text', 'Position', [100, 70, 120, 20], ...
        'String', 'Vertical Zoom');

    % Create a button to toggle axis clipping
    selectionButton = uicontrol('Style', 'pushbutton', 'String', 'Change Clipping', ...
        'Position', [100, 100, 120, 20], ...
        'Callback', @(src, event) selectRegion(ax1));

%% Local UI callbacks
    function adjustVerticalZoom(ax, slider, ~)
        % Get the current slider value
        zoomFactor = slider.Value;

        % Adjust the upper limit of the z-axis based on the zoom factor
        zData = get(ax.Children, 'ZData');
        if iscell(zData)
            zData = cell2mat(zData);
        end
        zMax = max(zData(:));
        baseline = -zMax / zoomFactor / 3; % Ensure the baseline is always visible
        zlim(ax, [baseline, zMax / zoomFactor]);
    end

    function selectRegion(ax)
        if clipping_variable
            ax.Clipping = "off";
            clipping_variable = false;
        else
            ax.Clipping = "on";
            clipping_variable = true;
        end
    end
end
