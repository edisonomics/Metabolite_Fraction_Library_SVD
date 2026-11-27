function lib = preprocessFractions(X, ppm, noise_region, remove_regions, show_plot)
% PREPROCESSFRACTIONS Replaces unwanted spectral peaks with baseline noise and clips negatives.
%
% Removes peaks that span all fractions (e.g., DSS, formate, water) by
% replacing those regions with synthetic noise matching a peak-free baseline
% region and clips negative values to zero.
%
% External Functions Required:
%   displaySpectra3D (only if show_plot is true)
%
% Toolboxes Required:
%   None
%
% Inputs:
%   X              - Data matrix (rows = fractions/samples, columns = ppm points)
%   ppm            - Chemical shift vector (ppm)
%   noise_region   - 1x2 vector [start, end] defining a peak-free noise region
%   remove_regions - Nx2 matrix where each row [start, end] is a ppm region to replace
%   show_plot      - (Optional) Logical true/false to run displaySpectra3D (default: true)
%
% Outputs:
%   lib            - Processed data matrix

    if nargin < 5
        show_plot = true;
    end

    % Combine all target peak regions into a single logical mask
    target_noise_mask = false(size(ppm));
    for i = 1:size(remove_regions, 1)
        r_min = min(remove_regions(i, :));
        r_max = max(remove_regions(i, :));
        target_noise_mask = target_noise_mask | (ppm >= r_min & ppm <= r_max);
    end

    % Isolate the baseline noise region and compute its mean and standard deviation
    noise_region_mask = (ppm >= min(noise_region) & ppm <= max(noise_region));
    X_noise_subset = X(:, noise_region_mask);
    noise_mean = mean(X_noise_subset(:));
    noise_std = std(X_noise_subset(:));

    % Overwrite the specified peak regions across all rows with matching baseline noise
    lib = X;
    num_total_rows = size(X, 1);
    num_target_cols = sum(target_noise_mask);
    lib(:, target_noise_mask) = noise_mean + noise_std * randn(num_total_rows, num_target_cols);

    % Clip negative values to zero
    lib(lib <= 0) = 0.0;

    % Optionally plot the processed fractions
    if show_plot
        displaySpectra3D(lib, ppm);
    end
end
