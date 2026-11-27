function [X_small, ppm_small] = downsampleData(X, ppm, num_points)
% DOWNSAMPLEDATA Downsamples spectral data matrix and ppm vector via linear interpolation.
%
% External Functions Required:
%   None
%
% Toolboxes Required:
%   None
%
% Inputs:
%   X          - Data matrix (rows = samples, columns = ppm points)
%   ppm        - Original ppm vector
%   num_points - Desired number of data points in the downsampled ppm vector
%
% Outputs:
%   X_small    - Downsampled data matrix
%   ppm_small  - Downsampled ppm vector
    if num_points >= length(ppm)
        error('num_points (%d) must be smaller than the current ppm length (%d).', ...
              num_points, length(ppm));
    end

    ppm_small = linspace(ppm(1), ppm(end), num_points);
    X_small = zeros(size(X, 1), num_points);

    for i = 1:size(X, 1)
        X_small(i, :) = interp1(ppm, X(i, :), ppm_small, 'linear', 'extrap');
    end
end

