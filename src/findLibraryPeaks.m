function peaks_struct = findLibraryPeaks(lib, ppm, noise_threshold)
% FINDLIBRARYPEAKS Detects spectral peaks above a noise threshold across all fractions.
%
% Finds local maxima exceeding noise_threshold in each row of the spectral
% library. Choose a threshold high enough to avoid picking up baseline noise
% streaks across fractions (e.g., ~2.5x the baseline noise level).
%
% External Functions Required:
%   None
%
% Toolboxes Required:
%   Signal Processing Toolbox (findpeaks)
%
% Inputs:
%   lib             - Library data matrix (rows = fractions, columns = ppm points)
%   ppm             - Chemical shift vector (ppm)
%   noise_threshold - (Optional) Minimum peak height threshold (default: 0.05)
%
% Outputs:
%   peaks_struct    - 1xN struct array with fields: peak_locs, peaks, ppm

    if nargin < 3
        noise_threshold = 0.05;
    end

    num_fracs = size(lib, 1);
    peaks_struct = struct('peak_locs', cell(1, num_fracs), ...
                          'peaks',     cell(1, num_fracs), ...
                          'ppm',       cell(1, num_fracs));

    for i = 1:num_fracs
        [peaks, locs] = findpeaks(lib(i, :), "MinPeakHeight", noise_threshold);
        peaks_struct(i).peak_locs = locs;
        peaks_struct(i).peaks     = peaks;
        peaks_struct(i).ppm       = reshape(ppm(locs), size(locs));
    end
end
