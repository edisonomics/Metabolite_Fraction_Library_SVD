function Vtrn_clean = filterExtractedSpectrum(Vtrn, corr_spectrum)
% FILTEREXTRACTEDSPECTRUM Apply correlation-modulated intensity filter and edge smoothing.
%   Clips negative intensities, thresholds each spectral point using a dynamic
%   cutoff scaled by its correlation to the target elution profile, applies a
%   morphological Gaussian buffer to hard threshold edges, and L2-normalizes
%   the cleaned spectrum.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       Image Processing Toolbox (imdilate)
%
%   Inputs:
%       Vtrn          - 1xN raw extracted spectrum (e.g., nodes_svd(i).raw_spectrum).
%       corr_spectrum - 1xN correlation spectrum (e.g., nodes_svd(i).corr_spectrum).
%
%   Outputs:
%       Vtrn_clean    - 1xN cleaned and L2-normalized spectrum.

    % Fixed filtering and smoothing constants
    corr_floor = 0.7;
    corr_ceil = 0.9;
    base_cutoff = 0.05;
    bleed_cutoff = 0.15;
    dilate_width = 5;
    smooth_window = 7;

    v_raw = max(0, Vtrn);
    max_peak = max(v_raw);

    if max_peak > 0
        % 1. Scale correlation between 0 and 1
        corr_scaled = (corr_spectrum - corr_floor) ./ (corr_ceil - corr_floor);
        corr_scaled = max(0, min(1, corr_scaled));

        % 2. Compute dynamic percentage and absolute intensity thresholds
        dynamic_pct = bleed_cutoff - corr_scaled .* (bleed_cutoff - base_cutoff);
        threshold_vector = dynamic_pct * max_peak;

        % 3. Apply dynamic intensity zeroing
        v_hard = v_raw;
        v_hard(v_raw < threshold_vector) = 0;

        % 4. Smooth hard cliff edges using a morphological buffer
        is_signal = v_hard > 0;
        cliff_mask = imdilate(is_signal, ones(1, dilate_width)) & ~is_signal;

        v_smoothed = smoothdata(v_hard, 'gaussian', smooth_window);
        v_hard(cliff_mask) = v_smoothed(cliff_mask);

        % 5. Re-normalize with zero-norm guard
        norm_val = norm(v_hard);
        if norm_val > 0
            Vtrn_clean = v_hard / norm_val;
        else
            Vtrn_clean = v_hard;
        end
    else
        Vtrn_clean = v_raw;
    end
end
