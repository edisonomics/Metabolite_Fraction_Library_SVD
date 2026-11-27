function [specIn, corr_spectrum] = computeCovarianceTarget(matrix_loop, peakIn)
% COMPUTECOVARIANCETARGET Compute gated, L2-normalized single-target covariance spectrum.
%   Calculates the sample covariance and Pearson correlation between a target
%   elution profile and each chemical shift column of a spectral sub-matrix,
%   gates out poorly correlated or low-intensity regions, and L2-normalizes
%   the resulting target spectrum.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       Statistics and Machine Learning Toolbox (zscore)
%
%   Inputs:
%       matrix_loop - PxN spectral sub-matrix for a node's fractional span.
%       peakIn      - Px1 L2-normalized target elution profile.
%
%   Outputs:
%       specIn        - 1xN gated and L2-normalized target covariance spectrum.
%       corr_spectrum - 1xN Pearson correlation spectrum against peakIn.

    % Fixed gating constants
    corr_gate = 0.97;
    floor_ratio = 0.1;

    peakIn = peakIn(:);
    [n_points, n_ppm] = size(matrix_loop);

    % Guard against degenerate single-fraction inputs
    if n_points <= 1
        specIn = zeros(1, n_ppm);
        corr_spectrum = zeros(1, n_ppm);
        return;
    end

    % Sample covariance (using implicit expansion instead of repmat) and correlation
    covar = (peakIn - mean(peakIn))' * (matrix_loop - mean(matrix_loop, 1)) ./ (n_points - 1);
    corr_spectrum = (zscore(peakIn)' * zscore(matrix_loop)) ./ (n_points - 1);

    % Correlation gate, NaN cleanup, and relative intensity floor
    covar(corr_spectrum < corr_gate) = 0;
    covar(isnan(covar)) = 0;
    covar_max = max(covar);
    if covar_max > 0
        covar(covar < (floor_ratio * covar_max)) = 0;
    end

    % L2 normalization with zero-scale guard
    Scl_spec = sum(covar.^2);
    if Scl_spec > 0
        specIn = covar ./ sqrt(Scl_spec);
    else
        specIn = covar;
    end
end
