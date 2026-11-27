function [wt, wt2, converged, iter] = refineWeightsALS(U, Vsub, peakIn, specIn, row, Amax, wt_init)
% REFINEWEIGHTSALS Refine SVD elution and spectral weights via alternating least squares.
%   Alternates between frequency-domain spectral optimization and time-domain
%   elution optimization with Tikhonov coupling until the relative Frobenius
%   norm change in elution weights falls below tolerance or the iteration
%   limit is reached.
%
%   External Functions Required:
%       private/optimizeSpectra, private/optimizeElution
%
%   Toolboxes Required:
%       Optimization Toolbox (via optimizeSpectra and optimizeElution)
%
%   Inputs:
%       U       - PxR left singular vectors from economy SVD.
%       Vsub    - RxN transposed right singular vectors from economy SVD.
%       peakIn  - Px1 L2-normalized target elution profile.
%       specIn  - 1xN L2-normalized target covariance spectrum.
%       row     - Index of maximum intensity in peakIn.
%       Amax    - Minimum required elution amplitude at index row.
%       wt_init - Rx1 initial elution weights (e.g., from optimizeElution).
%
%   Outputs:
%       wt        - Rx1 refined elution weights.
%       wt2       - Rx1 refined spectral weights.
%       converged - Logical flag indicating whether tolerance was met.
%       iter      - Total number of ALS iterations executed.

    % Fixed ALS constants
    SclS = 0.1;
    max_iter = 10;
    tol = 1e-2;

    iter = 0;
    converged = false;
    wt = wt_init;
    wt2 = wt_init;
    wt_old = wt;

    while iter < max_iter && ~converged
        iter = iter + 1;

        % Step 1: Optimize spectra (frequency domain) -> updates wt2
        wt2 = optimizeSpectra(Vsub, specIn, wt, SclS);

        % Step 2: Optimize elution (time domain) -> updates wt
        wt = optimizeElution(U, peakIn, row, Amax, wt2, SclS);

        % Relative Frobenius change with zero-norm guard
        old_norm = norm(wt_old, 'fro');
        diff_norm = norm(wt - wt_old, 'fro');
        if old_norm > 0
            rel_change = diff_norm / old_norm;
        elseif diff_norm == 0
            rel_change = 0;
        else
            rel_change = Inf;
        end

        if rel_change < tol
            converged = true;
        end

        wt_old = wt;
    end
end
