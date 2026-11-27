function wt2 = optimizeSpectra(Vsub, specIn, wt, SclS)
% OPTIMIZESPECTRA Optimize SVD spectral weights against target covariance spectra.
%   Solves a constrained nonlinear optimization problem in SVD space combining
%   shape-tracking SSE against a target covariance spectrum, Tikhonov
%   tethering to elution weights, soft non-negativity, targeted baseline/peak
%   sparsity penalties, and an L2 normalization constraint.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       Optimization Toolbox (optimvar, optimproblem, optimoptions, solve, fmincon)
%
%   Inputs:
%       Vsub   - RxN transposed right singular vector basis matrix.
%       specIn - KxN matrix of L2-normalized target covariance spectra.
%       wt     - RxK current elution weights in SVD space.
%       SclS   - Scalar Tikhonov tether weight anchoring wt2 to wt.
%
%   Outputs:
%       wt2    - RxK optimized spectral weights in SVD space.

    % Fixed solver weights, tolerances, and fmincon options
    epsilon = 1e-6;
    neg_weight = 5000;
    baseline_penalty_weight = 5000;
    peak_penalty_weight = 0.1;
    target_weight = 0.15;
    norm_tol = 0.005;

    n_svd = size(wt, 1);
    n_phys = size(wt, 2);
    wt2 = zeros(n_svd, n_phys);

    options = optimoptions('fmincon', ...
        'MaxIterations', 10000, ...
        'MaxFunctionEvaluations', 5000, ...
        'Display', 'off', ...
        'Algorithm', 'interior-point');

    for j = 1:n_phys
        a = optimvar('a', 1, n_svd);
        x0.a = wt(:, j)';

        target_spec = specIn(j, :);

        % Silent dead-component bypass when target spectrum is identically zero
        if sum(abs(target_spec)) == 0
            wt2(:, j) = wt(:, j);
            continue;
        end

        prob = optimproblem;
        fullVec = a * Vsub;

        % 1. Soft non-negativity penalty
        smooth_abs = sqrt(fullVec.^2 + epsilon);
        negative_violations = (-fullVec + smooth_abs) / 2;
        neg_penalty = neg_weight * sum(negative_violations.^2);

        % 2. Targeted sparsity penalty (guarded against empty masks)
        baseline_mask = (target_spec == 0);
        if any(baseline_mask)
            baseline_penalty = baseline_penalty_weight * sum(smooth_abs(baseline_mask));
        else
            baseline_penalty = 0;
        end
        if any(~baseline_mask)
            peak_penalty = peak_penalty_weight * sum(smooth_abs(~baseline_mask));
        else
            peak_penalty = 0;
        end
        sparsity_penalty = baseline_penalty + peak_penalty;

        % 3. Shape-tracking target pull
        target_pull = target_weight * sum((fullVec - target_spec).^2);

        % 4. Tikhonov tether to elution weights
        tether_penalty = SclS * sum((a - wt(:, j)').^2);

        % 5. Combined objective
        prob.Objective = target_pull + tether_penalty + neg_penalty + sparsity_penalty;

        % 6. Normalization constraint
        prob.Constraints.cons1 = (dot(fullVec, fullVec) - 1)^2 <= norm_tol;

        % Solve with safe fallback on empty/failed solution
        sol = solve(prob, x0, 'Options', options);
        if isempty(sol) || ~isfield(sol, 'a') || isempty(sol.a)
            wt2(:, j) = x0.a';
        else
            wt2(:, j) = sol.a';
        end
    end
end
