function wt = optimizeElution(Usub, peakIn, row, Amax, wt2, SclS)
% OPTIMIZEELUTION Optimize SVD elution weights with optional Tikhonov coupling.
%   Solves a constrained nonlinear optimization problem in SVD space to fit
%   target elution profile(s) with weighted baseline error, soft
%   non-negativity penalty, peak-amplitude floor, L2 normalization, and
%   baseline AUC ceiling constraints.
%
%   Operates in two modes:
%     - Initializer mode (wt2 omitted or empty): constructs a projected
%       non-negative initial guess and omits Tikhonov regularization.
%     - Refinement mode (wt2 and SclS provided): warm-starts from wt2 and
%       adds a Tikhonov tether penalizing deviation from wt2.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       Optimization Toolbox (optimvar, optimproblem, optimoptions, solve, fmincon)
%
%   Inputs:
%       Usub   - PxR left singular vector basis matrix.
%       peakIn - PxK matrix of L2-normalized target elution profile(s).
%       row    - 1xK (or Kx1) peak apex row indices for each component.
%       Amax   - 1xK (or Kx1) minimum allowed amplitude at row(j).
%       wt2    - (Optional) RxK spectral weights for warm-start and Tikhonov
%                tethering. Omit or pass [] for initializer mode.
%       SclS   - (Optional) Scalar Tikhonov weight (required when wt2 is non-empty).
%
%   Outputs:
%       wt     - RxK optimized elution weights in SVD space.

    is_init_mode = (nargin < 5) || isempty(wt2);
    if ~is_init_mode && (nargin < 6 || isempty(SclS))
        error('optimizeElution:MissingSclS', ...
            'SclS must be provided when wt2 is supplied.');
    end

    % Fixed solver weights, tolerances, and fmincon options
    baseline_weight = 100;
    epsilon = 1e-6;
    neg_weight = 5000;
    norm_tol = 0.005;
    auc_factor = 0.005;

    n_svd = size(Usub, 2);
    n_phys = size(peakIn, 2);
    wt = zeros(n_svd, n_phys);

    options = optimoptions('fmincon', ...
        'MaxIterations', 10000, ...
        'MaxFunctionEvaluations', 5000, ...
        'Display', 'off', ...
        'Algorithm', 'interior-point');

    for j = 1:n_phys
        a = optimvar('a', 1, n_svd);
        target_profile = peakIn(:, j)';

        if is_init_mode
            % Project target onto SVD space, enforce non-negativity and L2 = 1
            raw_guess = target_profile * Usub;
            physical_guess = raw_guess * Usub';
            physical_guess(physical_guess < 0) = 0;

            Scl_guess = sum(physical_guess.^2);
            if Scl_guess > 0
                physical_guess = physical_guess ./ sqrt(Scl_guess);
            end
            x0.a = physical_guess * Usub;
        else
            x0.a = wt2(:, j)';
        end

        prob = optimproblem;
        pred_profile = a * Usub';

        % Weighted baseline SSE
        baseline_mask = (target_profile == 0);
        error_weights = ones(1, length(target_profile));
        error_weights(baseline_mask) = baseline_weight;
        global_error = sum(error_weights .* (pred_profile - target_profile).^2);

        % Soft non-negativity penalty
        smooth_abs = sqrt(pred_profile.^2 + epsilon);
        negative_violations = (-pred_profile + smooth_abs) / 2;
        neg_penalty = neg_weight * sum(negative_violations.^2);

        if is_init_mode
            prob.Objective = global_error + neg_penalty;
        else
            tikhonov_penalty = SclS * sum((a - wt2(:, j)').^2);
            prob.Objective = global_error + tikhonov_penalty + neg_penalty;
        end

        % Constraints
        prob.Constraints.cons1 = a * Usub(row(j), :)' >= Amax(j);
        prob.Constraints.cons2 = (dot(pred_profile, pred_profile) - 1)^2 <= norm_tol;

        max_allowed_area = auc_factor * max(target_profile) * sum(baseline_mask);
        prob.Constraints.cons4 = sum(pred_profile(baseline_mask)) <= max_allowed_area;


        % Solve with safe fallback on empty/failed solution
        sol = solve(prob, x0, 'Options', options);
        if isempty(sol) || ~isfield(sol, 'a') || isempty(sol.a)
            wt(:, j) = x0.a';
        else
            wt(:, j) = sol.a';
        end
    end
end
