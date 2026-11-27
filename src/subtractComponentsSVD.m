function [lib_residual, lib_delete, metrics, fig] = subtractComponentsSVD(lib, ppm, best_state, component_idx, method, ppm_window, frac_idx, show_plot, verbose)
% SUBTRACTCOMPONENTSSVD Subtract selected consensus SVD components from a spectral library.
%   Reconstructs the spectral contribution of a specified subset of consensus
%   components (elutions' * spectra), subtracts or masks them from the input
%   library using the selected method, calculates Frobenius-norm variance
%   explained/remaining metrics over the requested fractions, and optionally
%   visualizes the residual library in 3D via displaySpectra3D.
%
%   External Functions Required:
%       displaySpectra3D (only when show_plot is true)
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       lib           - MxP spectral library matrix (rows = fractions, cols = ppm).
%       ppm           - 1xP numeric chemical shift vector.
%       best_state    - Struct returned by optimizeClusterCutoff containing
%                       'consensus_spectra' (CxP) and 'consensus_elutions' (CxM).
%       component_idx - (Optional) Vector of component row indices in best_state
%                       to subtract (default: 1:size(best_state.consensus_spectra, 1)).
%       method        - (Optional) Subtraction method (default: 'zero_mask'):
%                       'zero_mask'      - Zeroes lib wherever lib_delete ~= 0.
%                       'window_mask'    - Zeroes lib within +/- ppm_window around
%                                          every non-zero point in lib_delete per fraction.
%                       'linear'         - Direct residual: lib - lib_delete.
%                       'linear_clipped' - Direct residual with negatives clipped to 0.
%       ppm_window    - (Optional) Half-width in ppm for 'window_mask' mode
%                       (default: 0.03).
%       frac_idx      - (Optional) Vector of fraction row indices over which to
%                       compute variance metrics (default: 1:size(lib, 1)).
%       show_plot     - (Optional) Logical flag to display the 3D residual via
%                       displaySpectra3D (default: true).
%       verbose       - (Optional) Logical flag to print variance metrics to the
%                       command window (default: true).
%
%   Outputs:
%       lib_residual  - MxP residual spectral library after subtraction.
%       lib_delete    - MxP reconstructed contribution of the subtracted components.
%       metrics       - Struct containing variance metrics over frac_idx:
%                       'frac_idx'                (fractions evaluated),
%                       'method'                  (subtraction method used),
%                       'variance_removed_pct'    (% variance removed by selected method),
%                       'variance_remaining_pct'  (% variance remaining in lib_residual),
%                       'variance_linear_pct'     (% variance explained by linear subtraction),
%                       'variance_full_model_pct' (% variance explained by all components).
%       fig           - Figure handle from displaySpectra3D (or [] if show_plot is false).

    if nargin < 4 || isempty(component_idx)
        component_idx = 1:size(best_state.consensus_spectra, 1);
    end
    if nargin < 5 || isempty(method)
        method = 'zero_mask';
    end
    if nargin < 6 || isempty(ppm_window)
        ppm_window = 0.03;
    end
    if nargin < 7 || isempty(frac_idx)
        frac_idx = 1:size(lib, 1);
    end
    if nargin < 8 || isempty(show_plot)
        show_plot = true;
    end
    if nargin < 9 || isempty(verbose)
        verbose = true;
    end

    % --- Input Validation ---
    if isempty(fieldnames(best_state)) || ...
            ~isfield(best_state, 'consensus_spectra') || ...
            ~isfield(best_state, 'consensus_elutions') || ...
            isempty(best_state.consensus_spectra)
        error('subtractComponentsSVD:EmptyState', ...
            'best_state must contain non-empty ''consensus_spectra'' and ''consensus_elutions''.');
    end

    component_idx = component_idx(:)';
    frac_idx      = frac_idx(:)';
    method_str    = lower(char(method));

    num_total_comp = size(best_state.consensus_spectra, 1);
    if any(component_idx < 1) || any(component_idx > num_total_comp)
        error('subtractComponentsSVD:InvalidComponentIndex', ...
            'component_idx values must be between 1 and %d.', num_total_comp);
    end

    if size(lib, 2) ~= numel(ppm)
        error('subtractComponentsSVD:PPMSizeMismatch', ...
            'Number of columns in lib (%d) must match length of ppm (%d).', ...
            size(lib, 2), numel(ppm));
    end

    % --- 1. Reconstruct Selected Components ---
    cons_spec  = best_state.consensus_spectra(component_idx, :);
    cons_elut  = best_state.consensus_elutions(component_idx, :);
    lib_delete = cons_elut' * cons_spec;

    % --- 2. Apply Subtraction Method ---
    switch method_str
        case 'zero_mask'
            lib_residual = lib;
            lib_residual(lib_delete ~= 0) = 0;

        case 'window_mask'
            lib_residual = lib;
            for r = 1:size(lib_delete, 1)
                indices = find(lib_delete(r, :) ~= 0);
                for i = 1:numel(indices)
                    center_idx = indices(i);
                    region = abs(ppm - ppm(center_idx)) <= ppm_window;
                    lib_residual(r, region) = 0;
                end
            end

        case 'linear'
            lib_residual = lib - lib_delete;

        case 'linear_clipped'
            lib_residual = lib - lib_delete;
            lib_residual(lib_residual < 0) = 0;

        otherwise
            error('subtractComponentsSVD:InvalidMethod', ...
                ['Unsupported method ''%s''. Choose ''zero_mask'', ''window_mask'', ' ...
                 '''linear'', or ''linear_clipped''.'], method_str);
    end

    % --- 3. Calculate Variance Metrics on Indexed Fractions ---
    lib_full_recon = best_state.consensus_elutions' * best_state.consensus_spectra;

    lib_sub       = lib(frac_idx, :);
    lib_res_sub   = lib_residual(frac_idx, :);
    lib_del_sub   = lib_delete(frac_idx, :);
    lib_recon_sub = lib_full_recon(frac_idx, :);

    lib_norm_sq = norm(lib_sub, 'fro')^2;

    if lib_norm_sq > 0
        residual_norm_sq       = norm(lib_res_sub, 'fro')^2;
        variance_removed_pct   = (1 - (residual_norm_sq / lib_norm_sq)) * 100;
        variance_remaining_pct = (residual_norm_sq / lib_norm_sq) * 100;

        before_sub_res_norm_sq  = norm(lib_sub - lib_recon_sub, 'fro')^2;
        variance_full_model_pct = (1 - (before_sub_res_norm_sq / lib_norm_sq)) * 100;

        linear_res_norm_sq  = norm(lib_sub - lib_del_sub, 'fro')^2;
        variance_linear_pct = (1 - (linear_res_norm_sq / lib_norm_sq)) * 100;
    else
        variance_removed_pct    = 100;
        variance_remaining_pct  = 0;
        variance_full_model_pct = 100;
        variance_linear_pct     = 100;
    end

    metrics = struct();
    metrics.frac_idx                = frac_idx;
    metrics.method                  = method_str;
    metrics.variance_removed_pct    = variance_removed_pct;
    metrics.variance_remaining_pct  = variance_remaining_pct;
    metrics.variance_linear_pct     = variance_linear_pct;
    metrics.variance_full_model_pct = variance_full_model_pct;

    % --- 4. Optional Console Report ---
    if verbose
        fprintf('\n--- Variance Metrics (Fractions %d to %d, %d components, method: %s) ---\n', ...
            frac_idx(1), frac_idx(end), numel(component_idx), method_str);
        fprintf('Total Variance Explained (All Components, Before Subtraction): %.2f%%\n', ...
            variance_full_model_pct);
        fprintf('Variance Removed (%s): %.2f%%\n', method_str, variance_removed_pct);
        fprintf('Variance Remaining in Library (%s): %.2f%%\n', method_str, variance_remaining_pct);
        fprintf('Variance Explained (Linear Subtraction): %.2f%%\n', variance_linear_pct);
        fprintf('Variance Remaining (Linear Subtraction): %.2f%%\n\n', 100 - variance_linear_pct);
    end

    % --- 5. Optional 3D Visualization ---
    if show_plot
        displaySpectra3D(lib_residual, ppm);
        fig = gcf;
    else
        fig = [];
    end
end
