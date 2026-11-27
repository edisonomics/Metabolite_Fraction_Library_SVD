function [total_conc, fraction_concs, details] = integrateComponentsSVD(X_raw, ppm, best_state, component_idx, nuclei_analyte, conc_ref, ref_bounds, nuclei_ref, fractions, verbose)
% INTEGRATECOMPONENTSSVD Quantify consensus SVD components via automated qNMR integration.
%   Reconstructs each requested consensus component across the specified
%   fractions (elution_scalar * consensus_spectrum), integrates the full
%   reconstructed component spectrum and the raw internal standard peak
%   (e.g., DSS) via trapezoidal integration, and computes per-fraction and
%   summed molar concentrations using the standard qNMR peak-area ratio.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       X_raw          - MxP raw spectral matrix (rows = fractions, cols = ppm)
%                        prior to reference-peak removal so the internal
%                        standard (e.g., DSS) is intact.
%       ppm            - 1xP numeric chemical shift vector (ascending or descending).
%       best_state     - Struct returned by optimizeClusterCutoff containing
%                        'consensus_spectra' (CxP) and 'consensus_elutions' (CxM).
%       component_idx  - 1xN vector of component row indices in best_state to integrate.
%       nuclei_analyte - 1xN vector (or scalar) of total non-exchangeable protons
%                        represented by each component's full spectrum.
%       conc_ref       - (Optional) Known concentration of the internal reference
%                        standard, e.g., in mM (default: 0.333).
%       ref_bounds     - (Optional) 1x2 [low, high] chemical shift window (ppm)
%                        for integrating the reference peak (default: [-0.02, 0.02]).
%       nuclei_ref     - (Optional) Number of protons represented by the reference
%                        peak (default: 9 for DSS).
%       fractions      - (Optional) Vector of fraction row indices to integrate
%                        and sum over (default: 1:size(X_raw, 1)).
%       verbose        - (Optional) Logical flag to print a summary table to the
%                        command window (default: true).
%
%   Outputs:
%       total_conc     - 1xN vector of summed concentrations across 'fractions'
%                        for each requested component (same units as conc_ref).
%       fraction_concs - FxN matrix of per-fraction concentrations, where
%                        F = numel(fractions) and N = numel(component_idx).
%       details        - Struct with diagnostic outputs:
%                        'ref_integrals'       (Fx1 reference peak integrals),
%                        'component_integrals' (FxN reconstructed component integrals),
%                        'ref_range'           (index vector for reference peak),
%                        'fractions'           (row vector of fraction indices used),
%                        'component_idx'       (row vector of component indices),
%                        'nuclei_analyte'      (1xN proton counts used).

    if nargin < 6 || isempty(conc_ref)
        conc_ref = 0.333;
    end
    if nargin < 7 || isempty(ref_bounds)
        ref_bounds = [-0.02, 0.02];
    end
    if nargin < 8 || isempty(nuclei_ref)
        nuclei_ref = 9;
    end
    if nargin < 9 || isempty(fractions)
        fractions = 1:size(X_raw, 1);
    end
    if nargin < 10 || isempty(verbose)
        verbose = true;
    end

    % --- Input Validation ---
    if isempty(fieldnames(best_state)) || ...
            ~isfield(best_state, 'consensus_spectra') || ...
            ~isfield(best_state, 'consensus_elutions') || ...
            isempty(best_state.consensus_spectra)
        error('integrateComponentsSVD:EmptyState', ...
            'best_state must contain non-empty ''consensus_spectra'' and ''consensus_elutions''.');
    end

    component_idx  = component_idx(:)';
    fractions      = fractions(:)';
    num_comp       = numel(component_idx);
    num_frac       = numel(fractions);

    if isscalar(nuclei_analyte)
        nuclei_analyte = repmat(nuclei_analyte, 1, num_comp);
    else
        nuclei_analyte = nuclei_analyte(:)';
        if numel(nuclei_analyte) ~= num_comp
            error('integrateComponentsSVD:DimensionMismatch', ...
                'nuclei_analyte must be a scalar or have the same length as component_idx (%d).', ...
                num_comp);
        end
    end

    if numel(ref_bounds) ~= 2
        error('integrateComponentsSVD:InvalidRefBounds', ...
            'ref_bounds must be a 2-element vector [low_ppm, high_ppm].');
    end

    num_total_comp = size(best_state.consensus_spectra, 1);
    if any(component_idx < 1) || any(component_idx > num_total_comp)
        error('integrateComponentsSVD:InvalidComponentIndex', ...
            'component_idx values must be between 1 and %d.', num_total_comp);
    end

    if size(X_raw, 2) ~= numel(ppm)
        error('integrateComponentsSVD:PPMSizeMismatch', ...
            'Number of columns in X_raw (%d) must match length of ppm (%d).', ...
            size(X_raw, 2), numel(ppm));
    end

    if size(best_state.consensus_spectra, 2) ~= numel(ppm)
        error('integrateComponentsSVD:SpectraSizeMismatch', ...
            'Number of columns in consensus_spectra (%d) must match length of ppm (%d).', ...
            size(best_state.consensus_spectra, 2), numel(ppm));
    end

    % --- Locate Reference Peak Integration Indices ---
    [~, ref_idx1] = min(abs(ppm - ref_bounds(1)));
    [~, ref_idx2] = min(abs(ppm - ref_bounds(2)));
    ref_range = min(ref_idx1, ref_idx2) : max(ref_idx1, ref_idx2);

    % --- Preallocate Outputs ---
    ref_integrals       = zeros(num_frac, 1);
    component_integrals = zeros(num_frac, num_comp);
    fraction_concs      = zeros(num_frac, num_comp);

    % --- Integrate Reference Peak per Fraction ---
    zero_ref_count = 0;
    for i = 1:num_frac
        row_idx = fractions(i);
        ref_integrals(i) = trapz(X_raw(row_idx, ref_range));
        if ref_integrals(i) == 0
            zero_ref_count = zero_ref_count + 1;
        end
    end

    if zero_ref_count > 0
        warning('integrateComponentsSVD:ZeroReferenceIntegral', ...
            ['Reference peak integral was zero in %d / %d fractions. ' ...
             'Ensure X_raw is the unpreprocessed spectral matrix containing the reference peak.'], ...
            zero_ref_count, num_frac);
    end

    % --- Reconstruct & Integrate Each Component Across Fractions ---
    for k = 1:num_comp
        comp_row = component_idx(k);
        spectrum_vec = best_state.consensus_spectra(comp_row, :);
        n_protons = nuclei_analyte(k);

        for i = 1:num_frac
            row_idx = fractions(i);
            elution_val = best_state.consensus_elutions(comp_row, row_idx);
            reconstructed_signal = elution_val * spectrum_vec;

            int_met = trapz(reconstructed_signal);
            component_integrals(i, k) = int_met;

            if ref_integrals(i) == 0
                fraction_concs(i, k) = NaN;
            else
                fraction_concs(i, k) = ratio_calc( ...
                    int_met, ref_integrals(i), n_protons, conc_ref, nuclei_ref);
            end
        end
    end

    % --- Sum Per-Fraction Concentrations ---
    total_conc = sum(fraction_concs, 1);

    % --- Pack Diagnostic Struct ---
    details = struct();
    details.ref_integrals       = ref_integrals;
    details.component_integrals = component_integrals;
    details.ref_range           = ref_range;
    details.fractions           = fractions;
    details.component_idx       = component_idx;
    details.nuclei_analyte      = nuclei_analyte;

    % --- Optional Summary Output ---
    if verbose
        fprintf('\n--- Automated SVD Integration Summary (%d components, fractions %d:%d) ---\n', ...
            num_comp, fractions(1), fractions(end));
        fprintf('Reference standard: %.4g (units), %d protons, window [%.4f, %.4f] ppm\n', ...
            conc_ref, nuclei_ref, ref_bounds(1), ref_bounds(2));
        fprintf('%-12s %-10s %-16s\n', 'Component', 'Protons', 'Total Conc');
        for k = 1:num_comp
            fprintf('%-12d %-10g %-16.6g\n', component_idx(k), nuclei_analyte(k), total_conc(k));
        end
        fprintf('\n');
    end
end

function conc_x = ratio_calc(integral_analyte, integral_ref, nuclei_analyte, conc_ref, nuclei_ref)
% RATIO_CALC Compute analyte molar concentration from qNMR peak-area ratio.
    conc_x = (integral_analyte / integral_ref) * (nuclei_ref / nuclei_analyte) * conc_ref;
end
