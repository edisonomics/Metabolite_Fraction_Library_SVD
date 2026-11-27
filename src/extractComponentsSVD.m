function nodes_out = extractComponentsSVD(lib, nodes, opts)
% EXTRACTCOMPONENTSSVD Extract pure component spectra and elution profiles via SVD-ALS.
%   Iterates over each traced peak node, computes an economy SVD over the
%   node's fractional span, initializes elution weights, builds a
%   single-target covariance spectrum, refines weights via alternating least
%   squares (ALS), projects final spectral/elution components, and applies a
%   dynamic correlation-modulated intensity filter.
%
%   External Functions Required:
%       computeCovarianceTarget, refineWeightsALS, filterExtractedSpectrum,
%       private/optimizeElution, private/optimizeSpectra
%
%   Toolboxes Required:
%       Optimization Toolbox, Statistics and Machine Learning Toolbox,
%       Image Processing Toolbox, Parallel Computing Toolbox (only if
%       opts.use_parallel is true)
%
%   Inputs:
%       lib   - MxN spectral data matrix (rows = fractions, cols = ppm points).
%       nodes - 1xK struct array of traced peaks (e.g., nodes_filtered) with
%               fields 'fracs' and 'peaks'.
%       opts  - (Optional) Struct of runtime options:
%                 .use_parallel - Logical toggle to run with parfor
%                                 (default: false).
%                 .verbose      - Logical toggle for start/completion status
%                                 messages (default: true).
%
%   Outputs:
%       nodes_out - Copy of nodes augmented with fields:
%                     .clean_spectrum    - 1xN filtered, L2-normalized spectrum.
%                     .optimized_elution - Px1 non-negative, L2-normalized elution profile.
%                     .singular_value    - Scalar extracted core singular value.
%                     .raw_spectrum      - 1xN unfiltered extracted spectrum (Vtrn2).
%                     .corr_spectrum     - 1xN correlation spectrum against peakIn.
%                     .als_converged     - Logical ALS convergence flag.
%                     .als_iters         - Number of ALS iterations executed.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if isfield(opts, 'use_parallel') && ~isempty(opts.use_parallel)
        use_parallel = logical(opts.use_parallel);
    else
        use_parallel = false;
    end

    if isfield(opts, 'verbose') && ~isempty(opts.verbose)
        verbose = logical(opts.verbose);
    else
        verbose = true;
    end

    % Fixed algorithm constant
    AmaxScl = 0.7;

    n_nodes = numel(nodes);
    nodes_out = nodes;

    if n_nodes == 0
        return;
    end

    % Pre-extract sliced inputs for clean parfor classification
    frac_starts = zeros(n_nodes, 1);
    frac_ends = zeros(n_nodes, 1);
    raw_peaks = cell(n_nodes, 1);
    for i = 1:n_nodes
        frac_starts(i) = nodes(i).fracs(1);
        frac_ends(i) = nodes(i).fracs(end);
        raw_peaks{i} = nodes(i).peaks(:);
    end

    % Preallocate sliced output containers
    clean_spectra = cell(n_nodes, 1);
    opt_elutions = cell(n_nodes, 1);
    singular_vals = cell(n_nodes, 1);
    raw_spectra = cell(n_nodes, 1);
    corr_spectra = cell(n_nodes, 1);
    als_conv_flags = false(n_nodes, 1);
    als_iter_counts = zeros(n_nodes, 1);

    t_start = tic;
    if verbose
        if use_parallel
            fprintf('Starting SVD-ALS extraction for %d nodes (parallel mode)...\n', n_nodes);
        else
            fprintf('Starting SVD-ALS extraction for %d nodes (serial mode)...\n', n_nodes);
        end
    end

    if use_parallel
        parfor i = 1:n_nodes
            [clean_spectra{i}, opt_elutions{i}, singular_vals{i}, ...
             raw_spectra{i}, corr_spectra{i}, als_conv_flags(i), ...
             als_iter_counts(i)] = processSingleNode( ...
                lib(frac_starts(i):frac_ends(i), :), raw_peaks{i}, AmaxScl);
        end
    else
        for i = 1:n_nodes
            [clean_spectra{i}, opt_elutions{i}, singular_vals{i}, ...
             raw_spectra{i}, corr_spectra{i}, als_conv_flags(i), ...
             als_iter_counts(i)] = processSingleNode( ...
                lib(frac_starts(i):frac_ends(i), :), raw_peaks{i}, AmaxScl);
        end
    end

    % Assemble output struct
    for i = 1:n_nodes
        nodes_out(i).clean_spectrum = clean_spectra{i};
        nodes_out(i).optimized_elution = opt_elutions{i};
        nodes_out(i).singular_value = singular_vals{i};
        nodes_out(i).raw_spectrum = raw_spectra{i};
        nodes_out(i).corr_spectrum = corr_spectra{i};
        nodes_out(i).als_converged = als_conv_flags(i);
        nodes_out(i).als_iters = als_iter_counts(i);
    end

    elapsed_time = toc(t_start);
    if verbose
        fprintf('Completed SVD-ALS extraction for %d nodes in %.2f seconds.\n', ...
            n_nodes, elapsed_time);
    end
end

function [Vtrn4, Utrn2, Strn, Vtrn2, corr_spectrum, converged, iter] = ...
    processSingleNode(matrix_loop, raw_peak, AmaxScl)
% PROCESSSINGLENODE Execute SVD, ALS refinement, projection, and filtering for one node.

    [U, S, V] = svd(matrix_loop, 'econ');
    n_comp = min(size(matrix_loop, 1), size(V, 2));
    Vsub = V(:, 1:n_comp)';

    % L2 normalization of target elution trace with zero-scale guard
    Scl = sum(raw_peak.^2);
    if Scl > 0
        peakIn = raw_peak ./ sqrt(Scl);
    else
        peakIn = raw_peak;
    end

    [~, row] = max(peakIn);
    Amax = peakIn(row) * AmaxScl;

    % Initial elution optimization (initializer mode)
    wt_init = optimizeElution(U, peakIn, row, Amax);

    % Single-target covariance spectrum
    [specIn, corr_spectrum] = computeCovarianceTarget(matrix_loop, peakIn);

    % Alternating least squares refinement
    [wt, wt2, converged, iter] = refineWeightsALS( ...
        U, Vsub, peakIn, specIn, row, Amax, wt_init);

    % 1. Frequency domain extraction
    if norm(wt2) > 0
        Vtrn2 = wt2 \ Vsub;
    else
        Vtrn2 = zeros(1, size(Vsub, 2));
    end

    % 2. Time domain extraction
    Utrn2 = U * wt;
    Utrn2(Utrn2 < 0) = 0;
    norm_val = norm(Utrn2);
    if norm_val > 0
        Utrn2 = Utrn2 / norm_val;
    end

    % 3. Core singular value extraction
    if norm(wt) > 0
        Strn = abs(wt \ (S * wt2));
    else
        Strn = 0;
    end

    % 4. Dynamic correlation-modulated intensity filter
    Vtrn4 = filterExtractedSpectrum(Vtrn2, corr_spectrum);
end
