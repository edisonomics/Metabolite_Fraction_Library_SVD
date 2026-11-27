function out = gissmo_match_spectra(spec_matrix, ppm_target, gissmo_source, options)
% GISSMO_MATCH_SPECTRA Match consensus spectra against a GISSMO NMR library.
%
%   out = gissmo_match_spectra(spec_matrix, ppm_target, gissmo_source)
%   out = gissmo_match_spectra(..., 'Name', Value)
%
%   External Functions Required:
%       gissmo_load_library, gissmo_plot_matches (only if PlotMatches is true)
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       spec_matrix   - [ns x n_ppm] matrix of consensus spectra (one per row).
%       ppm_target    - [1 x n_ppm] chemical shift axis for spec_matrix.
%       gissmo_source - Either a file path to GISSMO_Lib.mat or a struct
%                       returned by gissmo_load_library.
%
%   Name-Value Options:
%       'MinFrac'     - Peak removal threshold as fraction of max peak (0..1, default: 0.0)
%       'ShiftMax'    - Max point shift in either direction (default: 50)
%       'ShiftStep'   - Point step size between tested shifts (default: 5)
%       'SearchMode'  - "global" (uniform shift) or "split" (independent L/R halves)
%                       (default: "global")
%       'PlotMatches' - Automatically call gissmo_plot_matches (default: false)
%       'Verbose'     - Print formatted summary tables to console (default: true)
%
%   Output:
%       out - Struct with fields:
%             .results             - Table ordered by Target_ID
%             .sorted_results      - Table ranked by Score (descending)
%             .ppm_lib             - [1 x pts] GISSMO ppm axis
%             .shifted_spectra     - [pts x ns] shifted & normalized target spectra
%             .matched_lib_spectra - [pts x ns] best-matching normalized GISSMO spectra
%             .options             - Options used for this run

    arguments
        spec_matrix (:,:) double {mustBeReal, mustBeFinite}
        ppm_target  (1,:) double {mustBeReal, mustBeFinite}
        gissmo_source
        options.MinFrac     (1,1) double  {mustBeReal, mustBeFinite, ...
            mustBeGreaterThanOrEqual(options.MinFrac, 0), ...
            mustBeLessThanOrEqual(options.MinFrac, 1)} = 0.0
        options.ShiftMax    (1,1) double  {mustBeInteger, mustBeNonnegative} = 50
        options.ShiftStep   (1,1) double  {mustBeInteger, mustBePositive}    = 5
        options.SearchMode  (1,1) string  {mustBeMember(options.SearchMode, ["global", "split"])} = "global"
        options.PlotMatches (1,1) logical = false
        options.Verbose     (1,1) logical = true
    end

    if size(spec_matrix, 2) ~= numel(ppm_target)
        error('gissmo_match_spectra:DimensionMismatch', ...
            'spec_matrix has %d columns, but ppm_target has %d elements.', ...
            size(spec_matrix, 2), numel(ppm_target));
    end

    if options.Verbose
        fprintf('\n=========================================================================\n');
        fprintf(' MATCHING CONSENSUS SPECTRA AGAINST GISSMO LIBRARY (%s mode)\n', upper(options.SearchMode));
        fprintf('=========================================================================\n');
    end

    % 1. Load or validate GISSMO library
    if isstruct(gissmo_source)
        lib = gissmo_source;
        if ~all(isfield(lib, {'spectra_norm', 'ppm', 'names'}))
            error('gissmo_match_spectra:InvalidLibrary', ...
                'Library struct must be created by gissmo_load_library.');
        end
    elseif ischar(gissmo_source) || isstring(gissmo_source)
        lib = gissmo_load_library(string(gissmo_source));
    else
        error('gissmo_match_spectra:InvalidSource', ...
            'gissmo_source must be a file path or a library struct.');
    end

    ns  = size(spec_matrix, 1);
    n   = size(lib.spectra_norm, 1);
    pts = size(lib.spectra_norm, 2);

    max_allowed_shift = floor(pts / 2);
    if options.SearchMode == "split"
        max_allowed_shift = floor(floor(pts / 2) / 2);
    end
    if options.ShiftMax > max_allowed_shift
        error('gissmo_match_spectra:ShiftTooLarge', ...
            'ShiftMax (%d) exceeds maximum allowed shift (%d) for %s mode.', ...
            options.ShiftMax, max_allowed_shift, options.SearchMode);
    end

    % 2. Interpolate, threshold, and L2-normalize targets onto GISSMO ppm axis
    if options.Verbose
        disp('Aligning and normalizing consensus spectra onto GISSMO ppm axis...');
    end
    norm_targets = preprocess_targets(spec_matrix, ppm_target, lib.ppm, options.MinFrac);

    % 3. Scan GISSMO library using original shift grid j = 1 : ShiftStep : 2*ShiftMax
    if options.Verbose
        fprintf('Scanning %d GISSMO compounds for %d consensus spectra...\n', n, ns);
    end

    score_high      = zeros(ns, 1);
    lib_best        = zeros(ns, 1);
    offset_l        = zeros(ns, 1);
    offset_r        = zeros(ns, 1);
    shifted_spectra = zeros(pts, ns);
    matched_lib     = zeros(pts, ns);

    for i = 1:ns
        if options.SearchMode == "global"
            [score_high(i), lib_best(i), offset_l(i), offset_r(i)] = ...
                search_global(lib.spectra_norm, norm_targets(:, i), options.ShiftMax, options.ShiftStep);
        else
            [score_high(i), lib_best(i), offset_l(i), offset_r(i)] = ...
                search_split(lib.spectra_norm, norm_targets(:, i), options.ShiftMax, options.ShiftStep);
        end
        shifted_spectra(:, i) = apply_plot_offset(norm_targets(:, i), offset_l(i), offset_r(i));
        matched_lib(:, i)     = lib.spectra_norm(lib_best(i), :).';
    end

    % 4. Build results tables
    target_id   = (1:ns).';
    match_names = lib.names(lib_best);

    results = table(target_id, match_names, lib_best, score_high, offset_l, offset_r, ...
        'VariableNames', {'Target_ID', 'Match_Name', 'Match_Index', 'Score', 'Offset_L', 'Offset_R'});
    sorted_results = sortrows(results, 'Score', 'descend');

    if options.Verbose
        print_summary(results, sorted_results);
    end

    % 5. Package output
    out = struct();
    out.results             = results;
    out.sorted_results      = sorted_results;
    out.ppm_lib             = lib.ppm;
    out.shifted_spectra     = shifted_spectra;
    out.matched_lib_spectra = matched_lib;
    out.options             = options;

    % 6. Optional automated plotting
    if options.PlotMatches
        gissmo_plot_matches(out);
    end
end

%% --- LOCAL HELPER FUNCTIONS ---

function norm_targets = preprocess_targets(spec_matrix, ppm_target, ppm_lib, min_frac)
    ns  = size(spec_matrix, 1);
    pts = numel(ppm_lib);
    norm_targets = zeros(pts, ns);

    for i = 1:ns
        aligned = interp1(ppm_target, spec_matrix(i, :), ppm_lib, 'linear', 0.0);

        % Always apply peak threshold (when min_frac == 0, clips negative values to 0)
        max_pk = max(aligned);
        aligned(aligned < min_frac * max_pk) = 0.0;

        scl = norm(aligned);
        if scl > 0
            norm_targets(:, i) = (aligned ./ scl).';
        end
    end
end

function [best_score, best_idx, off_l, off_r] = search_global(lib_norm, spec, shft_max, shft_step)
    sz = numel(spec);

    if shft_max == 0
        scores = lib_norm * spec;
        [best_score, best_idx] = max(scores);
        off_l = 0;
        off_r = 0;
        return;
    end

    j_vals      = 1 : shft_step : 2 * shft_max;
    num_j       = numel(j_vals);
    test_shft   = zeros(sz, 1);
    shifted_mat = zeros(sz, num_j);
    src_slice   = spec(shft_max : sz - shft_max);

    % Preserve persistent testShft buffer across j iterations to match libSearchOff
    for idx = 1:num_j
        j = j_vals(idx);
        test_shft(j : sz - 2 * shft_max + j) = src_slice;
        shifted_mat(:, idx) = test_shft;
    end

    scores = lib_norm * shifted_mat; % [n_lib x num_j]
    [best_score, lin_idx] = max(scores(:));
    [best_idx, best_col]  = ind2sub(size(scores), lin_idx);

    off_l = j_vals(best_col) - shft_max;
    off_r = off_l;
end

function [best_score, best_idx, off_l, off_r] = search_split(lib_norm, spec, shft_max, shft_step)
    szh = floor(numel(spec) / 2);

    if shft_max == 0
        scores = lib_norm(:, 1 : 2 * szh) * spec(1 : 2 * szh);
        [best_score, best_idx] = max(scores);
        off_l = 0;
        off_r = 0;
        return;
    end

    j_vals = 1 : shft_step : 2 * shft_max;
    num_j  = numel(j_vals);

    % Left-half shifted states (persistent buffer across j iterations)
    buf_l     = zeros(szh, 1);
    shifted_l = zeros(szh, num_j);
    src_l     = spec(shft_max : szh - shft_max);
    for idx = 1:num_j
        j = j_vals(idx);
        buf_l(j : szh - 2 * shft_max + j) = src_l;
        shifted_l(:, idx) = buf_l;
    end

    % Right-half shifted states (persistent buffer across k and j iterations):
    % j_idx == 1 starts from zeros; j_idx >= 2 starts from end state of j_idx == 1.
    buf_r           = zeros(szh, 1);
    shifted_r_first = zeros(szh, num_j);
    shifted_r_rest  = zeros(szh, num_j);
    src_r           = spec(szh + shft_max : 2 * szh - shft_max);
    for idx = 1:num_j
        k = j_vals(idx);
        buf_r(k : szh - 2 * shft_max + k) = src_r;
        shifted_r_first(:, idx) = buf_r;
    end
    for idx = 1:num_j
        k = j_vals(idx);
        buf_r(k : szh - 2 * shft_max + k) = src_r;
        shifted_r_rest(:, idx) = buf_r;
    end

    scores_l       = lib_norm(:, 1:szh)         * shifted_l;       % [n_lib x num_j]
    scores_r_first = lib_norm(:, szh+1 : 2*szh) * shifted_r_first; % [n_lib x num_j]
    scores_r_rest  = lib_norm(:, szh+1 : 2*szh) * shifted_r_rest;  % [n_lib x num_j]

    % Match exact (j, k, i) iteration order and strict '>' tie-breaking of libSearchMid
    best_score = -1.0;
    best_idx   = 1;
    off_l      = 0;
    off_r      = 0;

    for j_idx = 1:num_j
        if j_idx == 1
            sr = scores_r_first;
        else
            sr = scores_r_rest;
        end
        sl_col = scores_l(:, j_idx);
        for k_idx = 1:num_j
            [score_max, score_idx] = max(sl_col + sr(:, k_idx));
            if score_max > best_score
                best_score = score_max;
                best_idx   = score_idx;
                off_l      = j_vals(j_idx) - shft_max;
                off_r      = j_vals(k_idx) - shft_max;
            end
        end
    end
end

function test_off = apply_plot_offset(spec, off_l, off_r)
    pts      = numel(spec);
    mid      = floor(pts / 2);
    test_off = spec;

    % Offset test spectrum leaving midpoint and boundary segments unchanged
    for ii = abs(off_l) + 1 : mid - abs(off_l)
        test_off(ii) = spec(ii - off_l);
    end
    for ii = mid + abs(off_r) : pts - abs(off_r)
        test_off(ii) = spec(ii - off_r);
    end
end

function print_summary(results, sorted_results)
    disp('-------------------------------------------------------------------------');
    fprintf('Target # | GISSMO Match Name (Index)              | Score | offL | offR \n');
    disp('-------------------------------------------------------------------------');
    for i = 1:height(results)
        fprintf('   %2d    | %-30s (%4d) | %1.3f |  %3d |  %3d \n', ...
            results.Target_ID(i), results.Match_Name(i), results.Match_Index(i), ...
            results.Score(i), results.Offset_L(i), results.Offset_R(i));
    end
    disp('-------------------------------------------------------------------------');

    fprintf('\n=========================================================================\n');
    disp(' TOP GISSMO MATCHES (RANKED BY SCORE) ');
    fprintf('=========================================================================\n');
    fprintf('Rank | Target | GISSMO Match Name              | Score | GISSMO Index \n');
    disp('-------------------------------------------------------------------------');
    for k = 1:height(sorted_results)
        fprintf(' %2d  |   %2d   | %-30s | %1.3f |     %4d \n', ...
            k, sorted_results.Target_ID(k), sorted_results.Match_Name(k), ...
            sorted_results.Score(k), sorted_results.Match_Index(k));
    end
    disp('-------------------------------------------------------------------------');
end
