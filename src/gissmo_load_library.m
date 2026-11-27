function lib = gissmo_load_library(mat_path)
% GISSMO_LOAD_LIBRARY Load, validate, and L2-normalize a GISSMO reference library.
%
%   lib = gissmo_load_library(mat_path)
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       None
%
%   Input:
%       mat_path - Path to GISSMO .mat file containing lib_combined, ppm_lib,
%                  and combined_names.
%
%   Output:
%       lib      - Struct with fields:
%                  .spectra_norm  [n_lib x n_pts] L2-normalized reference spectra
%                  .ppm           [1 x n_pts]     chemical shift axis
%                  .names         [n_lib x 1]     string array of compound names
%                  .source_path   [1 x 1]         string path to loaded file

    arguments
        mat_path (1,1) string
    end

    if ~isfile(mat_path)
        error('gissmo_load_library:FileNotFound', ...
            'GISSMO library file not found: %s', mat_path);
    end

    raw = load(mat_path);
    required_fields = {'lib_combined', 'ppm_lib', 'combined_names'};
    for k = 1:numel(required_fields)
        if ~isfield(raw, required_fields{k})
            error('gissmo_load_library:MissingField', ...
                'Expected variable "%s" not found in %s.', ...
                required_fields{k}, mat_path);
        end
    end

    lib_spectra = double(raw.lib_combined);
    ppm_lib     = double(raw.ppm_lib(:)).';
    names       = string(raw.combined_names(:));

    if size(lib_spectra, 2) ~= numel(ppm_lib)
        error('gissmo_load_library:DimensionMismatch', ...
            'lib_combined has %d columns, but ppm_lib has %d elements.', ...
            size(lib_spectra, 2), numel(ppm_lib));
    end

    if size(lib_spectra, 1) ~= numel(names)
        error('gissmo_load_library:DimensionMismatch', ...
            'lib_combined has %d rows, but combined_names has %d entries.', ...
            size(lib_spectra, 1), numel(names));
    end

    % L2 normalize each reference spectrum
    norms = vecnorm(lib_spectra, 2, 2);
    norms(norms == 0) = 1;

    lib = struct();
    lib.spectra_norm = lib_spectra ./ norms;
    lib.ppm          = ppm_lib;
    lib.names        = names;
    lib.source_path  = mat_path;
end
