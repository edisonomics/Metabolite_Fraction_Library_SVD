function figs = plotSteppedReconstructions(lib, ppm, best_state, k_steps, include_full, include_raw)
% PLOTSTEPPEDRECONSTRUCTIONS Display 3D library reconstructions across component counts.
%   Sorts winning consensus components by descending maximum elution intensity
%   and calls displaySpectra3D for each component count in k_steps, followed
%   optionally by the full optimal_k reconstruction and the raw library.
%
%   External Functions Required:
%       displaySpectra3D
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       lib          - MxP original spectral library matrix (rows = fractions, cols = ppm).
%       ppm          - 1xP chemical shift vector.
%       best_state   - Struct returned by optimizeClusterCutoff containing
%                      'optimal_k', 'consensus_spectra', and 'consensus_elutions'.
%       k_steps      - (Optional) Row vector of component counts to reconstruct
%                      (default: [5, 10, 25]).
%       include_full - (Optional) Logical flag to plot the full optimal_k
%                      intensity-sorted reconstruction (default: true).
%       include_raw  - (Optional) Logical flag to plot the raw library (default: true).
%
%   Outputs:
%       figs - Struct containing figure handles:
%              'steps'            (1xK array of handles for plotted k_steps),
%              'step_ks'          (1xK vector of k values actually plotted),
%              'full_reconstruct' (figure handle or []),
%              'raw'              (figure handle or []).

    if nargin < 4 || isempty(k_steps)
        k_steps = [5, 10, 25];
    end
    if nargin < 5 || isempty(include_full)
        include_full = true;
    end
    if nargin < 6 || isempty(include_raw)
        include_raw = true;
    end

    figs = struct('steps', gobjects(0), 'step_ks', [], 'full_reconstruct', [], 'raw', []);

    if isempty(fieldnames(best_state)) || isempty(best_state.consensus_spectra)
        warning('plotSteppedReconstructions:EmptyState', ...
            'best_state is empty; skipping stepped 3D reconstructions.');
        return;
    end

    % Sort components by highest elution intensity (descending)
    opt_k    = best_state.optimal_k;
    spectra  = best_state.consensus_spectra;
    elutions = best_state.consensus_elutions;

    elution_max_ints = max(elutions, [], 2);
    [~, sort_idx] = sort(elution_max_ints, 'descend');

    spectra  = spectra(sort_idx, :);
    elutions = elutions(sort_idx, :);

    num_avail = size(spectra, 1);
    max_k     = min(num_avail, opt_k);

    % Plot requested stepped reconstructions
    for i = 1:length(k_steps)
        k = k_steps(i);
        if k > num_avail
            warning('plotSteppedReconstructions:ExceedsComponents', ...
                'Skipping k = %d because only %d components are available.', k, num_avail);
            continue;
        end
        displaySpectra3D(elutions(1:k, :)' * spectra(1:k, :), ppm);
        title(num2str(k));
        figs.steps(end + 1)   = gcf;
        figs.step_ks(end + 1) = k;
    end

    % Plot full optimal_k intensity-sorted reconstruction
    if include_full
        displaySpectra3D(elutions(1:max_k, :)' * spectra(1:max_k, :), ppm);
        title('Full Reconstruct');
        figs.full_reconstruct = gcf;
    end

    % Plot raw library
    if include_raw
        displaySpectra3D(lib, ppm);
        title('Raw');
        figs.raw = gcf;
    end
end
