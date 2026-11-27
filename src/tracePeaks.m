function nodes = tracePeaks(peaks_struct, ppm_tol, min_length)
% TRACEPEAKS Links detected peaks across consecutive fractions into traces.
%
% Connects peaks across adjacent fractions within ppm_tol, pads the start
% and end of each trace with zero-intensity base points where bounds allow,
% and removes any resulting traces shorter than min_length.
%
% External Functions Required:
%   None
%
% Toolboxes Required:
%   None
%
% Inputs:
%   peaks_struct - Struct array from findLibraryPeaks (fields: peak_locs, peaks, ppm)
%   ppm_tol      - (Optional) Maximum ppm shift between adjacent fractions (default: 0.003)
%   min_length   - (Optional) Minimum number of points to keep a trace (default: 3)
%
% Outputs:
%   nodes        - Struct array of peak traces with fields: fracs, ppms, locs, peaks

    if nargin < 2
        ppm_tol = 0.003;
    end
    if nargin < 3
        min_length = 3;
    end

    num_fracs = size(peaks_struct, 2);
    nodes = struct('fracs', {}, 'ppms', {}, 'locs', {}, 'peaks', {});
    peaks_work = peaks_struct;
    big_counter = 1;

    % Trace peaks forward across consecutive fractions
    for i = 1:num_fracs - 1
        for j = 1:length(peaks_work(i).ppm)
            locs  = peaks_work(i).peak_locs(j);
            ppms  = peaks_work(i).ppm(j);
            peaks = peaks_work(i).peaks(j);
            fracs = i;
            counter = 1;

            for k = i + 1:num_fracs
                if min(abs(peaks_work(k).ppm - ppms(end))) <= ppm_tol
                    [~, index] = min(abs(peaks_work(k).ppm - ppms(counter)));
                    fracs(counter + 1) = k;
                    ppms(counter + 1)  = peaks_work(k).ppm(index);
                    locs(counter + 1)  = peaks_work(k).peak_locs(index);
                    peaks(counter + 1) = peaks_work(k).peaks(index);
                    counter = counter + 1;

                    % Remove matched peak so it is not reused by later traces
                    peaks_work(k).ppm(index)       = [];
                    peaks_work(k).peak_locs(index) = [];
                    peaks_work(k).peaks(index)     = [];
                else
                    break
                end
            end

            nodes(big_counter).fracs = fracs;
            nodes(big_counter).ppms  = ppms;
            nodes(big_counter).locs  = locs;
            nodes(big_counter).peaks = peaks;
            big_counter = big_counter + 1;
        end
    end

    % Pad the beginning and end of each trace with zero-intensity anchor points
    for i = 1:length(nodes)
        if nodes(i).fracs(1) > 1
            nodes(i).fracs = [nodes(i).fracs(1) - 1, nodes(i).fracs];
            nodes(i).locs  = [nodes(i).locs(1), nodes(i).locs];
            nodes(i).ppms  = [nodes(i).ppms(1), nodes(i).ppms];
            nodes(i).peaks = [0, nodes(i).peaks];
        end

        if nodes(i).fracs(end) < num_fracs
            nodes(i).fracs = [nodes(i).fracs, nodes(i).fracs(end) + 1];
            nodes(i).locs  = [nodes(i).locs, nodes(i).locs(end)];
            nodes(i).ppms  = [nodes(i).ppms, nodes(i).ppms(end)];
            nodes(i).peaks = [nodes(i).peaks, 0];
        end
    end

    % Delete traces shorter than min_length
    indexer = false(1, length(nodes));
    for i = 1:length(nodes)
        if length(nodes(i).peaks) < min_length
            indexer(i) = true;
        end
    end
    nodes(indexer) = [];
end
