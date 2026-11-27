function [consensus_spectra, consensus_elutions, best_member_selected, cluster_sizes] = optimizeClusterElutions(lib, spectra_matrix, all_elutions, elution_norms_sq, cluster_idx, max_iters, snap_spectra_to_winner)
% OPTIMIZECLUSTERELUTIONS Initialize consensus clusters and optimize member elutions.
%   Computes initial median consensus spectra and elutions for each cluster,
%   then runs iterative alternating member selection to choose the member
%   elution profile per multi-member cluster that minimizes global Frobenius
%   reconstruction error against lib. Optionally snaps consensus spectra to
%   the winning member of each cluster.
%
%   External Functions Required:
%       None
%
%   Toolboxes Required:
%       None
%
%   Inputs:
%       lib                    - MxP spectral library matrix.
%       spectra_matrix         - NxP matrix of extracted clean spectra.
%       all_elutions           - NxM matrix of pre-padded, singular-value-scaled
%                                node elution profiles.
%       elution_norms_sq       - Nx1 vector of squared L2 norms of all_elutions rows
%                                (retained in signature for interface compatibility; unused).
%       cluster_idx            - Nx1 cluster assignment vector (1..C).
%       max_iters              - Maximum number of alternating selection passes.
%       snap_spectra_to_winner - Logical flag to replace median consensus spectra
%                                with the winning member's spectrum after convergence.
%
%   Outputs:
%       consensus_spectra    - CxP matrix of cluster consensus spectra.
%       consensus_elutions   - CxM matrix of optimized cluster elution profiles.
%       best_member_selected - Cx1 vector of winning node indices per cluster.
%       cluster_sizes        - Cx1 vector of member counts per cluster.

    num_total_clusters = max(cluster_idx);
    num_spectral_points = size(spectra_matrix, 2);
    num_elution_points = size(all_elutions, 2);

    consensus_spectra    = zeros(num_total_clusters, num_spectral_points);
    consensus_elutions   = zeros(num_total_clusters, num_elution_points);
    cluster_sizes        = zeros(num_total_clusters, 1);
    best_member_selected = zeros(num_total_clusters, 1);
    cluster_members      = cell(num_total_clusters, 1);

    % 1. Initialize consensus spectra and elutions per cluster
    for c = 1:num_total_clusters
        members = find(cluster_idx == c);
        cluster_members{c} = members;
        c_size = length(members);
        cluster_sizes(c) = c_size;

        if c_size > 1
            consensus_spectra(c, :)  = median(spectra_matrix(members, :), 1);
            consensus_elutions(c, :) = median(all_elutions(members, :), 1);
            best_member_selected(c)  = 0;
        else
            consensus_spectra(c, :)  = spectra_matrix(members, :);
            consensus_elutions(c, :) = all_elutions(members, :);
            best_member_selected(c)  = members(1);
        end
    end

    % 2. Alternating member selection minimizing global Frobenius reconstruction error
    for iter = 1:max_iters
        members_changed = false;

        for c = 1:num_total_clusters
            members = cluster_members{c};
            if length(members) == 1
                continue;
            end

            best_err = inf;
            best_elution_for_c = zeros(1, num_elution_points);
            best_m = 0;

            for m_idx = 1:length(members)
                node_idx = members(m_idx);
                test_elution = all_elutions(node_idx, :);

                consensus_elutions(c, :) = test_elution;
                reconstructed_lib = consensus_elutions' * consensus_spectra;
                err = norm(lib - reconstructed_lib, 'fro');

                if err < best_err
                    best_err = err;
                    best_elution_for_c = test_elution;
                    best_m = node_idx;
                end
            end

            consensus_elutions(c, :) = best_elution_for_c;
            if best_member_selected(c) ~= best_m
                best_member_selected(c) = best_m;
                members_changed = true;
            end
        end

        if ~members_changed
            break;
        end
    end


    % 3. Optionally snap consensus spectra to winning node
    if snap_spectra_to_winner
        for c = 1:num_total_clusters
            winning_node = best_member_selected(c);
            if winning_node > 0
                consensus_spectra(c, :) = spectra_matrix(winning_node, :);
            end
        end
    end
end
