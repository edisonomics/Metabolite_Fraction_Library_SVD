%% Perform SVD Analysis on mFL
%% Activate Functions
activate_metabolomics_toolbox
%% Read in data
%{
    The next two section have two different ways to read in the NMR data.
    1. Use the first section if the data is in NMRPipe Format
    2. Use the second section if the data is in Bruker Processed Format
%}
%% READ IN NMRPipe DATA
% path should be to a folder that containes only ft1 files
folder = pwd;
cd 'path_to_ft1_files'
spectra = load1DPipe();
cd(folder)
[X,ppm_ft]=Setup1D(spectra);
clear spectra folder
%% READ IN Bruker DATA
% path should be to folder that contains numbered TopSpin spectral folders
path = 'path_to_bruker_folders';
spectra = load1DTopspin(path,'bruker');
[X,ppm_ft]=Setup1D(spectra);
clear spectra path
%% Optional Downsampling
%{
    Use this section to downsample the data. This will speed up analysis
    time
%}
num_of_points_smaller = 5000;
ppm_small = linspace(ppm_ft(1),ppm_ft(end),num_of_points_smaller);
X_small = zeros(size(X,1),num_of_points_smaller);
for i = 1:size(X_small,1)
    X_small(i,:) = interp1(ppm_ft,X(i,:),ppm_small,'linear','extrap');
end
ppm_ft = ppm_small;
X = X_small;
clear num_of_points_smaller X_small ppm_small
%% Visualize the data
displaySpectra3D(X,ppm_ft)
%% Preprocess the fractions
%{
    This section will remove solvent or other peaks from the library and
    replace them with noise. Peaks to remove would be any peaks that span
    all the fractions. Ex/ DSS, Formate, Water.

    

%}
%%%
%%%
%%%
% Find region where this nothing but noise across all of the fractions
ppm_noise_start = 0.3;
ppm_noise_end = 0.5;

%{

More regions can be added
ex/ 
region4_mask = (ppm_ft >= # & ppm_ft <= #);
target_noise_mask = region1_mask | region2_mask | region3_mask | region4_mask;
         
%}
% Define the target ppm regions to be turned to noise
region1_mask = (ppm_ft >= 8.305 & ppm_ft <= 8.589);%formate
region2_mask = (ppm_ft >= 4.767 & ppm_ft <= 4.843);%water
region3_mask = (ppm_ft >= -0.0835 & ppm_ft <= 0.0808);%DSS

% Combine the regions into a single logical mask
target_noise_mask = region1_mask | region2_mask | region3_mask;
%%%
%%%
%%%

X_noisy = X;
% Isolate the columns corresponding to the 0.2 to 0.6 ppm noise region in X
noise_region_mask = (ppm_ft >= ppm_noise_start & ppm_ft <= ppm_noise_end);
X_noise_subset = X(:, noise_region_mask);

% Estimate the baseline noise mean and standard deviation from this specific region
noise_mean = mean(X_noise_subset(:));
noise_std = std(X_noise_subset(:)); 

% Get the total number of rows in the expanded noisy matrix and the number of target columns
[num_total_rows, ~] = size(X);
num_target_cols = sum(target_noise_mask);

% Overwrite the specified regions across all rows with matching baseline noise
X_noisy(:, target_noise_mask) = noise_mean + noise_std * randn(num_total_rows, num_target_cols);

X_noisy(X_noisy<= 0) = 0.0;
X_noisy=X_noisy/max(X_noisy,[],"all");
% Plot Fractions and Draw Bounding Boxes
displaySpectra3D(X_noisy, ppm_ft);
clear *noise* num* region*
%% Make final variable
% All other variables will be deleted to save memory
ppm = ppm_ft;
lib = X_noisy;
clearvars -except lib ppm
%% Find noise threshold
figure,plot(ppm,lib)
set(gca,'XDir','rev')
%% Visualize all of the peaks
%{
    This section will show all the peaks in the library. Best reults will
    have a noise threshold that is not extremely low. Choose a threshold
    that does not result in major streaks of pink circles that span the
    whole library near the noise. For example, a threshold that is 2.5x the
    noise level may be a good starting point.

%}
noise_threshold = 0.05;
peaks_struct = struct;
for i = 1:size(lib,1)
    [~,peaks_struct(i).peak_locs] = findpeaks(lib(i,:),"MinPeakHeight",noise_threshold);
    [peaks_struct(i).peaks,peaks_struct(i).ppm] = findpeaks(lib(i,:),ppm,"MinPeakHeight",noise_threshold);
end
% 3D great figure
clear x y z
for i = 1:size(lib,1)
    x(i,:) = ppm;
end
y = lib;
for i = 1:size(lib,1)
    z(i,:) = (zeros(1,size(lib,2)) + i);
end
colors = winter(size(lib,1));
figure
hold on
for i = 1:size(lib,1)
    plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
    plot3(x(i,peaks_struct(i).peak_locs),z(i,peaks_struct(i).peak_locs),y(i,peaks_struct(i).peak_locs),"o","Color","m")
end
hold off
set(gca,'XDir','rev');
view(0,70)
ax = gca;
ax.Clipping = "off";
clear colors ax i x y z
%% Now do the tracing
nodes = struct;
peaks_struct2 = peaks_struct;
big_counter = 1;
for i = 1:size(peaks_struct2,2)-1
    for j = 1:length(peaks_struct2(i).ppm)
        ppms = [];
        fracs = [];
        locs = [];
        peaks = [];
        locs(1) = peaks_struct2(i).peak_locs(j);
        ppms(1) = peaks_struct2(i).ppm(j);
        peaks(1) = peaks_struct2(i).peaks(j);
        fracs(1) = i;
        counter = 1;
        for k = i+1:size(peaks_struct2,2)
            if min(abs(peaks_struct2(k).ppm-ppms(end))) <= 0.003
                [~,index] = min(abs(peaks_struct2(k).ppm-ppms(counter)));
                fracs(counter+1) = k;
                ppms(counter+1) = peaks_struct2(k).ppm(index);
                locs(counter+1) = peaks_struct2(k).peak_locs(index);
                peaks(counter+1) = peaks_struct2(k).peaks(index);
                counter = counter + 1;
                peaks_struct2(k).ppm(index) = [];
                peaks_struct2(k).peak_locs(index) = [];
                peaks_struct2(k).peaks(index) = [];
            else
                break
            end
        end
        nodes(big_counter).fracs = fracs;
        nodes(big_counter).ppms = ppms;
        nodes(big_counter).locs = locs;
        nodes(big_counter).peaks = peaks;
        big_counter = big_counter + 1;
    end
end

% --- NEW CODE: Add extra points at the beginning and end of traces ---
for i = 1:length(nodes)
    % Pad the beginning if it doesn't exceed the first fraction
    if nodes(i).fracs(1) > 1
        nodes(i).fracs = [nodes(i).fracs(1) - 1, nodes(i).fracs];
        nodes(i).locs  = [nodes(i).locs(1), nodes(i).locs];
        nodes(i).ppms  = [nodes(i).ppms(1), nodes(i).ppms];
        nodes(i).peaks = [0, nodes(i).peaks]; % Base level set to 0
    end
    
    % Pad the end if it doesn't exceed the last fraction
    if nodes(i).fracs(end) < size(lib,1)
        nodes(i).fracs = [nodes(i).fracs, nodes(i).fracs(end) + 1];
        nodes(i).locs  = [nodes(i).locs, nodes(i).locs(end)];
        nodes(i).ppms  = [nodes(i).ppms, nodes(i).ppms(end)];
        nodes(i).peaks = [nodes(i).peaks, 0]; % Base level set to 0
    end
end
% ---------------------------------------------------------------------
% Delete any nodes that are less than 3
indexer = zeros(1,length(nodes));
for i = 1:length(nodes)
    if length(nodes(i).peaks) < 3
        indexer(i) = 1;
    end
end
nodes(logical(indexer)) = [];

clear big_counter counter fracs ppms i j k index locs peaks_struct2 peaks_struct
clear peaks indexer
% Plot how the tracing did
clear x y z
for i = 1:size(lib,1)
    x(i,:) = ppm;
end
y = lib;
for i = 1:size(lib,1)
    z(i,:) = (zeros(1,length(ppm)) + i);
end
colors = winter(size(lib,1));
fig = figure('Position',[10 10 900 750]);
ax1 = axes(fig);
hold on
for i = 1:size(lib,1)
    plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
end
view([0 70])
for i = 1:size(nodes,2)
    plot3(ppm(nodes(i).locs),nodes(i).fracs,nodes(i).peaks,'LineWidth',1.5)
end
hold off
set(gca,'XDir','rev');
ax1.Clipping = "off";
clear x y z colors i fig ax1 originalZLim
%% Remove Redundant Traces
% Group traces by their exact fractional spans, find those >= 0.99 correlated, 
% and keep the highest intensity

% 1. Create a string representation of the exact fractions each trace spans
frac_strs = cell(1, length(nodes));
for i = 1:length(nodes)
    frac_strs{i} = mat2str(nodes(i).fracs);
end

% 2. Find the unique fractional spans
[unique_spans, ~, ic] = unique(frac_strs);
indices_to_keep = [];

for i = 1:length(unique_spans)
    % Find all traces that share this exact sequence of fractions
    idx = find(ic == i);
    
    if isscalar(idx)
        indices_to_keep = [indices_to_keep, idx];
        continue;
    end
    
    % Extract peaks to a matrix for correlation
    L = length(nodes(idx(1)).fracs);
    P = zeros(length(idx), L);
    for j = 1:length(idx)
        P(j, :) = nodes(idx(j)).peaks;
    end
    
    % Determine highest absolute intensity for each trace
    max_ints = max(abs(P), [], 2);
    
    % Sort traces by maximum absolute intensity in descending order
    [~, sort_order] = sort(max_ints, 'descend');
    sorted_idx = idx(sort_order);
    P_sorted = P(sort_order, :);
    
    % Boolean mask to track which traces we are keeping
    keep_mask = true(length(sorted_idx), 1);
    
    for m = 1:length(sorted_idx)
        if keep_mask(m) && m < length(sorted_idx)
            % Compute Pearson correlation of the current highest intensity trace
            % with all remaining traces of lower intensity
            v1 = P_sorted(m, :)';
            V_others = P_sorted(m+1:end, :)';
            
            C = corr(v1, V_others);
            C(isnan(C)) = 0; % Handle any potential flat/zero-variance arrays safely
            
            % Find traces that are >= 0.99 correlated with the current one
            corr_idx = find(C >= 0.99) + m;
            
            % Mark them for removal
            keep_mask(corr_idx) = false;
        end
    end
    
    % Append the indices we are keeping from this exact span group
    indices_to_keep = [indices_to_keep, reshape(sorted_idx(keep_mask), 1, [])];
end

% Create the new variable with the filtered traces
nodes_filtered = nodes(indices_to_keep);

clear frac_strs unique_spans ic indices_to_keep i idx L P j max_ints
clear sort_order sorted_idx P_sorted keep_mask m v1 V_others C corr_idx
% Plot how the tracing did after trace removal
clear x y z
for i = 1:size(lib,1)
    x(i,:) = ppm;
end
y = lib;
for i = 1:size(lib,1)
    z(i,:) = (zeros(1,length(ppm)) + i);
end
colors = winter(size(lib,1));
fig = figure('Position',[10 10 900 750]);
ax1 = axes(fig);
hold on

% Plot original spectra
for i = 1:size(lib,1)
    plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
end
view([0 70])

% Plot filtered traces
for i = 1:size(nodes_filtered,2)
    plot3(ppm(nodes_filtered(i).locs), nodes_filtered(i).fracs, nodes_filtered(i).peaks, 'LineWidth', 1.5)
end

hold off
set(gca,'XDir','rev');
ax1.Clipping = "off";
clear x y z colors i fig ax1
%% Reorder by highest intensity nodes
% Reorder nodes_filtered to 
index_max = zeros(length(nodes_filtered),1);
for i = 1:length(index_max)
    index_max(i) = max(nodes_filtered(i).peaks);
end
[~,idx] = sort(index_max,'descend');
nodes_filtered_reordered = nodes_filtered(idx);
clear i index_max idx
%% Do the SVD Analysis
%{
    The sections below will do the main SVD analysis. The section directly
    below will take aroun 10 - 40 minutes
%}
%% Begin with the loop
AmaxScl = 0.7;
for i = 1:length(nodes_filtered_reordered)
%for i = 1:1
    matrix_loop = lib(nodes_filtered_reordered(i).fracs(1):nodes_filtered_reordered(i).fracs(end),:);
    [U,S,V] = svd(matrix_loop,'econ');
    Vsub = V(:,1:size(matrix_loop,1))';
    raw_peak = nodes_filtered_reordered(i).peaks(:); % Force as column vector
    Scl = sum(raw_peak.^2);
    peakIn = raw_peak ./ sqrt(Scl); % L2 Normalization
    [~, row] = max(peakIn);
    Amax = peakIn(row)*AmaxScl;
    tic
    wt = optimizeElution1(U,peakIn,row,Amax);
    toc
    % --- NEW: Single-Target Covariance ---
    n_points = size(matrix_loop, 1);

    % Calculate covariance and correlation directly against matrix_loop
    covar = (peakIn - mean(peakIn))' * (matrix_loop - repmat(mean(matrix_loop), n_points, 1)) ./ (n_points - 1);
    corr_spectrum = (zscore(peakIn)' * zscore(matrix_loop)) ./ (n_points - 1);

    % Apply your exact same correlation gates and intensity floors
    covar(corr_spectrum < 0.97) = 0;
    covar(isnan(covar)) = 0;
    covar_max = max(covar);
    if covar_max > 0
        covar(covar < (0.1 * covar_max)) = 0;
    end

    % L2 Normalize to create specIn
    Scl_spec = sum(covar.^2);
    if Scl_spec > 0
        specIn = covar ./ sqrt(Scl_spec);
    else
        specIn = covar;
    end
    tic
    SclS = 0.1; % scale factor for difference between wt and wt2

    max_iter = 10;   % Failsafe to prevent infinite loops
    tol = 1e-2;      % Convergence threshold (1% change)
    iter = 0;
    converged = false;

    % Initialize with the current elution weights (wt) before entering the loop
    wt_old = wt;

    disp('Starting ALS refinement loop...');

    while iter < max_iter && ~converged
        iter = iter + 1;

        % Step 1: Optimize Spectra (Frequency Domain) -> Updates wt2
        % Now uses the normalized covariance targets (specIn)
        wt2 = optimizeSpectra(Vsub, specIn, wt, SclS);

        % Step 2: Optimize Elution (Time Domain) -> Updates wt
        % Uses the newly calculated spectral weights (wt2) from Step 1
        wt = optimizeElution(U, peakIn, row, Amax, wt2, SclS);

        % --- CONVERGENCE METRIC ---
        % How much did the Elution weights (wt) change from the last loop?
        rel_change = norm(wt - wt_old, 'fro') / norm(wt_old, 'fro');

        fprintf('Iteration %d: Relative Change = %.6f\n', iter, rel_change);

        if rel_change < tol
            converged = true;
            fprintf('ALS converged successfully after %d iterations!\n', iter);
        end

        % Save the current wt to become the "old" wt for the next loop
        wt_old = wt;
    end

    if ~converged
        fprintf('Warning: ALS reached maximum iterations (%d) without fully converging.\n', max_iter);
    end
    toc
    % ==========================================================
    % --- FINAL EXTRACTION & TRANSFORMATION ---
    % ==========================================================

    % 1. Frequency Domain Extraction
    % Since wt2 is on the left side of the unknown space, use the Left Divide (\)
    Vtrn2 = wt2 \ Vsub;

    % 2. Time Domain Extraction (Using the full U matrix)
    % This is a direct projection, so it remains standard multiplication
    Utrn2 = U * wt;
    Utrn2(Utrn2 < 0) = 0;
    % 6. Re-normalize to create Utrn2
    norm_val = norm(Utrn2);
    if norm_val > 0
        Utrn2 = Utrn2 / norm_val;
    end

    % 3. Core Matrix (Optional - using full S matrix)
    Strn = abs(wt \ (S * wt2));

    % ==========================================================
    % --- DYNAMIC CORRELATION-MODULATED INTENSITY FILTER ---
    % ==========================================================

    v_raw = max(0, Vtrn2);
    max_peak = max(v_raw);

    if max_peak > 0
        corr_floor = 0.7; % Maximum intensity penalty applied here
        corr_ceil  = 0.9; % Minimum intensity penalty applied here

        base_cutoff  = 0.05; % 5% cutoff for highly correlated peaks
        bleed_cutoff = 0.15; % 15% cutoff for poorly correlated bleed

        % 1. Scale correlation (Reuse corr_spectrum calculated earlier!)
        corr_scaled = (corr_spectrum - corr_floor) ./ (corr_ceil - corr_floor);
        corr_scaled = max(0, min(1, corr_scaled)); % Clamp between 0 and 1

        % 2. Calculate the dynamic percentage threshold
        dynamic_pct = bleed_cutoff - corr_scaled .* (bleed_cutoff - base_cutoff);

        % 3. Create the 1xN absolute intensity threshold vector
        threshold_vector = dynamic_pct * max_peak;

        % 4. Apply the filter
        v_hard = v_raw;
        v_hard(v_raw < threshold_vector) = 0;   % Dynamic intensity zeroing ONLY

        % 5. Smooth the remaining cliff edges (Morphological Buffer)
        is_signal = v_hard > 0;
        cliff_mask = imdilate(is_signal, ones(1, 5)) & ~is_signal;

        v_smoothed = smoothdata(v_hard, 'gaussian', 7);
        v_hard(cliff_mask) = v_smoothed(cliff_mask);

        % 6. Re-normalize to create Vtrn4
        norm_val = norm(v_hard);
        if norm_val > 0
            Vtrn4 = v_hard / norm_val;
        else
            Vtrn4 = v_hard;
        end
    else
        Vtrn4 = v_raw; % Failsafe if the vector was entirely negative
    end

    % Save the final cleaned spectrum and optimized elution profile
    nodes_filtered_reordered(i).clean_spectrum = Vtrn4;
    nodes_filtered_reordered(i).optimized_elution = Utrn2;
    nodes_filtered_reordered(i).singular_value = Strn;
end
clear matrix_loop U S V Vsub raw_peak Scl peakIn row Amax 
clear wt n_points covar corr_spectrum covar_max Scl_spec 
clear specIn SclS max_iter tol iter converged wt_old wt2 
clear rel_change Vtrn2 Utrn2 Strn v_raw max_peak corr_floor 
clear corr_ceil base_cutoff bleed_cutoff corr_scaled 
clear dynamic_pct threshold_vector v_hard is_signal cliff_mask 
clear v_smoothed norm_val Vtrn4 i 

%%%
%%%
%%%
%%%
%%%
%%%
%%%
%%%
%%%
%% Quality Control
%% This section removes spectra that did not match well to target in SVD

% Preallocate arrays
num_nodes = length(nodes_filtered_reordered);
spec_correlations = zeros(num_nodes, 1);

for i = 1:num_nodes
    % 1. Recalculate the Target Covariance Spectrum (specIn)
    % (Requires 'lib' to be present in the workspace)
    matrix_loop = lib(nodes_filtered_reordered(i).fracs(1):nodes_filtered_reordered(i).fracs(end), :);
    n_points = size(matrix_loop, 1);

    raw_peak = nodes_filtered_reordered(i).peaks(:);
    Scl = sum(raw_peak.^2);
    peakIn = raw_peak ./ sqrt(Scl);

    % Calculate raw covariance and correlation
    covar = (peakIn - mean(peakIn))' * (matrix_loop - repmat(mean(matrix_loop), n_points, 1)) ./ (n_points - 1);
    corr_spectrum = (zscore(peakIn)' * zscore(matrix_loop)) ./ (n_points - 1);

    % Apply correlation gates and intensity floors identically to previous loop
    covar(corr_spectrum < 0.97) = 0;
    covar(isnan(covar)) = 0;
    covar_max = max(covar);
    if covar_max > 0
        covar(covar < (0.1 * covar_max)) = 0;
    end

    % Normalize to get specIn
    Scl_spec = sum(covar.^2);
    if Scl_spec > 0
        specIn = covar ./ sqrt(Scl_spec);
    else
        specIn = covar;
    end

    % 2. Extract the cleaned spectrum
    clean_spec = nodes_filtered_reordered(i).clean_spectrum(:);
    specIn = specIn(:);

    % 3. Calculate Correlation between Target and Cleaned Spectrum
    R = corrcoef(specIn, clean_spec);
    
    if size(R, 1) > 1 && ~isnan(R(1,2))
        spec_correlations(i) = R(1,2);
    else
        spec_correlations(i) = NaN; % Handle flatline/zero-variance edge cases
    end
end

% 4. Visualization
figure('Name', 'Spectral Correlation Analysis', 'Position', [100, 100, 800, 600]);

% Panel 1: Scatter plot of correlation by node
subplot(2, 1, 1);
scatter(1:num_nodes, spec_correlations, 35, 'filled', 'MarkerFaceColor', '#7E2F8E');
yline(0.8, 'r--', 'Target Threshold (0.8)', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);
title('Similarity: Clean Spectrum vs. Target Covariance Spectrum');
xlabel('Node Index');
ylabel('Correlation Coefficient (r)');
ylim([min([0, min(spec_correlations) - 0.1]) 1.1]);
grid on;

% Panel 2: Histogram of correlation distribution
subplot(2, 1, 2);
histogram(spec_correlations, 30, 'FaceColor', '#EDB120', 'EdgeColor', 'w');
title('Distribution of Spectral Correlations');
xlabel('Correlation Coefficient (r)');
ylabel('Frequency (Number of Nodes)');
grid on;

% Clear all temporary variables used in this visualization block
clear poor_corr_idx num_poor cols rows t k idx matrix_loop n_points ...
      raw_peak Scl peakIn covar corr_spectrum covar_max Scl_spec specIn clean_spec;

% 5. Clear all temporary variables created in this section
clear num_nodes  i matrix_loop n_points raw_peak Scl peakIn ...
      covar corr_spectrum covar_max Scl_spec specIn clean_spec R;

%% Remove any peaks that have a correlation below threshold or NaN
nodes_filtered_reordered(spec_correlations < 0.2 | isnan(spec_correlations)) = [];
spec_correlations_filt = spec_correlations(spec_correlations > 0.2 & ~isnan(spec_correlations));
%% This section combines redundant spetra

disp('Extracting Vtrn4 spectra for hierarchical clustering...');

% 1. Extract all purified spectra from the struct into a single 2D matrix
num_nodes = length(nodes_filtered_reordered);
num_spectral_points = length(nodes_filtered_reordered(1).clean_spectrum);
Vtrn4_matrix = zeros(num_nodes, num_spectral_points);

for i = 1:num_nodes
    Vtrn4_matrix(i, :) = nodes_filtered_reordered(i).clean_spectrum;
end

% 2. Calculate the Pairwise Correlation Matrix
% corr() calculates the standard Pearson correlation between the spectra.
% We transpose Vtrn4_matrix so columns represent individual target spectra.
corr_matrix = corr(Vtrn4_matrix');

% Clean up any NaNs (e.g., if a spectrum was perfectly flat/zeroed out)
corr_matrix(isnan(corr_matrix)) = 0;

% 3. Calculate Hierarchical Linkage 
% pdist calculates the distance (1 - correlation), which linkage requires.
% 'weighted' distance matches the methodology used in your alignment functions.
dist_metric = pdist(Vtrn4_matrix, 'correlation');
clustering_tree = linkage(dist_metric, 'complete');

% ==========================================================
% --- 2D VISUALIZATION: DENDROGRAM & SORTED HEATMAP ---
% ==========================================================
disp('Generating 2D clustering visualization...');

figCluster = figure('Name', 'Vtrn4 Hierarchical Correlation Clustering', ...
                    'Position', [150, 150, 1100, 800], 'Color', 'w');

% --- Left Panel: The Dendrogram ---
% We calculate the dendrogram first to extract the optimal leaf order (outperm),
% which we will use to sort the 2D matrix so highly correlated spectra group together.
ax1 = subplot(1, 5, 1);
[H, T, cluster_order] = dendrogram(clustering_tree, 0, 'Orientation', 'left');
set(H, 'LineWidth', 1.5, 'Color', 'k');
set(ax1, 'YDir', 'reverse', 'YTick', []); % Align with the matrix rows
title('Dendrogram', 'FontSize', 12, 'FontWeight', 'bold');
box off; axis off;

% --- Right Panel: The Sorted 2D Correlation Matrix ---
% Reorder the raw correlation matrix using the dendrogram's cluster order
sorted_corr_matrix = corr_matrix(cluster_order, cluster_order);

ax2 = subplot(1, 5, 2:5);
imagesc(sorted_corr_matrix);

% Apply colormap and correlation limits
colormap(ax2, jet(256));
cbar = colorbar;
cbar.Label.String = 'Pearson Correlation (r)';
cbar.Label.FontWeight = 'bold';
cbar.Label.FontSize = 12;
caxis([0 1]); % Assuming we are most interested in positive correlations

% Format the matrix axes to display the original node indices
title('2D Correlation Matrix (Hierarchically Clustered)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Clustered Node Index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Clustered Node Index', 'FontSize', 12, 'FontWeight', 'bold');

% Apply the sorted labels to the X and Y axes
set(gca, 'XTick', 1:num_nodes, 'XTickLabel', cluster_order, 'XTickLabelRotation', 90);
set(gca, 'YTick', 1:num_nodes, 'YTickLabel', cluster_order);

% Force square pixels for a clean grid look
axis square; 
grid on;

disp('Clustering complete.');
clear i num_nodes num_spectral_points dist_metric ...
      figCluster ax1 ax2 H T cluster_order sorted_corr_matrix cbar;

%% This section makes the final variable of the analysis
%{
    The final variable (best_state) will contain the basis set in
    best_state.consensus_spectra.

%}
% ==========================================================
% --- DEEP SWEEP: CORRELATION CUTOFF OPTIMIZATION ---
% ==========================================================
disp('Starting Deep Sweep for Optimal Correlation Cutoff...');

% --- SWEEP & PURITY CONTROLS ---
cutoff_range     = 0.5:0.02:0.95; % Range of correlation cutoffs to test
weight_covar     = 0.8;            % Purity: Weight for target covariance
weight_sim       = 0.0;            % Purity: Weight for intra-cluster similarity
weight_intensity = 0.2;            % Purity: Weight for initial library peak intensity
max_iters        = 20;             % Max iterations for elution optimization

% Setup loop trackers
num_sweeps         = length(cutoff_range);
sweep_min_errors   = zeros(num_sweeps, 1);
sweep_optimal_ks   = zeros(num_sweeps, 1);
sweep_num_clusters = zeros(num_sweeps, 1);

lib_norm = norm(lib, 'fro');

% Track the absolute best state to avoid rerunning at the end
global_best_error = inf;
best_state = struct();

% Dimension variables
num_spectral_points = size(Vtrn4_matrix, 2);
if isscalar(lib)
    num_elution_points = lib;
else
    num_elution_points = size(lib, 1); 
end

% === MASTER SWEEP LOOP ===
for s = 1:num_sweeps
    current_cutoff = cutoff_range(s);
    distance_cutoff = 1 - current_cutoff;
    
    % 1. CUT THE HIERARCHICAL TREE
    cluster_idx = cluster(clustering_tree, 'cutoff', distance_cutoff, 'criterion', 'distance');
    num_total_clusters = max(cluster_idx);
    sweep_num_clusters(s) = num_total_clusters;
    
    % 2. INITIALIZE SPECTRA & ELUTIONS
    consensus_spectra = zeros(num_total_clusters, num_spectral_points);
    consensus_elutions = zeros(num_total_clusters, num_elution_points);
    cluster_sizes = zeros(num_total_clusters, 1);
    best_member_selected = zeros(num_total_clusters, 1);
    
    for c = 1:num_total_clusters
        members = find(cluster_idx == c);
        cluster_sizes(c) = length(members);
        
        cluster_spectra = Vtrn4_matrix(members, :);
        if cluster_sizes(c) > 1
            consensus_spectra(c, :) = median(cluster_spectra, 1);
            best_member_selected(c) = 0; 
        else
            consensus_spectra(c, :) = cluster_spectra; 
            best_member_selected(c) = members(1); 
        end
        
        cluster_elutions = zeros(cluster_sizes(c), num_elution_points);
        for m_idx = 1:cluster_sizes(c)
            node_idx = members(m_idx);
            padded_elution = zeros(1, num_elution_points);
            padded_elution(nodes_filtered_reordered(node_idx).fracs) = ...
                nodes_filtered_reordered(node_idx).optimized_elution * ...
                nodes_filtered_reordered(node_idx).singular_value;
            cluster_elutions(m_idx, :) = padded_elution;
        end
        
        if cluster_sizes(c) > 1
            consensus_elutions(c, :) = median(cluster_elutions, 1);
        else
            consensus_elutions(c, :) = cluster_elutions;
        end
    end
    
    % 3. GLOBAL OPTIMIZATION LOOP (Iterative Alternating Selection)
    for iter = 1:max_iters
        members_changed = false;
        for c = 1:num_total_clusters
            members = find(cluster_idx == c);
            if length(members) == 1
                continue;
            end
            
            best_err = inf;
            best_elution_for_c = zeros(1, num_elution_points);
            best_m = 0;
            
            for m_idx = 1:length(members)
                node_idx = members(m_idx);
                test_elution = zeros(1, num_elution_points);
                test_elution(nodes_filtered_reordered(node_idx).fracs) = ...
                    nodes_filtered_reordered(node_idx).optimized_elution * ...
                    nodes_filtered_reordered(node_idx).singular_value;
                
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
    
    % Uncomment if want to use consensus spectra to winning node
    % 4. SNAP CONSENSUS SPECTRA TO WINNING NODE
    % for c = 1:num_total_clusters
    %     winning_node = best_member_selected(c);
    %     if winning_node > 0 
    %         consensus_spectra(c, :) = Vtrn4_matrix(winning_node, :);
    %     end
    % end
    
    % 5. PURITY CONFIDENCE SCORING & RANKING
    score_covar = zeros(num_total_clusters, 1);
    score_sim = zeros(num_total_clusters, 1);
    score_intensity = zeros(num_total_clusters, 1);
    
    for c = 1:num_total_clusters
        members = find(cluster_idx == c);
        num_members = length(members);
        
        score_covar(c) = median(spec_correlations_filt(members));
        if num_members > 1
            sub_corr = corr_matrix(members, members);
            pairwise_vals = sub_corr(triu(true(num_members), 1));
            score_sim(c) = mean(pairwise_vals);
        else
            score_sim(c) = 1.0;
        end
        temp_intensities = zeros(num_members, 1);
        for m = 1:num_members
            temp_intensities(m) = max(nodes_filtered_reordered(members(m)).peaks);
        end
        score_intensity(c) = median(temp_intensities);
    end
    
    score_covar(score_covar < 0) = 0;
    score_sim(score_sim < 0) = 0;     
    max_intensity_overall = max(score_intensity);
    if max_intensity_overall > 0
        score_intensity = score_intensity ./ max_intensity_overall;
    end
    
    total_weight = weight_covar + weight_sim + weight_intensity;
    purity_score = ((score_covar .* weight_covar) + ...
                    (score_sim .* weight_sim) + ...
                    (score_intensity .* weight_intensity)) ./ total_weight;
    
    [sorted_scores, rank_idx] = sort(purity_score, 'descend');
    
    % 6. INCREMENTAL RECONSTRUCTION (Find Optimal k)
    reordered_spectra = consensus_spectra(rank_idx, :);
    reordered_elutions = consensus_elutions(rank_idx, :);
    
    cumulative_rel_error = zeros(num_total_clusters, 1);
    for k = 1:num_total_clusters
        partial_recon = reordered_elutions(1:k, :)' * reordered_spectra(1:k, :);
        cumulative_rel_error(k) = (norm(lib - partial_recon, 'fro') / lib_norm) * 100;
    end
    
    [min_error, optimal_k] = min(cumulative_rel_error);
    
    % Record sweep results
    sweep_min_errors(s) = min_error;
    sweep_optimal_ks(s) = optimal_k;
    
    fprintf('Cutoff %.2f | Clusters: %3d | Opt Components: %3d | Rel Error: %.2f%%\n', ...
            current_cutoff, num_total_clusters, optimal_k, min_error);
            
    % 7. SAVE STATE IF BEST OVERALL
    if min_error < global_best_error
        global_best_error = min_error;
        
        best_state.cutoff = current_cutoff;
        best_state.cluster_idx = cluster_idx;
        best_state.num_total_clusters = num_total_clusters;
        best_state.consensus_spectra = consensus_spectra;
        best_state.consensus_elutions = consensus_elutions;
        best_state.best_member_selected = best_member_selected;
        best_state.cluster_sizes = cluster_sizes;
        best_state.rank_idx = rank_idx;
        best_state.optimal_k = optimal_k;
        best_state.purity_score = purity_score;
        best_state.score_covar = score_covar;
        best_state.score_sim = score_sim;
        best_state.score_intensity = score_intensity;
        best_state.optimal_reconstructed_lib = reordered_elutions(1:optimal_k, :)' * reordered_spectra(1:optimal_k, :);
        best_state.cumulative_rel_error = cumulative_rel_error;
    end
end

fprintf('\n>>> DEEP SWEEP COMPLETE <<<\n');
fprintf('Winning Cutoff: r >= %.2f (Error: %.2f%% using %d components)\n\n', ...
        best_state.cutoff, global_best_error, best_state.optimal_k);

%
% ==========================================================
% --- VISUALIZATION 1: THE SWEEP CURVE ---
% ==========================================================
figure('Name', 'Cutoff Sweep Results', 'Position', [100, 100, 800, 500], 'Color', 'w');
yyaxis left;
plot(cutoff_range, sweep_min_errors, '-o', 'LineWidth', 2, 'MarkerSize', 6);
ylabel('Minimum Relative Error (%)', 'FontWeight', 'bold');
yline(global_best_error, 'r--', 'LineWidth', 1.5);

yyaxis right;
plot(cutoff_range, sweep_optimal_ks, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'Color', '#7E2F8E');
ylabel('Optimal Components (k)', 'FontWeight', 'bold');

title('Deep Sweep: Correlation Cutoff Optimization', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Correlation Cutoff (r)', 'FontSize', 12, 'FontWeight', 'bold');
grid on; box on;
xline(best_state.cutoff, 'k--', sprintf('Best Cutoff = %.2f', best_state.cutoff), ...
      'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'LineWidth', 2);
legend('Min Error (%)', 'Global Minimum', 'Optimal Components (k)', 'Location', 'northwest');

%
% ==========================================================
% --- UNPACK BEST STATE FOR FINAL VISUALIZATIONS ---
% ==========================================================
% Extract winning variables to the workspace
cluster_idx = best_state.cluster_idx;
num_total_clusters = best_state.num_total_clusters;
consensus_spectra = best_state.consensus_spectra;
consensus_elutions = best_state.consensus_elutions;
best_member_selected = best_state.best_member_selected;
cluster_sizes = best_state.cluster_sizes;
rank_idx = best_state.rank_idx;
optimal_k = best_state.optimal_k;

%
% ==========================================================
% --- VISUALIZATION 2: DENDROGRAM FOR WINNING CUTOFF ---
% ==========================================================
distance_cutoff = 1 - best_state.cutoff;
figDendro = figure('Name', sprintf('Dendrogram - Optimal Cutoff %.2f', distance_cutoff), ...
                    'Position', [200, 150, 1000, 600], 'Color', 'w');
dendrogram(clustering_tree, 0); 
title(sprintf('Optimal Dendrogram Cutoff (r \\geq %.2f)', best_state.cutoff), ...
      'FontSize', 14, 'FontWeight', 'bold');
ylabel('Distance (1 - r)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Spectra Index', 'FontSize', 12, 'FontWeight', 'bold');
yline(distance_cutoff, 'r--', sprintf('Distance Cutoff = %.2f', distance_cutoff), ...
      'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'FontSize', 11);
box on;

%
% ==========================================================
% --- VISUALIZATION 3: BEFORE & AFTER CONSENSUS (WINNER) ---
% ==========================================================
[~, sorted_cluster_indices] = sort(cluster_sizes, 'descend');
num_to_plot = min(20, num_total_clusters);
top_clusters = sorted_cluster_indices(1:num_to_plot);

% SPECTRA
figConsensus = figure('Name', 'Winning Cutoff: Spectra', 'Position', [250, 250, 1200, 800], 'Color', 'w');
t_cons = tiledlayout(ceil(num_to_plot/2), min(2, num_to_plot), 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_cons, sprintf('Consensus Spectra (Winner: Cutoff %.2f)', best_state.cutoff), 'FontSize', 16, 'FontWeight', 'bold');

for k = 1:num_to_plot
    c_idx = top_clusters(k);
    members = find(cluster_idx == c_idx);
    nexttile; hold on;
    for m = 1:length(members)
        plot(ppm, Vtrn4_matrix(members(m), :), 'Color', [0.7 0.7 0.7 0.4], 'LineWidth', 1);
    end
    plot(ppm, consensus_spectra(c_idx, :), 'Color', '#D95319', 'LineWidth', 2.5);
    hold off;
    set(gca, 'XDir', 'rev'); 
    title(sprintf('Cluster %d (Selected Node: %d)', c_idx, best_member_selected(c_idx)), 'FontSize', 12);
    xlabel('ppm', 'FontWeight', 'bold'); ylabel('Intensity', 'FontWeight', 'bold');
    grid on; box on;
end

% ELUTIONS
figElutions = figure('Name', 'Winning Cutoff: Elutions', 'Position', [300, 300, 1200, 800], 'Color', 'w');
t_elut = tiledlayout(ceil(num_to_plot/2), min(2, num_to_plot), 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_elut, sprintf('Consensus Elutions (Winner: Cutoff %.2f)', best_state.cutoff), 'FontSize', 16, 'FontWeight', 'bold');

for k = 1:num_to_plot
    c_idx = top_clusters(k);
    members = find(cluster_idx == c_idx);
    nexttile; hold on;
    for m_idx = 1:length(members)
        node_idx = members(m_idx);
        padded_elution = zeros(1, num_elution_points);
        padded_elution(nodes_filtered_reordered(node_idx).fracs) = ...
            nodes_filtered_reordered(node_idx).optimized_elution * ...
            nodes_filtered_reordered(node_idx).singular_value;
        plot(padded_elution, 'Color', [0.7 0.7 0.7 0.4], 'LineWidth', 1);
    end
    plot(consensus_elutions(c_idx, :), 'Color', '#0072BD', 'LineWidth', 2.5);
    hold off;
    title(sprintf('Cluster %d (Selected Node: %d)', c_idx, best_member_selected(c_idx)), 'FontSize', 12);
    xlabel('Fraction / Index', 'FontWeight', 'bold'); ylabel('Intensity (S-Scaled)', 'FontWeight', 'bold');
    grid on; box on;
end

%
% ==========================================================
% --- VISUALIZATION 4: PURITY RANKING & ERROR CURVE ---
% ==========================================================
figScores = figure('Name', 'Optimal Purity & Reconstruct Error', 'Position', [150, 150, 1200, 600], 'Color', 'w');
t_scores = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Bar Chart
nexttile;
num_bars = min(20, num_total_clusters); 
total_weight = weight_covar + weight_sim + weight_intensity;
bar_data = [best_state.score_covar(rank_idx(1:num_bars)) .* (weight_covar / total_weight), ...
            best_state.score_sim(rank_idx(1:num_bars)) .* (weight_sim / total_weight), ...
            best_state.score_intensity(rank_idx(1:num_bars)) .* (weight_intensity / total_weight)];

b = bar(1:num_bars, bar_data, 'stacked');
b(1).FaceColor = '#0072BD'; b(2).FaceColor = '#D95319'; b(3).FaceColor = '#EDB120';
title(sprintf('Top Ranked Components (Cutoff %.2f)', best_state.cutoff), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Cluster Rank', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Purity Score', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 1.05]); xticks(1:num_bars);
xticklabels(arrayfun(@num2str, rank_idx(1:num_bars), 'UniformOutput', false));
legend({'Covar Match', 'Similarity', 'Intensity'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
grid on;

% Error Curve
nexttile;
plot(1:num_total_clusters, best_state.cumulative_rel_error, 'LineWidth', 2.5, 'Color', '#0072BD');
hold on;
plot(optimal_k, global_best_error, 'p', 'MarkerSize', 15, 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'k');
xline(optimal_k, 'r--', sprintf('Optimal Cutoff = %d', optimal_k), ...
      'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5, 'FontSize', 11);
hold off;
title('Incremental Reconstruction Accuracy', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Number of Ranked Components Included', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Error (%)', 'FontSize', 12, 'FontWeight', 'bold');
grid on; box on;

%
% ==========================================================
% --- VISUALIZATION 5: 3D SPECTRA RECONSTRUCTION ---
% ==========================================================
displaySpectra3D(consensus_spectra(rank_idx,:), ppm);
title('Ranked Consensus Spectra (Winner)', 'FontSize', 14);

displaySpectra3D(lib, ppm);
title('Original NMR Fraction Library', 'FontSize', 14);

displaySpectra3D(best_state.optimal_reconstructed_lib, ppm);
title(sprintf('Optimal Reconstruction (Top %d Components) - Error: %.2f%%', optimal_k, global_best_error), 'FontSize', 14);



% Cleanup sweep variables
clear s current_cutoff distance_cutoff num_sweeps c members num_members m_idx node_idx iter max_iters ...
      members_changed best_err best_elution_for_c best_m test_elution reconstructed_lib err winning_node ...
      score_covar score_sim score_intensity temp_intensities m sub_corr pairwise_vals max_intensity_overall ...
      total_weight purity_score sorted_scores reordered_spectra reordered_elutions cumulative_rel_error ...
      k partial_recon min_error num_to_plot top_clusters sorted_cluster_indices c_idx padded_elution ...
      figDendro figConsensus figElutions figScores t_cons t_elut t_scores num_bars bar_data b ...
      AmaxScl  cluster* cutoff_range ...
      num* optimal_k sweep* weight* ;

%% Make figure to plot the variance explained by scree plot
% 
% ==========================================================
% --- VARIANCE EXPLAINED BY INTENSITY-SORTED COMPONENTS ---
% ==========================================================
disp('Calculating cumulative variance explained up to optimal k...');

% 1. Extract variables from the struct for cleaner code
opt_k    = best_state.optimal_k;
spectra  = best_state.consensus_spectra;
elutions = best_state.consensus_elutions;

% 2. Sort components by highest elution intensity (descending)
% Find the maximum intensity for each elution profile across all fractions
elution_max_ints = max(elutions, [], 2);

% Get the sorting indices in descending order
[~, sort_idx] = sort(elution_max_ints, 'descend');

% Reorder the spectra and elutions matrices based on this intensity ranking
spectra  = spectra(sort_idx, :);
elutions = elutions(sort_idx, :);

% Ensure we don't loop beyond available dimensions if opt_k was manually altered
max_k = min(size(spectra, 1), opt_k);

% 3. Calculate the total squared norm of the original library
lib_norm_sq = norm(lib, 'fro')^2;

% Preallocate the array for variance explained
variance_explained_pct = zeros(max_k, 1);

% 4. Incrementally rebuild the library and calculate variance explained
for k = 1:max_k
    % Reconstruct using only the top 'k' intensity-sorted components
    partial_reconstruction = elutions(1:k, :)' * spectra(1:k, :);
    
    % Calculate the squared norm of the residual (unexplained variance)
    residual_norm_sq = norm(lib - partial_reconstruction, 'fro')^2;
    
    % Variance Explained = 1 - (Residual Variance / Total Variance)
    variance_explained_pct(k) = (1 - (residual_norm_sq / lib_norm_sq)) * 100;
end

% ==========================================================
% --- VISUALIZATION: PUBLICATION QUALITY PLOT ---
% ==========================================================

figVariance = figure('Name', 'Cumulative Variance Explained (Intensity Sorted)', ...
                     'Position', [150, 150, 600, 600], 'Color', 'w');
hold on;

% 1. Calculate Y-Limits early so your vertical lines know how far to stretch
min_y = max(0, min(variance_explained_pct) - 5);
max_y = min(100, max(variance_explained_pct) + 5);

% 2. Draw vertical lines in the BACKGROUND using standard 'plot'
for i = [5, 10, 25]
    % plot([x1 x2], [y1 y2])
    plot([i, i], [min_y, max_y], ':', 'Color', 'k', 'LineWidth', 1.5); 
end

% 3. Plot the curve ON TOP of the dotted lines
% The 'MarkerFaceColor', 'w' will now properly mask the lines behind it
plot(1:max_k, variance_explained_pct, '-o', 'LineWidth', 2.5, 'MarkerSize', 7, ...
     'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'k');

% 4. Aesthetic tweaks
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on', ...
         'TickDir', 'out', 'Layer', 'top'); 

xlim([0, max_k + 1]);
ylim([min_y, max_y]);

hold off;


xlim([0, max_k + 1]);

% Dynamic Y-Limits based on the data range
min_y = max(0, min(variance_explained_pct) - 5);
max_y = min(100, max(variance_explained_pct) + 5);
ylim([min_y, max_y]);

fprintf('Variance explained at optimal cutoff (k = %d): %.2f%%\n', max_k, variance_explained_pct(max_k));

% ==========================================================
% --- CLEANUP TEMPORARY VARIABLES ---
% ==========================================================
clear opt_k spectra elutions elution_max_ints sort_idx max_k lib_norm_sq ...
      variance_explained_pct k partial_reconstruction residual_norm_sq ...
      figVariance min_y max_y;
%% Reconstruct library with different number of elements
% ==========================================================
% --- VARIANCE EXPLAINED BY INTENSITY-SORTED COMPONENTS ---
% ==========================================================
disp('Calculating cumulative variance explained up to optimal k...');

% 1. Extract variables from the struct for cleaner code
opt_k    = best_state.optimal_k;
spectra  = best_state.consensus_spectra;
elutions = best_state.consensus_elutions;

% 2. Sort components by highest elution intensity (descending)
% Find the maximum intensity for each elution profile across all fractions
elution_max_ints = max(elutions, [], 2);

% Get the sorting indices in descending order
[~, sort_idx] = sort(elution_max_ints, 'descend');

% Reorder the spectra and elutions matrices based on this intensity ranking
spectra  = spectra(sort_idx, :);
elutions = elutions(sort_idx, :);

% Ensure we don't loop beyond available dimensions if opt_k was manually altered
max_k = min(size(spectra, 1), opt_k);

% ==========================================================
% --- VISUALIZATION: PUBLICATION QUALITY PLOT ---
% ==========================================================

displaySpectra3D(elutions(1:5, :)' * spectra(1:5, :),ppm)
title("5")
displaySpectra3D(elutions(1:10, :)' * spectra(1:10, :),ppm)
title("10")
displaySpectra3D(elutions(1:25, :)' * spectra(1:25, :),ppm)
title("25")
displaySpectra3D(elutions(1:opt_k, :)' * spectra(1:opt_k, :),ppm)
title("Full Reconstruct")
displaySpectra3D(lib,ppm)
title("Raw")

clear opt_k spectra elutions elution_max_ints sort_idx max_k lib_norm_sq ...
      variance_explained_pct k partial_reconstruction residual_norm_sq ...
      figVariance min_y max_y;
%%%
%% ========================================================================
% HELPER FUNCTIONS (Must remain at the very bottom of your .m file)
% =========================================================================
%% Hybrid Target-Driven Elution Optimization (The Initializer)
function wt = optimizeElution1(Usub, peakIn, row, Amax)
    % --- DUAL DIMENSION AWARENESS ---
    n_svd = size(Usub, 2);     
    n_phys = size(peakIn, 2);  
    
    wt = zeros(n_svd, n_phys); 
    
    options = optimoptions('fmincon', 'MaxIterations', 10000, ...
                           'MaxFunctionEvaluations', 5000, 'Display', 'off',...
                           'Algorithm','interior-point');

    for j = 1:n_phys
        a = optimvar('a', 1, n_svd);
        target_profile = peakIn(:, j)'; 

        % --- THE INITIALIZATION FIX ---
        % 1. Project the ideal shape onto SVD components
        raw_guess = target_profile * Usub; 
        
        % 2. Calculate physical shape and force it to be legal
        physical_guess = raw_guess * Usub';
        physical_guess(physical_guess < 0) = 0; % Floor negatives
        
        % 3. Force normalization (L2 = 1.0)
        Scl_guess = sum(physical_guess.^2);
        if Scl_guess > 0
            physical_guess = physical_guess ./ sqrt(Scl_guess);
        end
        
        % 4. Project the clean, legal shape back to SVD space for the starting point
        x0.a = physical_guess * Usub; 
        
        prob = optimproblem;
        pred_profile = a * Usub';
        
        % --- THE BASELINE IRON ---
        error_weights = ones(1, length(target_profile));
        error_weights(target_profile == 0) = 100; 
        global_error = sum(error_weights .* (pred_profile - target_profile).^2);
        
        % --- THE SPEED UPGRADE: SOFT NON-NEGATIVITY PENALTY ---
        % Replaces 1000+ hard constraints with one smooth math equation
        epsilon = 1e-6;
        smooth_abs = sqrt(pred_profile.^2 + epsilon);
        negative_violations = (-pred_profile + smooth_abs) / 2;
        neg_penalty = 5000 * sum(negative_violations.^2);

        % Combine for the final objective
        prob.Objective = global_error + neg_penalty;
        
        % --- CONSTRAINTS ---
        prob.Constraints.cons1 = a * Usub(row(j), :)' >= Amax(j);
        prob.Constraints.cons2 = (dot(pred_profile, pred_profile) - 1)^2 <= 0.005;

        % --- THE IRON CEILING (AUC Fix) ---
        baseline_mask = target_profile == 0;
        max_allowed_area = 0.005 * max(target_profile) * sum(baseline_mask);
        prob.Constraints.cons4 = sum(pred_profile(baseline_mask)) <= max_allowed_area;

        % Solve
        sol = solve(prob, x0, 'Options', options);
        wt(:, j) = sol.a';
    end
end
%% Hybrid Target-Driven Elution Optimization with Tikhonov Regularization
function wt = optimizeElution(Usub, peakIn, row, Amax, wt2, SclS)
    % --- DUAL DIMENSION AWARENESS ---
    n_svd = size(Usub, 2);     
    n_phys = size(peakIn, 2);  
    
    wt = zeros(n_svd, n_phys);
    
    options = optimoptions('fmincon', 'MaxIterations', 10000, ...
                           'MaxFunctionEvaluations', 5000, 'Display', 'off',...
                           'Algorithm','interior-point');

    for j = 1:n_phys
        a = optimvar('a', 1, n_svd);
        
        % Inherit the previous iteration's rectangular weights
        x0.a = wt2(:, j)';
        
        prob = optimproblem;
        pred_profile = a * Usub'; 
        target_profile = peakIn(:, j)'; 
        
        % --- THE BASELINE IRON ---
        error_weights = ones(1, length(target_profile));
        error_weights(target_profile == 0) = 100; 
        global_error = sum(error_weights .* (pred_profile - target_profile).^2);
        
        % Tikhonov Regularization (ALS Stability)
        tikhonov_penalty = SclS * sum((a - wt2(:, j)').^2);
        
        % --- THE SPEED UPGRADE: SOFT NON-NEGATIVITY PENALTY ---
        epsilon = 1e-6;
        smooth_abs = sqrt(pred_profile.^2 + epsilon);
        negative_violations = (-pred_profile + smooth_abs) / 2;
        neg_penalty = 5000 * sum(negative_violations.^2);
        
        % Combine all forces
        prob.Objective = global_error + tikhonov_penalty + neg_penalty;
        
        % --- CONSTRAINTS ---
        prob.Constraints.cons1 = a * Usub(row(j), :)' >= Amax(j);
        
        proj = a * Usub'; 
        prob.Constraints.cons2 = (dot(proj, proj) - 1)^2 <= 0.005;

        % --- THE IRON CEILING (AUC Fix) ---
        baseline_mask = target_profile == 0;
        max_allowed_area = 0.005 * max(target_profile) * sum(baseline_mask);
        prob.Constraints.cons4 = sum(pred_profile(baseline_mask)) <= max_allowed_area;
   
        % Solve
        sol = solve(prob, x0, 'Options', options);
        wt(:, j) = sol.a';
    end
end
%%
%{
    Objective - Minimizes the Sum of Squared Errors against a target covariance spectrum.
    Tether - Tikhonov Regularization (anchors to elution weights).
    Penalties (Soft Constraints) - Non-negativity & Sparsity.
    Constraint - Normalization (L2 = 1.0).
%}
function wt2 = optimizeSpectra(Vsub, specIn, wt, SclS) 
    % --- DUAL DIMENSION AWARENESS ---
    n_svd = size(wt, 1);    
    n_phys = size(wt, 2);   

    wt2 = zeros(n_svd, n_phys);

    % --- Optimization options (fmincon required for the normalization constraint) ---
    options = optimoptions('fmincon', 'MaxIterations', 10000, ...
                           'MaxFunctionEvaluations', 5000, 'Display', 'off', ...
                           'Algorithm', 'interior-point');

    for j = 1:n_phys
        a = optimvar('a', 1, n_svd);
        x0.a = wt(:, j)'; 

        % Isolate the physical target spectrum for this specific component
        target_spec = specIn(j, :);

        % --- THE DEAD COMPONENT BYPASS ---
        % If the target spectrum has no variance (sum of absolute values is 0)
        if sum(abs(target_spec)) == 0
            fprintf('Notice: Component %d has an empty target spectrum. Bypassing optimization.\n', j);
            wt2(:, j) = wt(:, j); 
            continue; 
        end

        prob = optimproblem;

        % The predicted spectrum for this iteration (driven by the abstract SVD space)
        fullVec = a * Vsub; 

        % --- 1. THE NON-NEGATIVITY PENALTY (Soft Constraint) ---
        % Creates a smooth absolute value to isolate negative dips
        epsilon = 1e-6; 
        smooth_abs = sqrt(fullVec.^2 + epsilon);
        negative_violations = (-fullVec + smooth_abs) / 2;

        % We multiply the negative violations by a massive weight. 
        % This forces the baseline up without crashing the solver.
        neg_weight = 5000; 
        neg_penalty = neg_weight * sum(negative_violations.^2);

        % % --- 2. THE SPARSITY PENALTY (The Mixture Bleed Killer) ---
        % % By penalizing the sum of the absolute values of the whole spectrum, 
        % % we force the solver to prefer flat baselines and drop ghost peaks.
        % sparsity_weight = 5000; 
        % sparsity_penalty = sparsity_weight * sum(smooth_abs);

        % --- 2. TARGETED SPARSITY (The Smart Bleed Killer) ---
        % Identify where the STOCSY target is zero (the true baseline)
        baseline_mask = (target_spec == 0);
        
        % Apply a MASSIVE penalty (5000) only to the baseline regions to kill ghost peaks
        baseline_penalty = 5000 * sum(smooth_abs(baseline_mask));
        
        % Apply a GENTLE penalty (0.1) to the actual peak regions to preserve their natural shape
        peak_penalty = 0.1 * sum(smooth_abs(~baseline_mask));
        
        % Combine them
        sparsity_penalty = baseline_penalty + peak_penalty;

        %%%
        %%%

        % --- 3. THE SHAPE-TRACKING TARGET PULL ---
        % Minimizes the Sum of Squared Errors (SSE) between the math and the physics
        % Because both vectors are L2 normalized, the maximum possible error is strictly bounded!
        target_weight = 0.15; 
        target_pull = target_weight * sum((fullVec - target_spec).^2);

        % --- 4. THE TETHER ---
        tether_penalty = SclS * sum((a - wt(:,j)').^2);

        % --- 5. THE HYBRID OBJECTIVE FUNCTION ---
        prob.Objective = target_pull + tether_penalty + neg_penalty + sparsity_penalty; 

        % --- 6. HARD CONSTRAINT ---
        % Normalization (keeps the vector mathematically legal in SVD space)
        prob.Constraints.cons1 = (dot(fullVec, fullVec) - 1)^2 <= 0.005; 

        % Solve
        sol = solve(prob, x0, 'Options', options);
        wt2(:, j) = sol.a';
    end
end