%% --- MATCH CONSENSUS SPECTRA TO GISSMO LIBRARY ---
% The variable best_state from the SVD analysis must be in the workspace. 
% The GISSMO database is loaded on line 15. Provide the path to it.

% 1. SETTINGS & GISSMO FILE PATH
minFrac         = 0.0;  % Decimal fraction of highest peak for peak removal
shftMax         = 50;   % Max point shift each way (out of 16384 pts, in steps of 5)
plot_matches    = false; % Set to false to skip generating figures

fprintf('\n=========================================================================\n');
fprintf(' MATCHING CONSENSUS SPECTRA AGAINST GISSMO LIBRARY \n');
fprintf('=========================================================================\n');

% 2. LOAD GISSMO LIBRARY & EXTRACT TARGETS
load('GISSMO_Lib.mat');

spec_matrix = best_state.consensus_spectra;
ppm_target  = ppm(:)';     
ppm_lib     = ppm_lib(:)'; 

ns  = size(spec_matrix, 1);  
n   = size(lib_combined, 1);
pts = size(lib_combined, 2);

% L2 Normalization of the Reference GISSMO Library
lib_norms = vecnorm(lib_combined, 2, 2);
lib_norms(lib_norms == 0) = 1; % Prevent division by zero
lib_combined_norm = lib_combined ./ lib_norms;

% 3. INTERPOLATE & NORMALIZE CONSENSUS SPECTRA
disp('Aligning and normalizing consensus spectra onto GISSMO ppm axis...');
test_3 = zeros(pts, ns);

for i = 1:ns
    current_spec = spec_matrix(i, :);
    
    % Interpolate 65536-pt consensus spectrum onto 16384-pt GISSMO ppm_lib axis
    aligned_spec = interp1(ppm_target, current_spec, ppm_lib, 'linear', 0.0);
    
    % Remove peaks with < minFrac intensity
    MaxPk = max(aligned_spec);
    aligned_spec(aligned_spec < minFrac * MaxPk) = 0.0;
    
    % L2 Normalize
    scl = norm(aligned_spec);
    if scl > 0
        test_3(:, i) = aligned_spec ./ scl;
    else
        test_3(:, i) = 0.0;
    end
end

% 4. SCAN GISSMO LIBRARY
fprintf('Scanning %d GISSMO compounds for %d consensus spectra...\n', n, ns);
scoreHigh = zeros(ns, 1);
libBest   = zeros(ns, 1);
offsetL   = zeros(ns, 1);
offsetR   = zeros(ns, 1);

for i = 1:ns
    score = zeros(n, 1);
    [scoreHigh(i), libBest(i), offsetL(i), offsetR(i)] = ...
        libSearchOff(lib_combined_norm, test_3(:, i), n, score, shftMax);
end

% 5. OUTPUT TEXT SCORES & BUILD RESULTS TABLES
disp('-------------------------------------------------------------------------');
fprintf('Target # | GISSMO Match Name (Index)              | Score | offL | offR \n');
disp('-------------------------------------------------------------------------');

Target_ID        = (1:ns)';
Match_Name_Array = strings(ns, 1);

for i = 1:ns
    match_name          = combined_names(libBest(i));
    Match_Name_Array(i) = string(match_name);
    
    formatSpec = '   %2d    | %-30s (%4d) | %1.3f |  %3d |  %3d \n';
    fprintf(formatSpec, i, match_name, libBest(i), scoreHigh(i), offsetL(i), offsetR(i));
end
disp('-------------------------------------------------------------------------');

% Store unsorted (by Target_ID) and sorted (by Score) results tables
gissmo_results = table(Target_ID, Match_Name_Array, libBest, scoreHigh, offsetL, offsetR, ...
    'VariableNames', {'Target_ID', 'Match_Name', 'Match_Index', 'Score', 'Offset_L', 'Offset_R'});

sorted_gissmo_results = sortrows(gissmo_results, 'Score', 'descend');

% Display Ranked Summary
fprintf('\n=========================================================================\n');
disp(' TOP GISSMO MATCHES (RANKED BY SCORE) ');
fprintf('=========================================================================\n');
fprintf('Rank | Target | GISSMO Match Name              | Score | GISSMO Index \n');
disp('-------------------------------------------------------------------------');

for k = 1:height(sorted_gissmo_results)
    fprintf(' %2d  |   %2d   | %-30s | %1.3f |     %4d \n', ...
        k, sorted_gissmo_results.Target_ID(k), sorted_gissmo_results.Match_Name(k), ...
        sorted_gissmo_results.Score(k), sorted_gissmo_results.Match_Index(k));
end
disp('-------------------------------------------------------------------------');

% 6. AUTOMATED OUTPUT PLOTS
if plot_matches
    disp('Generating GISSMO match figures...');
    testOff = test_3;

    for id = 1:ns
        % Offset test spectrum leaving midpoint segment unchanged
        for ii = abs(offsetL(id))+1 : floor(pts/2)-abs(offsetL(id))
            testOff(ii, id) = test_3(ii - offsetL(id), id);
        end
        for ii = floor(pts/2)+abs(offsetR(id)) : pts-abs(offsetR(id))
            testOff(ii, id) = test_3(ii - offsetR(id), id);
        end
        
        name3      = ['Target-', num2str(id)];
        match_name = combined_names(libBest(id));
        
        fig_name = sprintf('GISSMO Match: Target %d (%s)', id, match_name);
        fig = figure('Name', fig_name, 'Color', 'w', ...
                     'Units', 'inches', 'Position', [2, 2, 7.5, 5]);
        
        p1 = plot(ppm_lib, lib_combined_norm(libBest(id), :), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5);
        hold on;
        p2 = plot(ppm_lib, testOff(:, id), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
        
        ax = gca;
        ax.XDir = 'reverse';
        ax.FontName = 'Arial';
        ax.FontSize = 14;
        ax.LineWidth = 1.5;
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.XColor = 'k';
        ax.YColor = 'k';
        
        xlabel('Chemical Shift (ppm)', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Intensity (a.u.)', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
        
        title_str = sprintf('%s vs. %s (#%d) | Score: %.3f', name3, match_name, libBest(id), scoreHigh(id));
        title(title_str, 'Interpreter', 'none', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
        
        legend_str = {sprintf('GISSMO: %s', match_name), 'Consensus Spectrum'};
        lgd = legend([p1, p2], legend_str, 'Location', 'northwest', 'Interpreter', 'none');
        lgd.FontName = 'Arial';
        lgd.FontSize = 14;
        legend('boxoff');
        
        hold off;
        drawnow;
    end
end

% 7. EXPLICIT WORKSPACE CLEANUP
% Explicitly removes only the variables created or loaded within this section
clear gissmo_mat_path minFrac shftMax plot_matches ...
      combined_names lib_combined ppm_lib ...
      spec_matrix ppm_target ns n pts lib_norms lib_combined_norm ...
      test_3 i current_spec aligned_spec MaxPk scl ...
      scoreHigh libBest offsetL offsetR score ...
      Target_ID Match_Name_Array match_name formatSpec k ...
      testOff id ii name3 fig_name fig p1 p2 ax title_str legend_str lgd;

disp('GISSMO matching complete. Temporary section variables cleared.');
%%
function [scoreHigh, libBest, offsetL, offsetR] = libSearchMid(lib_combined, test_3, n, score, shftMax)
    scoreHigh = -1.0;
    libBest   = 1;
    offsetL   = 0;
    offsetR   = 0;
    szh = size(test_3, 1) / 2; 
    testShft = zeros(2 * szh, 1);
    for j = 1:5:2*shftMax  
        testShft(j:szh-2*shftMax+j) = test_3(shftMax:szh-shftMax); 
        for k = 1:5:2*shftMax
            testShft(szh+k:2*szh-2*shftMax+k) = test_3(szh+shftMax:2*szh-shftMax); 
            for i = 1:n
                score(i) = dot(lib_combined(i, :)', testShft);
            end
            [scoreMax, scoreIdx] = max(score);
            if scoreMax > scoreHigh 
                scoreHigh = scoreMax;
                libBest   = scoreIdx;
                offsetL   = j - shftMax;
                offsetR   = k - shftMax;
            end
        end
    end
end

function [scoreHigh, libBest, offsetL, offsetR] = libSearchOff(lib_combined, test_3, n, score, shftMax)
    scoreHigh = -1.0;
    libBest   = 1;
    offsetL   = 0;
    offsetR   = 0;
    sz = size(test_3, 1); 
    testShft = zeros(size(test_3));
    for j = 1:5:2*shftMax
        testShft(j:sz-2*shftMax+j) = test_3(shftMax:sz-shftMax); 
        for i = 1:n
            score(i) = dot(lib_combined(i, :)', testShft);
        end
        [scoreMax, scoreIdx] = max(score);
        if scoreMax > scoreHigh 
            scoreHigh = scoreMax;
            libBest   = scoreIdx;
            offsetL   = j - shftMax;
            offsetR   = j - shftMax;
        end
    end
end

