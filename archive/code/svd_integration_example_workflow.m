%% Archived Automated & Manual SVD Component Integration Workflow
% =========================================================================
% Reproduces the manuscript qNMR concentration quantification comparing
% automated SVD integration (via integrateComponentsSVD) against manual
% single-peak integration across 23 metabolites using the archived
% manuscript decomposition (archive/result/output.mat) and raw spectral
% library with DSS intact (test_data/crassa.mat).
%
% Pipeline Stages:
%   1. Load archived best_state & ppm (archive/result/output.mat) and raw
%      spectral matrix X with DSS intact (test_data/crassa.mat)
%   2. Quantify all 23 metabolites via integrateComponentsSVD
%   3. Compute manual peak-bound integration benchmark and R^2
%   4. Export concentration summary CSV and 300-DPI PNG grouped bar chart
%      under cfg.results_dir (default: archive/result)
%
% External Functions Required:
%   integrateComponentsSVD
%
% Required MATLAB Toolboxes:
%   None (compatible with MATLAB R2019b or later)
% =========================================================================

%% 1. VERIFY PATH SETUP
% Run from the repository root directory (the folder containing setup.m and src/).
if isempty(which('integrateComponentsSVD'))
    if isfile('setup.m')
        setup;
    else
        error('svd_integration_example_workflow:PathNotSet', ...
            'Run this workflow from the repository root directory (the folder containing setup.m and src/).');
    end
else
    disp('Repository is on MATLAB path and archived SVD integration workflow is ready.');
end

%% 2. USER CONFIGURATION
% Edit the parameters below to match your dataset and metabolite assignments.

% --- Output Directory & Exports ---
% Defaults to archive/result; set to any relative or absolute path if desired.
cfg.results_dir  = fullfile('archive', 'result');
cfg.save_outputs = true;      % Export metabolite concentrations CSV and PNG bar chart

% --- Input Data ---
% Archived manuscript best_state & ppm, paired with raw crassa.mat (DSS intact).
cfg.svd_results_file = fullfile('archive', 'result', 'output.mat');
cfg.raw_data_file    = fullfile('test_data', 'crassa.mat');

% --- Target Metabolites (23 quantified metabolites, Glucose #8 excluded) ---
cfg.metabolite_names = { ...
    'Trehalose', 'Alanine', 'Valine', 'Glutamic Acid', 'Succinic Acid', ...
    'Glycine', 'Arginine', 'Lactate', 'Glycerol', 'Glutamine', ...
    'Lysine', 'Leucine', 'GABA', 'Histidine', 'Malic Acid', ...
    'Asparagine', 'Proline', 'GSSG', 'Isoleucine', 'AMP', ...
    'KIV', '1GP', 'Uridine' ...
};

% Row indices in archive/result/output.mat best_state.consensus_spectra
cfg.component_idx = [ ...
    21, 73, 36, 20, 47, ...
    65, 28, 16, 17, 39, ...
    26, 57, 53, 71, 35, ...
    43, 68, 29, 74, 31, ...
    49, 15, 37 ...
];

% Total non-exchangeable protons per metabolite (full-spectrum SVD integration)
cfg.nuclei_analyte = [ ...
    14,  4,  8,  5,  4, ...
     2,  7,  4,  5,  5, ...
     9, 10,  6,  5,  3, ...
     3,  7, 20, 10,  7, ...
     7,  5,  8 ...
];

% --- Manual Single-Peak Integration Benchmark Parameters ---
cfg.compute_manual_benchmark = true;
cfg.manual_fractions = { ...
    120:124,  74:99,   74:82,   92:100,   4:16, ...
     93:100, 107:114,   2:14,   33:37,   99:106, ...
    111:119,  55:66,   68:73,  110:117,  27:32, ...
    103:108,  80:82,  134:136,  64:72,  129:133, ...
      4:8,   125:127,  46:49 ...
};
cfg.manual_ppm_bounds = [ ...
    5.11813, 5.22689; ... % Trehalose
    1.44946, 1.48810; ... % Alanine
    1.02350, 1.04258; ... % Valine
    2.29424, 2.38201; ... % Glutamic Acid
    2.38440, 2.40252; ... % Succinic Acid
    3.54100, 3.55450; ... % Glycine
    1.61117, 1.67700; ... % Arginine
    1.28823, 1.34786; ... % Lactate
    3.74625, 3.79682; ... % Glycerol
    2.40062, 2.48982; ... % Glutamine
    2.99258, 3.02931; ... % Lysine
    1.65267, 1.77621; ... % Leucine
    2.25847, 2.31809; ... % GABA
    7.79223, 7.86188; ... % Histidine
    2.62719, 2.68921; ... % Malic Acid
    2.81991, 2.87667; ... % Asparagine
    4.09829, 4.14265; ... % Proline
    2.47551, 2.59905; ... % GSSG
    0.98915, 1.01348; ... % Isoleucine
    8.57262, 8.63034; ... % AMP
    1.10220, 1.12367; ... % KIV
    3.64704, 3.69760; ... % 1GP
    7.83516, 7.89431  ... % Uridine
];
cfg.manual_nuclei = [ ...
    2, 3, 3, 2, 4, ...
    2, 1, 3, 1, 2, ...
    2, 3, 2, 1, 1, ...
    1, 1, 4, 3, 1, ...
    6, 2, 1 ...
];

% --- Internal Reference Standard (DSS) & Fraction Settings ---
cfg.conc_ref   = 0.333;          % Known DSS concentration (mM)
cfg.ref_bounds = [-0.02, 0.02];  % DSS peak integration window (ppm)
cfg.nuclei_ref = 9;              % Number of protons for DSS methyl singlet
cfg.fractions  = 1:140;          % Fractions to sum over for SVD integration
cfg.verbose    = true;           % Print summary to Command Window

% --- Visualization Options ---
cfg.show_conc_plot = true;       % Plot publication-style log-scale bar chart

%% 3. LOAD & VALIDATE INPUTS
total_timer = tic;

if ~isfile(cfg.svd_results_file)
    error('svd_integration_example_workflow:MissingSVDResults', ...
        'Archived SVD results file not found: %s', cfg.svd_results_file);
end
fprintf('Loading archived SVD results from: %s\n', cfg.svd_results_file);
loaded_svd = load(cfg.svd_results_file);

if ~isfield(loaded_svd, 'best_state')
    error('svd_integration_example_workflow:MissingBestState', ...
        'MAT file must contain variable ''best_state''.');
end
best_state = loaded_svd.best_state;

if isfield(loaded_svd, 'ppm')
    ppm_target = loaded_svd.ppm;
elseif isfield(loaded_svd, 'ppm_work')
    ppm_target = loaded_svd.ppm_work;
else
    error('svd_integration_example_workflow:MissingPPM', ...
        'MAT file must contain variable ''ppm'' or ''ppm_work''.');
end

if ~isfile(cfg.raw_data_file)
    error('svd_integration_example_workflow:MissingRawData', ...
        'Raw spectral data file not found: %s', cfg.raw_data_file);
end
fprintf('Loading raw spectral data (with DSS intact) from: %s\n', cfg.raw_data_file);
loaded_raw = load(cfg.raw_data_file);

if isfield(loaded_raw, 'X')
    X_raw = loaded_raw.X;
elseif isfield(loaded_raw, 'X_work')
    X_raw = loaded_raw.X_work;
else
    error('svd_integration_example_workflow:MissingX', ...
        'Raw MAT file must contain variable ''X''.');
end

if size(X_raw, 2) ~= numel(ppm_target)
    error('svd_integration_example_workflow:DimensionMismatch', ...
        'Columns of X_raw (%d) do not match length of ppm (%d).', ...
        size(X_raw, 2), numel(ppm_target));
end

fprintf('Loaded raw spectra (%d fractions x %d ppm points) for %d target metabolites.\n', ...
    size(X_raw, 1), numel(ppm_target), numel(cfg.component_idx));

%% 4. AUTOMATED SVD INTEGRATION
fprintf('\n--- Running Automated SVD Integration ---\n');

[overall_conc_svd, fraction_concs_svd, integration_details] = integrateComponentsSVD( ...
    X_raw, ppm_target, best_state, cfg.component_idx, cfg.nuclei_analyte, ...
    cfg.conc_ref, cfg.ref_bounds, cfg.nuclei_ref, cfg.fractions, cfg.verbose);

%% 5. MANUAL INTEGRATION BENCHMARK & R^2
num_met = numel(cfg.component_idx);
if cfg.compute_manual_benchmark
    overall_conc_manual = zeros(1, num_met);
    [~, dss_idx1] = min(abs(ppm_target - cfg.ref_bounds(1)));
    [~, dss_idx2] = min(abs(ppm_target - cfg.ref_bounds(2)));
    dss_range = min(dss_idx1, dss_idx2) : max(dss_idx1, dss_idx2);

    for k = 1:num_met
        fracs_k = cfg.manual_fractions{k};
        [~, m_idx1] = min(abs(ppm_target - cfg.manual_ppm_bounds(k, 1)));
        [~, m_idx2] = min(abs(ppm_target - cfg.manual_ppm_bounds(k, 2)));
        met_range = min(m_idx1, m_idx2) : max(m_idx1, m_idx2);

        frac_concs_k = zeros(numel(fracs_k), 1);
        for f = 1:numel(fracs_k)
            row_f   = fracs_k(f);
            int_dss = trapz(X_raw(row_f, dss_range));
            int_met = trapz(X_raw(row_f, met_range));
            frac_concs_k(f) = (int_met / int_dss) * (cfg.nuclei_ref / cfg.manual_nuclei(k)) * cfg.conc_ref;
        end
        overall_conc_manual(k) = sum(frac_concs_k);
    end

    SS_res    = sum((overall_conc_manual(:) - overall_conc_svd(:)).^2);
    SS_tot    = sum((overall_conc_manual(:) - mean(overall_conc_manual(:))).^2);
    R_squared = 1 - (SS_res / SS_tot);
    fprintf('Manual vs. SVD Integration R^2 = %.4f\n\n', R_squared);

    integration_results = table( ...
        cfg.metabolite_names(:), ...
        cfg.component_idx(:), ...
        cfg.nuclei_analyte(:), ...
        overall_conc_manual(:), ...
        overall_conc_svd(:), ...
        'VariableNames', {'Metabolite', 'Component_ID', 'Nuclei_Analyte', ...
                          'Manual_Concentration_mM', 'SVD_Concentration_mM'});
else
    overall_conc_manual = [];
    R_squared = NaN;
    integration_results = table( ...
        cfg.metabolite_names(:), ...
        cfg.component_idx(:), ...
        cfg.nuclei_analyte(:), ...
        overall_conc_svd(:), ...
        'VariableNames', {'Metabolite', 'Component_ID', 'Nuclei_Analyte', 'SVD_Concentration_mM'});
end

disp(integration_results);

if cfg.save_outputs
    table_dir = fullfile(cfg.results_dir, 'tables');
    if ~isfolder(table_dir), mkdir(table_dir); end
    writetable(integration_results, fullfile(table_dir, 'metabolite_concentrations.csv'));
    fprintf('Saved metabolite concentrations table to: %s\n', ...
        fullfile(table_dir, 'metabolite_concentrations.csv'));
end

%% 6. PUBLICATION-STYLE CONCENTRATION COMPARISON BAR CHART
if cfg.show_conc_plot
    fig_conc = figure('Name', 'Comparison of Metabolite Concentrations', ...
                      'Color', 'w', 'Position', [100, 100, 1600, 800]);

    if cfg.compute_manual_benchmark
        plot_data_matrix = [overall_conc_manual(:), overall_conc_svd(:)];
        bar_handles = bar(plot_data_matrix, 'grouped', 'EdgeColor', 'k', 'LineWidth', 1.5);
        set(gca, 'YScale', 'log');
        bar_handles(1).BaseValue = 0.1;
        bar_handles(2).BaseValue = 0.1;
        bar_handles(1).FaceColor = [0.2, 0.2, 0.2];
        bar_handles(2).FaceColor = 'r';
        bar_handles(1).EdgeColor = 'k';
        bar_handles(2).EdgeColor = 'k';
    else
        bar_handle = bar(overall_conc_svd(:), 'FaceColor', 'r', ...
                         'EdgeColor', 'k', 'LineWidth', 1.5);
        set(gca, 'YScale', 'log');
        bar_handle.BaseValue = 0.1;
    end

    set(gca, 'LineWidth', 2.5, ...
             'FontSize', 16, ...
             'FontWeight', 'bold', ...
             'FontName', 'Carlito', ...
             'TickDir', 'out', ...
             'XColor', 'k', ...
             'YColor', 'k');

    xlim([0, num_met + 1]);
    xticks(1:num_met);
    xticklabels(cfg.metabolite_names);
    xtickangle(45);

    ylim([0.1, 600]);
    yticks([0.1, 1, 10, 100, 500]);
    yticklabels({'0.1', '1', '10', '100', '500'});
    box off;

    xlabel('Metabolite', 'FontSize', 20, 'FontWeight', 'bold', ...
           'FontName', 'Carlito', 'Color', 'k');
    ylabel('Concentration (mM)', 'FontSize', 20, 'FontWeight', 'bold', ...
           'FontName', 'Carlito', 'Color', 'k');

    if cfg.compute_manual_benchmark
        title('Comparison of Metabolite Concentrations', ...
              'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Carlito');
        legend_handle = legend('Manual Integration', 'SVD Integration', ...
            'Location', 'northeast', 'NumColumns', 1);
        legend_handle.Box        = 'off';
        legend_handle.FontSize   = 16;
        legend_handle.FontWeight = 'bold';
        legend_handle.FontName   = 'Carlito';

        text(0.8185, 0.88, sprintf('R^2 = %.3f', R_squared), ...
            'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', ...
            'FontName', 'Carlito', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'none');
    else
        title('Automated SVD Metabolite Concentrations', ...
              'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Carlito');
    end

    if cfg.save_outputs && isgraphics(fig_conc)
        fig_dir = fullfile(cfg.results_dir, 'figures', 'integration');
        if ~isfolder(fig_dir), mkdir(fig_dir); end
        print(fig_conc, fullfile(fig_dir, 'metabolite_concentrations.png'), '-dpng', '-r300');
        fprintf('Saved concentration bar chart to: %s\n', ...
            fullfile(fig_dir, 'metabolite_concentrations.png'));
    end
end

fprintf('\nTotal archived SVD integration workflow runtime: %.1f s.\n', toc(total_timer));
