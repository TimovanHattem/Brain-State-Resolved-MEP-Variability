clear;          % Clear workspace variables
close all;      % Close all figure windows
clc;            % Clear command window

% Get full path of current script and move to its directory
filePath = matlab.desktop.editor.getActiveFilename;
mainpath = fileparts(filePath);   % Adjust to project root folder structure
cd(mainpath);

%% ----------- USER SETTINGS -----------
region              = 'M1'; %M1
response            = 'MEP'; %MEP %TEP

var_method          = 'SD'; %SD %CV %MAD %IQR %CQV %TTV_NIU
bin_mean_normalization   = false;
MEP_feature         = 'p2p_microvolts'; %p2p_microvolts %p2p_log_microvolts %p2p_zscored_log
TEP_feature         = 'p2p';        %p2p %AUC %GMFP %LMFP
TEP_metric_type     = 'raw';  %raw %log %z_log
TEP_window          = 'allearly';   %allearly %N15P30 %P30N45 %N45P60

site                = 'Tuebingen'; %Tuebingen %Aalto %combined

prestim_feature     = 'Phase'; %Power %Connectivity %Phase %PowerPhase %PowerConnectivity
bands_of_interest   = {'alpha'}; %{'theta','alpha','beta', 'gamma'};
if strcmpi(prestim_feature, 'Phase') || strcmpi(prestim_feature, 'PowerPhase') || strcmpi(prestim_feature, 'PowerConnectivity')
    bands_of_interest = {'alpha'};
end
splits              = 'tails25'; %tails25 %tertiles %median %quintiles %octiles %deciles '4phasebins' '8phasebins'
nPhaseBins          = 4;
if strcmpi(prestim_feature, 'Phase') && nPhaseBins == 4
    splits          = '4phasebins';
elseif strcmpi(prestim_feature, 'Phase') && nPhaseBins == 8
    splits          = '8phasebins';
end
nPerState_general   = 100;
if nPhaseBins == 4
    nPerState_phase     = 100;
else
    nPerState_phase = 50;
end
nPerState_combos    = 25;
nSurrogate          = 50;

do_bootstrap        = true;
n_boot_iter         = 10;

save_output         = false;

generate_plots      = true;

plots_show          = 'VARrel'; %VARrel %VAR %MEAN_AMP
if strcmpi(plots_show, 'MEAN_AMP')
    MEP_feature =  'p2p_zscored_log';
    TEP_metric_type = 'z_log';
end

if strcmpi(var_method, 'CV')
    MEP_feature = 'p2p_microvolts';
elseif strcmpi(var_method, 'SD')
    MEP_feature = 'p2p_zscored_log';
end

rng(51)

%% ----------- SELECT DATA INPUT -----------

% % REFTEP
reftep = true; screen = false;
DATAIN_SPLIT = '';
DATAOUT_VIS = '';
DATAIN_EP = '';
all_subjects = load('');

%% ----------- EXTRACT HJORTH FEATURES -----------

if strcmpi(response, 'MEP')
    load(fullfile(DATAIN_EP, 'alldata_MEP_master.mat'), 'alldata_MEP_master');
elseif strcmpi(response, 'TEP')

    roi_channels = {'C3'};
    %roi_channels = {'C3', 'CP3', 'C1'};

    % Load your TEP data if not already in workspace
    if ~exist('alldata_sensor_TEP_master','var')
        load(fullfile(DATAIN_EP, 'alldata_sensor_TEP_master.mat'), 'alldata_sensor_TEP_master');
    end
    load(fullfile(DATAIN_EP, 'times_eeg_post.mat'))
    load(fullfile(DATAIN_EP, 'chanlocs_eeg.mat'))
    timeVec = times_eeg_post;
    chanlocs = chanlocs_eeg;

    load(fullfile(DATAIN_EP, 'alldata_MEP_master.mat'), 'alldata_MEP_master'); % for synchronizing rejected trials
end

%REFTEP_TUEBINGEN=1:29 %REFTEP_AALTO=30:51 %SCREEN3=1:51 %SCREEN5=52:65
if reftep
    if strcmpi(site, 'Tuebingen')
        sub_range = 1:29;
    elseif strcmpi(site, 'Aalto')
        sub_range = 30:51;
    elseif strcmpi(site, 'combined')
        sub_range = 1:51;
    end
else
    sub_range = 1:numel(all_subjects);
end

count = 0;
resultsTable = struct([]);

results.sig   = {};   % subjects WITH significant phase effect
results.nosig = {};   % subjects WITHOUT significant phase effect
alpha = 0.05;         % significance threshold
reftep_sig_phase_subject_list= {};
reftep_insig_phase_subject_list = {};
reftep_sig_power_subject_list= {};
reftep_insig_power_subject_list = {};
reftep_sig_connectivity_subject_list= {};
reftep_insig_connectivity_subject_list = {};
phase_effects = struct();
for subidx = sub_range
    subID = all_subjects{subidx};
    subids = strrep(subID, 'REFTEP_', '');
    subids = strrep(strrep(subids, 'SCREEN3_', ''), 'SCREEN5_', '');
    if contains(subID, 'SCREEN5_')
        subids = num2str(str2num(subids) + 500);
    end
    sub_num = horzcat('sub-', subids);

    % --------- Load response data (MEP/TEP) ---------
    if strcmpi(response, 'MEP')
        all_data = alldata_MEP_master(subidx).(MEP_feature);
    end

    % --------- Load entire splits table for feature ---------
    splitsTbl = load_splits_table(DATAIN_SPLIT, sub_num, 1, prestim_feature, splits);

    % --------- Variability quantification ---------
    Tb = splitsTbl; % ALPHA ONLY

    %% Analysis: Circular to linear regression

    if strcmpi(prestim_feature, 'Phase')
        addpath(genpath('D:\2025-08 T2TVAR\Analysis\toolboxes\CircStat2012a'))

        all_phase = Tb.Phase_C3H;    % enter phase values
        all_meps = all_data;     % enter MEPs

        % --- Circular-linear correlation ---
        [rho, pval_circ] = circ_corrcl(all_phase, all_meps);
        fprintf('Circular-linear correlation: rho = %.4f, p = %.4f\n', rho, pval_circ);

        % --- Circular-linear regression ---
        mdl_linear.estphase     = all_phase(:);
        mdl_linear.excitability = all_meps(:);
        mdl_linear.fit = fitlm([cos(mdl_linear.estphase), sin(mdl_linear.estphase)], ...
            mdl_linear.excitability);
        mdl_linear.pval       = mdl_linear.fit.coefTest;
        mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
        mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;
        fprintf('Regression F-test: p = %.4f\n', mdl_linear.pval);
        fprintf('R² = %.4f, R²adj = %.4f\n', mdl_linear.R2ordinary, mdl_linear.R2adjusted);

        % --- Regression line: y = b0 + b1*cos(x) + b2*sin(x) ---
        x = linspace(-pi, pi, 1000);
        y = mdl_linear.fit.Coefficients.Estimate(1) ...
            + mdl_linear.fit.Coefficients.Estimate(2) .* cos(x) ...
            + mdl_linear.fit.Coefficients.Estimate(3) .* sin(x);
        mdl_linear.opt = x(y == max(y));   % optimal phase (MEP maximum)
        fprintf('Optimal phase: %.4f rad (%.2f deg)\n', mdl_linear.opt, rad2deg(mdl_linear.opt));

        phase_effects(subidx).subject = subID;
        phase_effects(subidx).rho = rho;
        phase_effects(subidx).pval_circ = pval_circ;
        phase_effects(subidx).mdl_linear_pval = mdl_linear.pval;
        phase_effects(subidx).mdl_linear_R2ordinary =mdl_linear.R2ordinary;
        phase_effects(subidx).mdl_linear_R2adjusted =mdl_linear.R2adjusted;
        phase_effects(subidx).mdl_linear_opt = mdl_linear.opt;

        % if generate_plots
        %     % --- Plot ---
        %     figure; hold on;
        %     % Raw data
        %     scatter(all_phase, mdl_linear.excitability, 10, [0.7 0.7 0.7], 'filled', ...
        %         'HandleVisibility', 'off');
        %     % Regression line
        %     plot(x, y, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
        %     % Optimal phase marker
        %     xline(mdl_linear.opt, '-', {'Maximum MEP'}, ...
        %         'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'Layer', 'bottom');
        % 
        %     xlabel('Phase (rad)');
        %     ylabel('MEP amplitude (z-scored)');
        %     title(sprintf('Phase effect: p = %.4f, optimal phase = %.4f', mdl_linear.pval, mdl_linear.opt));
        %     xlim([-pi pi]);
        %     xticks([-pi -pi/2 0 pi/2 pi]);
        %     xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        %     legend({'regression fit', 'optimal phase'});
        % end

        % Save significant and insignificant phase effect subjects
        if mdl_linear.pval < alpha
            reftep_sig_phase_subject_list{end+1} = subID;
            fprintf('%-15s  SIGNIFICANT    p = %.4f, R²adj = %.4f\n', subID, mdl_linear.pval, mdl_linear.R2adjusted);
        else
            reftep_insig_phase_subject_list{end+1} = subID;
            fprintf('%-15s  not significant p = %.4f, R²adj = %.4f\n', subID, mdl_linear.pval, mdl_linear.R2adjusted);
        end

        % Save as high-quality PNG (300 dpi)
        if save_output
            exportgraphics(gcf, fullfile(DATAOUT_VIS, 'Phase', sprintf('%s_Phase_effect.png', subID)), 'Resolution', 300);
            save(fullfile('')
            save(fullfile('')
        end

    elseif strcmpi(prestim_feature, 'Power')

        all_power = Tb.Hjorth_C3_Power(strcmpi(Tb.Freq_band, 'alpha'));    % enter phase values

        all_meps = all_data;     % enter MEPs

        % --- Linear correlation ---
        [rho, pval] = corr(all_power, all_meps, 'Rows', 'complete', 'Type', 'Spearman');
        fprintf('Linear correlation: rho = %.4f, p = %.4f\n', rho, pval);

        % --- Linear regression ---
        mdl_linear.power = all_power(:);
        mdl_linear.excitability = all_meps(:);
        mdl_linear.fit = fitlm(mdl_linear.power, mdl_linear.excitability);
        mdl_linear.pval = mdl_linear.fit.Coefficients.pValue(2);
        mdl_linear.slope = mdl_linear.fit.Coefficients.Estimate(2);
        mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
        mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;

        % Save significant and insignificant phase effect subjects
        if pval < alpha
            reftep_sig_power_subject_list{end+1} = subID;
            fprintf('%-15s  SIGNIFICANT    p = %.4f\n', subID, pval);
        else
            reftep_insig_power_subject_list{end+1} = subID;
            fprintf('%-15s  not significant p = %.4f\n', subID, pval);
        end

        power_effects(subidx).subject = subID;
        power_effects(subidx).rho = rho;
        power_effects(subidx).pval = pval;
        power_effects(subidx).mdl_linear_slope = mdl_linear.slope;
        power_effects(subidx).mdl_linear_pval = mdl_linear.pval;
        power_effects(subidx).mdl_linear_R2ordinary =mdl_linear.R2ordinary;
        power_effects(subidx).mdl_linear_R2adjusted =mdl_linear.R2adjusted;


    elseif strcmpi(prestim_feature, 'Connectivity')

        all_connectivity = Tb.Hjorth_C3C4_Connectivity(strcmpi(Tb.Freq_band, 'alpha'));    % enter phase values

        all_meps = all_data;     % enter MEPs

        % --- Linear correlation ---
        [rho, pval] = corr(all_connectivity, all_meps, 'Rows', 'complete', 'Type', 'Spearman');
        fprintf('Linear correlation: rho = %.4f, p = %.4f\n', rho, pval);

        mdl_linear.connectivity = all_connectivity(:);
        mdl_linear.excitability = all_meps(:);
        mdl_linear.fit = fitlm(mdl_linear.connectivity, mdl_linear.excitability);
        mdl_linear.pval = mdl_linear.fit.Coefficients.pValue(2);
        mdl_linear.slope = mdl_linear.fit.Coefficients.Estimate(2);
        mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
        mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;

        % Save significant and insignificant phase effect subjects
        if pval < alpha
            reftep_sig_connectivity_subject_list{end+1} = subID;
            fprintf('%-15s  SIGNIFICANT    p = %.4f\n', subID, pval);
        else
            reftep_insig_connectivity_subject_list{end+1} = subID;
            fprintf('%-15s  not significant p = %.4f\n', subID, pval);
        end

        connectivity_effects(subidx).subject = subID;
        connectivity_effects(subidx).rho = rho;
        connectivity_effects(subidx).pval = pval;
        connectivity_effects(subidx).mdl_linear_slope = mdl_linear.slope;
        connectivity_effects(subidx).mdl_linear_pval = mdl_linear.pval;
        connectivity_effects(subidx).mdl_linear_R2ordinary =mdl_linear.R2ordinary;
        connectivity_effects(subidx).mdl_linear_R2adjusted =mdl_linear.R2adjusted;

    end

end

%% Plot
if generate_plots
    if strcmpi(prestim_feature, 'Phase')
        phase_data = readtable(['']);

        new_high = (phase_data.VARrel_phasebins_1);
        new_low = (phase_data.VARrel_phasebins_3);
        new_diff = new_high-new_low;

        [r, p_val] = corr([phase_effects.rho]', new_diff, 'Rows', 'complete', 'Type', 'Spearman');
        
        figure;
        scatter([phase_effects.rho]', new_diff, 60, 'filled', ...
            'MarkerFaceColor', [0.3 0.6 0.8], 'MarkerFaceAlpha', 0.7);
        hold on;

        % Fit and plot regression line
        p = polyfit([phase_effects.rho]', new_diff, 1);
        xFit = linspace(min([phase_effects.rho]'), max([phase_effects.rho]'), 100);
        yFit = polyval(p, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 2);

        hold off;

        %ylabel('\ - Low phase bins (VARrel)', 'FontWeight', 'bold', 'FontSize', 14);
        %xlabel('Mu phase effect strength (\rho)', 'FontWeight', 'bold', 'FontSize', 14);
        %set(gca, 'FontSize', 12, 'Box', 'off');

        % Add correlation coefficient and p-value as text on the plot
        text(0.35, 0.99, sprintf('r = %.3f, p = %.3f', r, p_val), ...
            'Units', 'normalized', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

        %ylim([-40 40])
        %xlim([-0.15, 0.25])
        box off;
        axis square
        %set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
        ax = gca;
        set(ax, 'FontSize', 14);
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.LineWidth = 2;
        fontsize(gcf, scale=1.5)
        fontname('Arial')

    end

    if strcmpi(prestim_feature, 'Power')
        power_data = readtable(['']);

        new_diff = power_data.DeltaSpread_highlow;

        [r, p_val] = corr([power_effects.rho]', new_diff, 'Rows', 'complete', 'Type', 'Spearman');
        
        figure;
        scatter([power_effects.rho]', new_diff, 60, 'filled', ...
            'MarkerFaceColor', [0.3 0.6 0.8], 'MarkerFaceAlpha', 0.7);
        hold on;

        % Fit and plot regression line
        p = polyfit([power_effects.rho]', new_diff, 1);
        xFit = linspace(min([power_effects.rho]'), max([power_effects.rho]'), 100);
        yFit = polyval(p, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 2);

        hold off;

        %ylabel('\ - Low phase bins (VARrel)', 'FontWeight', 'bold', 'FontSize', 14);
        %xlabel('Mu phase effect strength (\rho)', 'FontWeight', 'bold', 'FontSize', 14);
        %set(gca, 'FontSize', 12, 'Box', 'off');

 
        if p_val < 0.001
            p_val = 0.001;
            text(0.35, 0.99, sprintf('r = %.3f, p < %.3f', r, p_val), ...
                'Units', 'normalized', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
        else
            % Add correlation coefficient and p-value as text on the plot
            text(0.35, 0.99, sprintf('r = %.3f, p = %.3f', r, p_val), ...
                'Units', 'normalized', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
        end

        %ylim([-40 40])
        %xlim([-0.15, 0.25])
        box off;
        axis square
        %set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
        ax = gca;
        set(ax, 'FontSize', 14);
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.LineWidth = 2;
        fontsize(gcf, scale=1.5)
        fontname('Arial')
    end

    if strcmpi(prestim_feature, 'Connectivity')
        connectivity_data = readtable(['']);

        new_diff = connectivity_data.DeltaSpread_highlow;

        [r, p_val] = corr([connectivity_effects.rho]', new_diff, 'Rows', 'complete', 'Type', 'Spearman');
        
        figure;
        scatter([connectivity_effects.rho]', new_diff, 60, 'filled', ...
            'MarkerFaceColor', [0.3 0.6 0.8], 'MarkerFaceAlpha', 0.7);
        hold on;

        % Fit and plot regression line
        p = polyfit([connectivity_effects.rho]', new_diff, 1);
        xFit = linspace(min([connectivity_effects.rho]'), max([connectivity_effects.rho]'), 100);
        yFit = polyval(p, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 2);

        hold off;

        %ylabel('\ - Low phase bins (VARrel)', 'FontWeight', 'bold', 'FontSize', 14);
        %xlabel('Mu phase effect strength (\rho)', 'FontWeight', 'bold', 'FontSize', 14);
        %set(gca, 'FontSize', 12, 'Box', 'off');

        % Add correlation coefficient and p-value as text on the plot
        text(0.35, 0.99, sprintf('r = %.3f, p = %.3f', r, p_val), ...
            'Units', 'normalized', 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

        %ylim([-40 40])
        %xlim([-0.15, 0.25])
        box off;
        axis square
        %set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
        ax = gca;
        set(ax, 'FontSize', 14);
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.LineWidth = 2;
        fontsize(gcf, scale=1.5)
        fontname('Arial')
    end

    if save_output
        exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR',  sprintf('VAR_%s_%s_%s_%s_%s_individualeffects.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site)), 'Resolution', 600);
    end

end

