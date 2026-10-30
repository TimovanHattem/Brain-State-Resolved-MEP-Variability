%% T2TVAR BRAIN-STATE-RESOLVED VARIABILITY QUANTIFICATION CONTINUOUS
% By Timo van Hattem (timo.van-hattem@uni-tuebingen.de)
% Brain Networks & Plasticity Lab, University of Tuebingen
%
% Description: Quantify brain-state-resolved variability in motor-evoked
% potentials (MEPs) using different EEG features and variability metrics.
% Analysis settings allow different quantification methods to be selected.

%% ----------- CLEAR WORKSPACE -----------

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

var_method          = 'CV'; %SD %CV %MAD %IQR %CQV %TTV_NIU
bin_mean_normalization   = false;
MEP_feature         = 'p2p_microvolts'; %p2p_microvolts %p2p_log_microvolts %p2p_zscored_log
TEP_feature         = 'p2p';        %p2p %AUC %GMFP %LMFP
TEP_metric_type     = 'raw';  %raw %log %z_log
TEP_window          = 'allearly';   %allearly %N15P30 %P30N45 %N45P60

site                = 'Tuebingen'; %Tuebingen %Aalto %combined

prestim_feature     = 'Power'; %Power %Connectivity %Power_relative
bands_of_interest   = {'alpha'}; %{'theta','alpha','beta', 'gamma'};
splits              = 'octiles'; %tails25 %tertiles %median %quintiles %octiles %deciles '4phasebins' '8phasebins'
if strcmpi(splits,'deciles')
    nBins = 10;
elseif strcmpi(splits,'octiles') || strcmpi(splits,'octiles_without_nonresponse_rejection') || strcmpi(splits,'octiles_APB') || strcmpi(splits,'octiles_FDI') || strcmpi(splits,'octiles_with_outlier_rejection') 
    nBins = 8;
elseif strcmpi(splits,'quintiles')
    nBins = 5;
elseif strcmpi(splits,'tails25')
    nBins = 4;
end
nPerState_general   = 50;
nSurrogate          = 50;

plots_show          = 'VARrel'; %VARrel %VAR %MEAN_AMP

if strcmpi(plots_show, 'MEAN_AMP')
    MEP_feature =  'p2p_zscored_log';
    TEP_metric_type = 'z_log';
end

generate_plots      = true;

save_output         = false;

rng(51);

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

% % REVISION BRAIN STIMULATION RECONSTRUCTED CHANNELS POWER
% load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\reftep_C3_affected_subjects_list.mat');
% alldata_MEP_master = alldata_MEP_master(1:29);
% alldata_MEP_master = alldata_MEP_master(~ismember({alldata_MEP_master.subID}, C3_affected_ids));
% all_subjects = {alldata_MEP_master.subID};
% sub_range = 1:17;
% REVISION BRAIN STIMULATION RECONSTRUCTED CHANNELS CONNECTIVITY
% load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\reftep_C3_affected_subjects_list.mat');
% load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\reftep_C4_affected_subjects_list.mat');
% C3C4_affected_ids = unique(horzcat(C3_affected_ids, C4_affected_ids));
% all_subjects = all_subjects(~ismember(all_subjects, C3C4_affected_ids));
% alldata_MEP_master = alldata_MEP_master(1:29);
% alldata_MEP_master = alldata_MEP_master(~ismember({alldata_MEP_master.subID}, C3C4_affected_ids));
% all_subjects = {alldata_MEP_master.subID};
% sub_range = 1:15;

% % REVISION BRAIN STIMULATION WITHOUT NON-RESPONSE REJECTION
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_without_nonresponse_rejection.mat'), 'alldata_MEP_master');
 
% % REVISION BRAIN STIMULATION WITH OUTLIER REJECTION
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_with_outlier_rejection.mat'), 'alldata_MEP_master');

% % REVISION BRAIN STIMULATION APB FDI separately
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_FDI.mat'), 'alldata_MEP_master');

count = 0;
resultsTable = struct([]);
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
        if ~strcmpi(var_method, 'TTV_NIU')
            all_data = alldata_MEP_master(subidx).(MEP_feature);
        elseif strcmpi(var_method, 'TTV_NIU')
            all_data = alldata_MEP_master(subidx).final_meps;
        end
    elseif strcmpi(response, 'TEP')
        all_all_data = alldata_sensor_TEP_master(subidx).all_trials_all_chans;

        % Remove rejected trials based on neuronavigation (+ MEPs, for
        % direct comparison)
        rejected_trials = alldata_MEP_master(subidx).total_trials_rejected;
        all_all_data(:, :, rejected_trials) = [];

        if ~strcmpi(var_method, 'TTV_NIU')
            % --- Compute TEP excitability metric ---
            [raw_tep_metric] = compute_TEP_metrics(all_all_data, 1000, timeVec, TEP_feature, TEP_window, roi_channels, chanlocs);
            % --- Log-transform ---
            log_tep_metric  = log(raw_tep_metric);
            % Z-score within subject
            z_log_tep_metric  = (log_tep_metric  - mean(log_tep_metric))  ./ std(log_tep_metric);

            switch TEP_metric_type
                case 'raw'
                    all_data = raw_tep_metric;
                case 'log'
                    all_data = log_tep_metric;
                case 'z_log'
                    all_data = z_log_tep_metric;
            end
        elseif strcmpi(var_method, 'TTV_NIU')
            % ROI average [nTrials x nTime]
            allLabels = {chanlocs.labels};
            roiIdx = find(ismember(allLabels, roi_channels));
            all_data = squeeze(mean(all_all_data(roiIdx,:,:), 1))';

            % Time window indices
            postWin_TEP = [20 60];
            t_post = find(timeVec >= postWin_TEP(1) & timeVec < postWin_TEP(2));

            % Extract post-stim data
            all_data = all_data(:, t_post);
        end
    end

    % --------- Load entire splits table for feature ---------
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity') || strcmpi(prestim_feature, 'Phase')  || strcmpi(prestim_feature, 'Power_relative') || strcmpi(prestim_feature, 'Power_IAF') || strcmpi(prestim_feature, 'Connectivity_IAF')
        splitsTbl = load_splits_table(DATAIN_SPLIT, sub_num, 1, prestim_feature, splits);
        if isempty(splitsTbl)
            warning('No splits found for %s', subID);
            return
        end
    end

    % --------- Variability quantification ---------
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity') || strcmpi(prestim_feature, 'Power_relative') || strcmpi(prestim_feature, 'Power_IAF') || strcmpi(prestim_feature, 'Connectivity_IAF')

        for fb = 1:numel(bands_of_interest)
            thisBand = bands_of_interest{fb};

            % Load splits table of band
            Tb = splitsTbl(strcmpi(splitsTbl.Freq_band, thisBand), :);
            if isempty(Tb)
                return
            end

            bin_idx = cell(nBins,1);
            for k = 1:nBins
                bin_idx{k} = double(Tb.Trial_index(Tb.Index == k));
            end

            toofew = 0;
            for k = 1:length(bin_idx)
                if length(bin_idx{k}) < nPerState_general
                    toofew = 1;
                end
            end
            if toofew
                disp(['Subject skipped. Not enough trials for ' subID])
                continue
            end

            VAR_bin  = nan(nBins,1);
            MEAN_bin = nan(nBins,1);
            MEAN_all = mean(all_data);

            for k = 1:nBins
                xk = all_data(bin_idx{k});
                MEAN_bin(k) = mean(xk);

                switch upper(var_method)
                    case 'SD'
                        VAR_all = std(all_data);
                        VAR_bin(k) = std(xk);

                        % SD of mean-normalize data
                        if bin_mean_normalization
                            VAR_all = std(all_data / MEAN_all);
                            VAR_bin(k)   = std(xk / MEAN_bin(k));
                        end
                    case 'CV'
                        VAR_all  = std(all_data) / MEAN_all;
                        VAR_bin(k) = std(xk) / mean(xk);
                    case 'MAD'
                        VAR_all  = 1.4826 * mad(all_data, 1);
                        VAR_bin(k) = 1.4826 * mad(xk, 1);
                    case 'IQR'
                        VAR_all  = iqr(all_data);
                        VAR_bin(k) = iqr(xk);
                    case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)
                        VAR_all  = calculate_CQV(all_data);
                        VAR_bin(k) = calculate_CQV(xk);
                    case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                        VAR_all = calculate_TTV_Niu(all_data, nSurrogate);
                        VAR_bin(k) = calculate_TTV_Niu(all_data(xk,:), nSurrogate);
                    otherwise
                        error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method);
                end
            end

            VARrel_bin = (VAR_bin / VAR_all) * 100;

            % Save results to table
            count = count + 1;
            resultsTable(count).subject = subID;
            resultsTable(count).site = site;
            resultsTable(count).type    = horzcat('Hjorth_', prestim_feature);
            resultsTable(count).band    = thisBand;
            resultsTable(count).nBins      = nBins;
            resultsTable(count).splits = splits;
            resultsTable(count).MEAN_all       = MEAN_all;
            resultsTable(count).VAR_all         = VAR_all;
            resultsTable(count).MEAN_bin   = MEAN_bin(:)';      % 1x10 row
            resultsTable(count).VAR_bin    = VAR_bin(:)';       % 1x10 row
            resultsTable(count).VARrel_bin = VARrel_bin(:)';    % 1x10 row
        end
    end
end
T_all = struct2table(resultsTable);

% Save output table
if save_output & ~strcmpi(plots_show, 'MEAN_AMP')
    outfile = sprintf('T_VAR_allblocks_continuous_%s_%s_%s_%s_%s_%s_%s.csv', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(var_method), site, splits, strjoin(bands_of_interest, ''));
    outpath = DATAIN_EP;
    writetable(T_all, fullfile(outpath, outfile));
elseif save_output & strcmpi(plots_show, 'MEAN_AMP')
    outfile = sprintf('T_MEAN_allblocks_continuous_%s_%s_%s_%s_%s_%s_%s.csv', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(var_method), site, splits, strjoin(bands_of_interest, ''));
    outpath = DATAIN_EP;
    writetable(T_all, fullfile(outpath, outfile));
end

%% ----------- PLOTTING -----------
if generate_plots
    if strcmpi(plots_show, 'VARrel')

        % ---- Bands + display names ----
        bandInfo = { ...
            'theta', 'theta';
            'alpha', 'mu';
            'beta',  'beta';
            'gamma', 'gamma'};
        bands     = string(bandInfo(:,1));
        plotNames = string(bandInfo(:,2));

        % ---- Select correct type ----
        thisType = horzcat('Hjorth_', prestim_feature);

        % ---- Define colors per band ----
        colMap = containers.Map;
        colMap('theta') = [120, 165, 170] / 255;
        colMap('alpha') = [150, 125, 160] / 255;
        colMap('beta')  = [190, 145, 95]  / 255;
        colMap('gamma') = [130, 165, 135] / 255;

        % ---- Determine nBins from the first available entry ----
        anySel = find(strcmp({resultsTable.type}, thisType) & isfield(resultsTable,'VARrel_bin'), 1, 'first');
        if isempty(anySel)
            warning('No VARrel_bin found in resultsTable for type %s', thisType);
            return
        end

        x = 1:nBins;

        figure('Color','w','Name',sprintf('%s | %s | %d bins', prestim_feature, upper(var_method), nBins));
        hold on;

        leg = strings(numel(bands),1);

        for fb = 1:numel(bands)
            thisBand = bands(fb);

            sel = strcmp({resultsTable.band}, char(thisBand)) & ...
                strcmp({resultsTable.type}, thisType) & ...
                isfield(resultsTable, 'VARrel_bin');

            if ~any(sel)
                continue
            end

            M = vertcat(resultsTable(sel).VARrel_bin);   % nSub x nBins

            % mean + SEM across subjects (per bin)
            mu = mean(M, 1, 'omitnan');
            nEff = sum(~isnan(M), 1);
            se = std(M, 0, 1, 'omitnan') ./ sqrt(max(nEff,1));

            c = colMap(char(thisBand));

            % % Shaded SEM band
            % xx = [x, fliplr(x)];
            % yy = [mu - se, fliplr(mu + se)];
            % fill(xx, yy, c, 'FaceAlpha', 0.10, 'EdgeColor','none', 'HandleVisibility','off');
            %
            % % Mean line
            % plot(x, mu, '-o', 'LineWidth', 4, 'MarkerSize', 6, 'Color', c, 'MarkerFaceColor', c);

            errorbar(x, mu, se, 'o', ...
                'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            % --- Regression line across bins ---
            pfit = polyfit(x, mu, 1);           % linear fit
            yhat = polyval(pfit, x);            % fitted values

            plot(x, yhat, '-', ...
                'Color', c, ...
                'LineWidth', 4, ...
                'HandleVisibility','off');


            leg(fb) = plotNames(fb);
        end

        yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

        % === Formatting ===
        xlim([0.5 nBins + 0.5]);
        xticks(1:nBins);
        xticklabels(arrayfun(@(x) sprintf('Q%d', x), 1:nBins, 'UniformOutput', false));

        if strcmpi(prestim_feature, 'Power_relative')
            legend({'Normalized mu power', 'Normalized beta power'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 14);
        end


        ylim([85 115])
        box off;

        set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
        ax = gca;
        set(ax, 'FontSize', 14);
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.LineWidth = 2;
        fontsize(gcf, scale=1.5)
        fontname('Arial')

        % xlabel(sprintf('%s bin', prestim_feature), 'FontSize', 16, 'FontWeight','bold');
        % ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 16, 'FontWeight','bold');
        %title(sprintf('%s | %s', prestim_feature, upper(var_method)), 'Interpreter','none');


        % % Adjust label positions
        % xlab = ax.XLabel;
        % xlab.Position(2) = xlab.Position(2) - 0.5;
        % ylab = ax.YLabel;
        % ylab.Position(1) = ylab.Position(1) - 0.1;

        % Save as high-quality PNG (300 dpi)
        if save_output
            if strcmpi(response, 'TEP')
                exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR',  sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s_continuous.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
            elseif strcmpi(response, 'MEP')
                exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR',  sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_continuous.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
            end
        end
   
    elseif strcmpi(plots_show, 'MEAN_AMP')
        % ---- Bands + display names ----
        bandInfo = { ...
            'theta', 'theta';
            'alpha', 'mu-alpha';
            'beta',  'beta';
            'gamma', 'gamma'};
        bands     = string(bandInfo(:,1));
        plotNames = string(bandInfo(:,2));

        % ---- Select correct type ----
        thisType = horzcat('Hjorth_', prestim_feature);

        % ---- Define colors per band ----
        colMap = containers.Map;
        colMap('theta') = [120, 165, 170] / 255;
        colMap('alpha') = [150, 125, 160] / 255;
        colMap('beta')  = [190, 145, 95]  / 255;
        colMap('gamma') = [130, 165, 135] / 255;

        % ---- Determine nBins from the first available entry ----
        anySel = find(strcmp({resultsTable.type}, thisType) & isfield(resultsTable,'MEAN_bin'), 1, 'first');
        if isempty(anySel)
            warning('No VARrel_bin found in resultsTable for type %s', thisType);
            return
        end

        x = 1:nBins;

        figure('Color','w','Name',sprintf('%s | %s | %d bins', prestim_feature, upper(var_method), nBins));
        hold on;

        leg = strings(numel(bands),1);

        for fb = 1:numel(bands)
            thisBand = bands(fb);

            sel = strcmp({resultsTable.band}, char(thisBand)) & ...
                strcmp({resultsTable.type}, thisType) & ...
                isfield(resultsTable, 'VARrel_bin');

            if ~any(sel)
                continue
            end

            M = vertcat(resultsTable(sel).MEAN_bin);   % nSub x nBins

            % mean + SEM across subjects (per bin)
            mu = mean(M, 1, 'omitnan');
            nEff = sum(~isnan(M), 1);
            se = std(M, 0, 1, 'omitnan') ./ sqrt(max(nEff,1));

            c = colMap(char(thisBand));

            % % Shaded SEM band
            % xx = [x, fliplr(x)];
            % yy = [mu - se, fliplr(mu + se)];
            % fill(xx, yy, c, 'FaceAlpha', 0.10, 'EdgeColor','none', 'HandleVisibility','off');
            %
            % % Mean line
            % plot(x, mu, '-o', 'LineWidth', 4, 'MarkerSize', 6, 'Color', c, 'MarkerFaceColor', c);

            errorbar(x, mu, se, 'o', ...
                'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            % --- Regression line across bins ---
            pfit = polyfit(x, mu, 1);           % linear fit
            yhat = polyval(pfit, x);            % fitted values

            plot(x, yhat, '-', ...
                'Color', c, ...
                'LineWidth', 4, ...
                'HandleVisibility','off');


            leg(fb) = plotNames(fb);
        end

        yline(0, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

        % === Formatting ===
        xlim([0.5 nBins + 0.5]);
        xticks(1:nBins);
        xticklabels(arrayfun(@(x) sprintf('Q%d', x), 1:nBins, 'UniformOutput', false));

        if strcmpi(prestim_feature, 'Power_relative')
            legend({'Normalized mu power', 'Beta (n.s.)'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 14);
        end


        ylim([-0.2 0.2])
        box off;

        set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
        ax = gca;
        set(ax, 'FontSize', 14);
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.LineWidth = 2;
        fontsize(gcf, scale=1.5)
        fontname('Arial')

        % xlabel(sprintf('%s bin', prestim_feature), 'FontSize', 16, 'FontWeight','bold');
        % ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 16, 'FontWeight','bold');
        %title(sprintf('%s | %s', prestim_feature, upper(var_method)), 'Interpreter','none');


        % % Adjust label positions
        % xlab = ax.XLabel;
        % xlab.Position(2) = xlab.Position(2) - 0.5;
        % ylab = ax.YLabel;
        % ylab.Position(1) = ylab.Position(1) - 0.1;

        % Save as high-quality PNG (300 dpi)
        if save_output
            if strcmpi(response, 'TEP')
                exportgraphics(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_AMP_%s_%s_%s_%s_%s_%s_%s_%s_%s_continuous.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
            elseif strcmpi(response, 'MEP')
                exportgraphics(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP',  sprintf('MEAN_AMP_%s_%s_%s_%s_%s_%s_%s_continuous.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
            end
        end
    end


    %% ----------- SUBJECT-LEVEL REGRESSION SPAGHETTI PLOT -----------
    x = (1:nBins)';

    figure('Color','w','Name',sprintf('%s | subject regression lines | %d bins', prestim_feature, nBins));
    hold on;

    for fb = 1:numel(bands)
        thisBand = bands(fb);

        sel = strcmp({resultsTable.band}, char(thisBand)) & ...
              strcmp({resultsTable.type}, thisType) & ...
              isfield(resultsTable, 'VARrel_bin');

        if ~any(sel), continue; end

        % Subject x bins matrix
        M = vertcat(resultsTable(sel).VARrel_bin);   % nSub x nBins

        c = colMap(char(thisBand));

        % --- Fit and plot each subject line ---
        nSub = size(M,1);
        subjSlopes = nan(nSub,1);
        subjInts   = nan(nSub,1);
        Yhat       = nan(nSub, nBins);

        nSubjects = 29;
        colors = parula(nSubjects);  % or lines(nSubjects), jet(nSubjects), turbo(nSubjects), etc.


        for s = 1:nSub
            y = M(s,:)';
            ok = ~isnan(y);

            if sum(ok) < 3
                continue
            end

            p = polyfit(x(ok), y(ok), 1);
            subjSlopes(s) = p(1);
            subjInts(s)   = p(2);

            yhat = polyval(p, x);
            Yhat(s,:) = yhat;

            % Faint subject regression line
            plot(x, yhat, '-', ...
                'Color', colors(s,:), ...
                'LineWidth', 3, ...
                'HandleVisibility','off');
        end

        % --- Group-average fitted line (mean of subject fits) ---
        mu_hat = mean(Yhat, 1, 'omitnan');
        plot(x, mu_hat, '-', ...
            'Color', 'k', ...
            'LineWidth', 4, ...
            'DisplayName', plotNames(fb));

        % Optional: print slope summary in console
        ms = mean(subjSlopes, 'omitnan');
        fprintf('%s: mean(subject slopes)=%.4f (N=%d)\n', thisBand, ms, sum(~isnan(subjSlopes)));
    end

    % === Formatting ===
    yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

    % === Formatting ===
    xlim([0.5 nBins + 0.5]);
    xticks(1:nBins);
    xticklabels(arrayfun(@(x) sprintf('Q%d', x), 1:nBins, 'UniformOutput', false));

    if strcmpi(prestim_feature, 'Power_relative')
        legend({'Normalized mu power', 'Normalized beta power'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 14);
    end


    ylim([80 120])
    box off;

    set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
    ax = gca;
    set(ax, 'FontSize', 14);
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.LineWidth = 2;
    fontsize(gcf, scale=1.5)
    fontname('Arial')


    % Save as high-quality PNG (300 dpi)
    if save_output
        exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_continuous_subjects.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
    end
    
end

%% TREND STATISTICS

if strcmpi(plots_show, 'VARrel')
    table_name_of_interest = 'VARrel_bin'
elseif strcmpi(plots_show, 'MEAN_AMP')
    table_name_of_interest = 'MEAN_bin'
end

binIndex = (1:nBins)';

thisType = ['Hjorth_' prestim_feature];

% -------- helper: Benjamini-Hochberg FDR --------
bh_fdr = @(p) local_bh_fdr(p, 0.05);  %#ok<NASGU> (defined at end)

% -------- collect slopes per band --------
TrendStats = table('Size',[numel(bands) 9], ...
    'VariableTypes', {'string','double','double', 'double','double','double','double','double','double'}, ...
    'VariableNames', {'Band','Nsub','MeanSlope', 'SDSlope', 'SEMSlope','df','t','p_twoSided','p_oneSided'});

all_p = nan(numel(bands),1);

for b = 1:numel(bands)
    band = bands{b};

    % Select entries for this band/type and quintiles (nBins==5)
    sel = strcmpi({resultsTable.band}, band) & strcmpi({resultsTable.type}, thisType) & isfield(resultsTable, table_name_of_interest);

    if ~any(sel)
        warning('No entries found for band=%s, type=%s', band, thisType);
        TrendStats.Band(b) = string(band);
        continue
    end

    R = resultsTable(sel);

    % Ensure quintiles (1x5). If mixed nBins exist, filter
    nBinsVec = arrayfun(@(x) numel(x.(table_name_of_interest)), R);
    R = R(nBinsVec == nBins);
    if isempty(R)
        warning('No quintile (5-bin) entries for band=%s, type=%s', band, thisType);
        TrendStats.Band(b) = string(band);
        continue
    end

    % Group by subject: compute one slope per subject
    subj = string({R.subject});
    uSubj = unique(subj, 'stable');
    slopes = nan(numel(uSubj),1);

    for s = 1:numel(uSubj)
        rs = R(subj == uSubj(s));

        % If multiple rows per subject for this band/type, average VARrel_bin first
        M = vertcat(rs.(table_name_of_interest));          % k x 5
        y = mean(M, 1, 'omitnan')';          % 5 x 1

        % Skip if too many NaNs
        ok = ~isnan(y);
        if sum(ok) < 3
            slopes(s) = NaN;
            continue
        end

        % Fit slope y ~ a + b*x using polyfit (or regress)
        pfit = polyfit(binIndex(ok), y(ok), 1);
        slopes(s) = pfit(1);  % b (expected > 0 if variability increases as bin goes from high->low)
    end

    slopes = slopes(~isnan(slopes));
    N = numel(slopes);

    if N < 3
        warning('Not enough subjects with valid slopes for band=%s (N=%d)', band, N);
        TrendStats.Band(b) = string(band);
        TrendStats.Nsub(b) = N;
        continue
    end

    % One-sample t-test for slope > 0
    [~, p_two, ci, stats] = ttest(slopes, 0, 'Tail','left');
    p_one = p_two;

    % Two-sided p for completeness
    [~, p2, ci, stats] = ttest(slopes, 0, 'Tail','both');
  
    % Get descriptive values
    mu = mean(slopes);
    sd = std(slopes);
    sem = std(slopes) / sqrt(N);

    % Create output
    TrendStats.Band(b)      = string(band);
    TrendStats.Nsub(b)      = N;
    TrendStats.MeanSlope(b) = mu;
    TrendStats.SDSlope(b)   = sd;
    TrendStats.SEMSlope(b)  = sem;
    TrendStats.df(b)        = stats.df;
    TrendStats.t(b)         = stats.tstat;
    TrendStats.p_twoSided(b)= p2;
    TrendStats.p_oneSided(b)= p_one;

    all_p(b) = p2;
end

% -------- FDR across bands (one-sided p) --------
[qvals, sig] = local_bh_fdr(all_p, 0.05);
TrendStats.FDR_q = qvals;
TrendStats.FDR_sig = sig;

disp('--- Trend (slope>0) stats per band ---');
disp(TrendStats);


%% Helper functions

function [qvals, sig] = local_bh_fdr(pvals, alpha)
if nargin < 2, alpha = 0.05; end
p = pvals(:);
m = numel(p);
qvals = nan(size(p));
sig = false(size(p));

ok = ~isnan(p);
p_ok = p(ok);
[ps, idx] = sort(p_ok, 'ascend');
m_ok = numel(ps);

% BH critical values
crit = (1:m_ok)'/m_ok * alpha;
k = find(ps <= crit, 1, 'last');
if ~isempty(k)
    sig_ok = false(m_ok,1);
    sig_ok(1:k) = true;
else
    sig_ok = false(m_ok,1);
end

% q-values (Benjamini-Hochberg adjusted p)
q = ps .* m_ok ./ (1:m_ok)';
% enforce monotonicity
for i = m_ok-1:-1:1
    q(i) = min(q(i), q(i+1));
end

% map back
qvals_ok = nan(m_ok,1);
qvals_ok(idx) = q;
sig_ok2 = false(m_ok,1);
sig_ok2(idx) = sig_ok;

qvals(ok) = qvals_ok;
sig(ok) = sig_ok2;
end
