%% T2TVAR BRAIN-STATE-RESOLVED VARIABILITY QUANTIFICATION SLIDING WINDOW
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

prestim_feature     = 'Power'; %Power %Connectivity
bands_of_interest   = {'alpha'}; %{'theta','alpha','beta', 'gamma'};
splits              = 'tails25'; %tails25 %tertiles %median %quintiles %octiles %deciles '4phasebins' '8phasebins'
nPhaseBins          = 4;
if strcmpi(prestim_feature, 'Phase') && nPhaseBins == 4
    splits          = '4phasebins';
elseif strcmpi(prestim_feature, 'Phase') && nPhaseBins == 8
    splits          = '8phasebins';
end
nPerState_general   = 100;
nPerState_phase     = 100;
nPerState_combos    = 25;
nSurrogate          = 50;

windowSize       = 300;    % number of trials per window
stepSize         = 50;     % how many trials to move the window each iteration

plots_show          = 'VARrel'; %VARrel %VAR

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
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity') || strcmpi(prestim_feature, 'Phase')
        splitsTbl = load_splits_table(DATAIN_SPLIT, sub_num, 1, prestim_feature, splits);
        if isempty(splitsTbl)
            warning('No splits found for %s', subID);
            return
        end
    end

    trialIdx = 1:length(all_data);
    nWindows = floor((length(trialIdx) - windowSize) / stepSize) + 1;

    % --------- Variability quantification ---------
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')

        for fb = 1:numel(bands_of_interest)
            thisBand = bands_of_interest{fb};

            % Load splits table of band
            Tb = splitsTbl(strcmpi(splitsTbl.Freq_band, thisBand), :);
            if isempty(Tb)
                return
            end

            % Trial-wise prestim values
            PrestimVals = nan(size(all_data));
            if strcmpi(prestim_feature, 'Power')
                PrestimVals(Tb.Trial_index) = Tb.Hjorth_C3_Power;
            elseif strcmpi(prestim_feature, 'Connectivity')
                PrestimVals(Tb.Trial_index) = Tb.Hjorth_C3C4_Connectivity;
            end

            MEAN_prestimVal_high_time = nan(1, nWindows);
            MEAN_prestimVal_low_time = nan(1, nWindows);
            winIdx_time = nan(nWindows, windowSize);
            MEAN_all_time  = nan(1, nWindows);
            MEAN_high_time   = nan(1, nWindows);
            MEAN_low_time   = nan(1, nWindows);
            MEAN_rand_time = nan(1, nWindows);
            VAR_all_time  = nan(1, nWindows);
            VAR_high_time   = nan(1, nWindows);
            VAR_low_time   = nan(1, nWindows);
            VAR_rand_time = nan(1, nWindows);
            VARrel_high_time = nan(1, nWindows);
            VARrel_low_time  = nan(1, nWindows);
            VARrel_rand_time = nan(1, nWindows);

            for w = 1:nWindows
                idxStart = (w - 1) * stepSize + 1;
                idxEnd   = idxStart + windowSize - 1;
                winIdx   = idxStart:idxEnd;

                ep_win = all_data(winIdx);
                prestimval_win = PrestimVals(winIdx);

                validIdx = ~isnan(prestimval_win) & ~isnan(ep_win);
                if sum(validIdx) < windowSize * 0.5
                    return
                end

                % Sort trials by power, split into quartiles
                [~, sortIdx] = sort(prestimval_win);
                nQ = floor(length(sortIdx)/4);
                qIdx = reshape(sortIdx(1:4*nQ), nQ, 4);

                high_idx = qIdx(:,4);
                low_idx = qIdx(:,1);
                rand_idx = randperm(length(ep_win), numel(low_idx));

                % Mean across trials
                MEAN_all  = mean(ep_win);
                POWERMEAN_high   = mean(ep_win(high_idx));
                POWERMEAN_low   = mean(ep_win(low_idx));
                MEAN_rand = mean(ep_win(rand_idx));

                % Compute variability per state
                switch upper(var_method)
                    case 'SD' % Standard deviation (assumes log-microvolts)
                        VAR_all = std(ep_win);
                        VAR_high = std(ep_win(high_idx));
                        VAR_low = std(ep_win(low_idx));
                        VAR_rand = std(ep_win(rand_idx));

                        % SD of mean-normalize data
                        if bin_mean_normalization
                            VAR_all = std((ep_win / MEAN_all));
                            VAR_high   = std((ep_win(high_idx) / POWERMEAN_high));
                            VAR_low   = std((ep_win(low_idx) / POWERMEAN_low));
                            VAR_rand = std((ep_win(rand_idx) / MEAN_rand));
                        end
                    case 'CV' % Coefficient of variation (assumes raw microvolts)
                        VAR_all  = std(ep_win) / mean(ep_win);
                        VAR_high   = std(ep_win(high_idx)) / mean(ep_win(high_idx));
                        VAR_low   = std(ep_win(low_idx)) / mean(ep_win(low_idx));
                        VAR_rand = std(ep_win(rand_idx)) / mean(ep_win(rand_idx));
                    case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                        VAR_all  = 1.4826 * mad(ep_win, 1);
                        VAR_high   = 1.4826 * mad(ep_win(high_idx), 1);
                        VAR_low   = 1.4826 * mad(ep_win(low_idx), 1);
                        VAR_rand = 1.4826 * mad(ep_win(rand_idx), 1);
                    case 'IQR' % Interquartile range (assumes log-microvolts)
                        VAR_all  = iqr(ep_win);
                        VAR_high   = iqr(ep_win(high_idx));
                        VAR_low   = iqr(ep_win(low_idx));
                        VAR_rand = iqr(ep_win(rand_idx));
                    case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)
                        VAR_all  = calculate_CQV(ep_win);
                        VAR_high = calculate_CQV(ep_win(high_idx));
                        VAR_low  = calculate_CQV(ep_win(low_idx));
                        VAR_rand = calculate_CQV(ep_win(rand_idx));
                    case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                        VAR_all = calculate_TTV_Niu(ep_win, nSurrogate);
                        VAR_high = calculate_TTV_Niu(ep_win(high_idx,:), nSurrogate);
                        VAR_low = calculate_TTV_Niu(ep_win(low_idx,:), nSurrogate);
                        VAR_rand = calculate_TTV_Niu(ep_win(rand_idx,:), nSurrogate);
                    otherwise
                        error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
                end

                % relative spreads
                VARrel_high = (VAR_high   / VAR_all) * 100;
                VARrel_low  = (VAR_low   / VAR_all) * 100;
                VARrel_rand = (VAR_rand / VAR_all) * 100;

                % Save time metrics
                MEAN_prestimVal_high_time(w) = mean(prestimval_win(high_idx));
                MEAN_prestimVal_low_time(w) = mean(prestimval_win(low_idx));
                winIdx_time(w,:) = winIdx;
                MEAN_all_time(w)  = MEAN_all;
                MEAN_high_time(w)   = POWERMEAN_high;
                MEAN_low_time(w)   = POWERMEAN_low;
                MEAN_rand_time(w) = MEAN_rand;
                VAR_all_time(w)  = VAR_all;
                VAR_high_time(w)   = VAR_high;
                VAR_low_time(w)   = VAR_low;
                VAR_rand_time(w) = VAR_rand;
                VARrel_high_time(w) = VARrel_high;
                VARrel_low_time(w)  = VARrel_low;
                VARrel_rand_time(w) = VARrel_rand;
            end

            % Store in results structure
            count = count + 1;
            resultsTable(count).subject = subID;
            resultsTable(count).type    = horzcat('Hjorth_', prestim_feature);
            resultsTable(count).band    = thisBand;
            resultsTable(count).winCenters = winIdx_time;
            resultsTable(count).MEAN_prestimVal_high_time = MEAN_prestimVal_high_time;
            resultsTable(count).MEAN_prestimVal_low_time = MEAN_prestimVal_low_time;
            resultsTable(count).MEAN_all_time       = MEAN_all_time;
            resultsTable(count).MEAN_high_time       = MEAN_high_time;
            resultsTable(count).MEAN_low_time       = MEAN_low_time;
            resultsTable(count).MEAN_rand_time       = MEAN_rand_time;
            resultsTable(count).VAR_all_time       = VAR_all_time;
            resultsTable(count).VAR_high_time       = VAR_high_time;
            resultsTable(count).VAR_low_time       = VAR_low_time;
            resultsTable(count).VAR_rand_time       = VAR_rand_time;
            resultsTable(count).VARrel_high_time       = VARrel_high_time;
            resultsTable(count).VARrel_low_time        = VARrel_low_time;
            resultsTable(count).VARrel_rand_time       = VARrel_rand_time;

            % ===== Normalize & interpolate each subject individually (matrix version) =====
            commonN = 100;  % number of normalized time points
            commonAxis = linspace(0, 1, commonN);  % normalized session timeline

            for i = 1:numel(resultsTable)
                % --- Compute window centers from window index matrix ---
                if isempty(resultsTable(i).winCenters) || all(isnan(resultsTable(i).winCenters(:)))
                    continue
                end
                winCenters_perWin = mean(resultsTable(i).winCenters, 2);  % nWindows × 1
                xNorm = (winCenters_perWin - min(winCenters_perWin)) / (max(winCenters_perWin) - min(winCenters_perWin));

                % --- Interpolate all time-resolved variables to the common axis ---
                vars = {'MEAN_prestimVal_high_time', 'MEAN_prestimVal_low_time', 'MEAN_all_time','MEAN_high_time','MEAN_low_time','MEAN_rand_time', 'VAR_all_time','VAR_high_time','VAR_low_time','VAR_rand_time','VARrel_high_time','VARrel_low_time','VARrel_rand_time'};
                for v = 1:numel(vars)
                    y = resultsTable(i).(vars{v});
                    if isempty(y) || all(isnan(y))
                        resultsTable(i).([vars{v} '_interp']) = nan(1, commonN);
                    else
                        resultsTable(i).([vars{v} '_interp']) = interp1(xNorm, y, commonAxis, 'linear', 'extrap');
                    end
                end

                % Store normalized axis (0–1 progression)
                resultsTable(i).xNorm = commonAxis;
            end
        end
    end
end

% Convert to tidy table
T_all = struct2table(resultsTable);
if save_output
    if strcmpi(response, 'TEP')
        outfile = sprintf('T_VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s_slidingwindow.csv', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits,  strjoin(bands_of_interest, ''));
        outpath = DATAIN_EP;
    elseif strcmpi(response, 'MEP')
        outfile = sprintf('T_VAR_%s_%s_%s_%s_%s_%s_%s_slidingwindow.csv',  upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(var_method), site, splits,  strjoin(bands_of_interest, ''));
        outpath = DATAIN_EP;
    end
    writetable(T_all, fullfile(outpath, outfile));
end

%% ====== Group-level sliding-window analysis ======
bands = bands_of_interest;
nBands = numel(bands);

% Define frequency bands in desired order + pretty names
allBandInfo = { ...
    'theta', 'theta';
    'alpha', 'mu-alpha';
    'beta', 'beta';
    'gamma', 'gamma'};

[~, idx] = ismember(bands_of_interest, allBandInfo(:,1));
idx = idx(idx > 0); % Remove zeros (bands not found)
bandInfo = allBandInfo(idx, :);
bands     = string(bandInfo(:,1));  % raw band names in desired order
plotNames = string(bandInfo(:,2));  % pretty labels

col_low = [59, 82, 139] / 255;    % Medium blue
col_high = [253, 70, 101] / 255;    % Vibrant red
zero_col = [0.3 0.3 0.3];          % gray line for zero

figure('Color','w','Name','Sliding-window');
%tiledlayout(1, nBands, "TileSpacing","loose", "Padding","compact")
xCommon = linspace(0,1,100);

for b = 1:nBands
    thisBand = bands{b};
    idx = strcmp({resultsTable.band}, thisBand);
    dat = resultsTable(idx);

    if isempty(dat)
        warning('No data for band: %s', thisBand);
        continue
    end

    % Normalized axis
    if iscell(dat(1).xNorm)
        x = dat(1).xNorm{1};
    else
        x = dat(1).xNorm;
    end

    nexttile;

    % ===== PLOT VARIABILITY (LEFT AXIS) =====
    VAR_high  = cat(1, dat(~cellfun(@isempty,{dat.([(plots_show), '_high_time_interp'])})).([(plots_show), '_high_time_interp']));
    VAR_low  = cat(1, dat(~cellfun(@isempty,{dat.([(plots_show), '_low_time_interp'])})).([(plots_show), '_low_time_interp']));
    VAR_high_mean = mean(VAR_high, 'omitnan');
    VAR_high_sem = std(VAR_high, 'omitnan') ./ sqrt(size(VAR_high,1));
    VAR_low_mean = mean(VAR_low, 'omitnan');
    VAR_low_sem = std(VAR_low, 'omitnan') ./ sqrt(size(VAR_low,1));  % Fixed typo: was VAR_high_low

    plot(xCommon, VAR_high_mean, '-', 'Color', col_high, 'LineWidth', 4);
    hold on;
    plot(xCommon, VAR_low_mean, '-', 'Color', col_low, 'LineWidth', 4);
    fill([xCommon fliplr(xCommon)], [VAR_high_mean-VAR_high_sem fliplr(VAR_high_mean+VAR_high_sem)], col_high, 'FaceAlpha',0.15, 'EdgeColor','none');
    fill([xCommon fliplr(xCommon)], [VAR_low_mean-VAR_low_sem fliplr(VAR_low_mean+VAR_low_sem)], col_low,  'FaceAlpha',0.15, 'EdgeColor','none');

    % % ===== PLOT MEAN MEP (RIGHT AXIS) =====
    % yyaxis right
    % MEPMEAN_high = cat(1, dat(~cellfun(@isempty,{dat.MEAN_high_time_interp})).MEAN_high_time_interp);
    % MEPMEAN_low = cat(1, dat(~cellfun(@isempty,{dat.MEAN_low_time_interp})).MEAN_low_time_interp);
    % MEPMEAN_high_mean = mean(MEPMEAN_high, 'omitnan');
    % MEPMEAN_low_mean = mean(MEPMEAN_low, 'omitnan');
    %
    % plot(xCommon, MEPMEAN_high_mean, '--', 'Color', col_high, 'LineWidth', 1.6, 'DisplayName', 'High/Random (Mean)');
    % plot(xCommon, MEPMEAN_low_mean, '--', 'Color', col_low, 'LineWidth', 1.6, 'DisplayName', 'Low/Random (Mean)');
    %
    % % ===== PLOT PRESTIM POWER (RESCALED TO RIGHT AXIS) =====
    % POWERMEAN_high = cat(1, dat(~cellfun(@isempty,{dat.MEAN_prestimVal_high_time_interp})).MEAN_prestimVal_high_time_interp);
    % POWERMEAN_low = cat(1, dat(~cellfun(@isempty,{dat.MEAN_prestimVal_low_time_interp})).MEAN_prestimVal_low_time_interp);
    % POWERMEAN_high_mean = mean(POWERMEAN_high, 'omitnan');
    % POWERMEAN_low_mean = mean(POWERMEAN_low, 'omitnan');
    %
    % % Rescale prestim power to match MEP amplitude range
    % current_ylim = ylim;
    % low_all = min([POWERMEAN_high_mean POWERMEAN_low_mean]);
    % high_all = max([POWERMEAN_high_mean POWERMEAN_low_mean]);
    % POWERMEAN_high_mean_rescaled = rescale(POWERMEAN_high_mean, current_ylim(1), current_ylim(2), 'InputMin', low_all, 'InputMax', high_all);
    % POWERMEAN_low_mean_rescaled = rescale(POWERMEAN_low_mean, current_ylim(1), current_ylim(2), 'InputMin', low_all, 'InputMax', high_all);
    %
    % plot(xCommon, POWERMEAN_high_mean_rescaled, ':', 'Color', col_high, 'LineWidth', 1.3, 'DisplayName', 'Prestim Power - High');
    % plot(xCommon, POWERMEAN_low_mean_rescaled, ':', 'Color', col_low, 'LineWidth', 1.3, 'DisplayName', 'Prestim Power - Low');
    %
    % ylabel('Mean MEP amplitude (a.u.)');
    % set(gca,'YColor','k');
    % if strcmpi(response, 'MEP')
    %     ylim([800 1500]);
    % elseif strcmpi(response, 'TEP')
    %     ylim([1, 3]);
    % end

    yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');
    
    legend({'High', 'Low'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 14);
    
    % --- Style ---
    %title(strrep(plotNames(b),'_','-'), 'Interpreter','none', 'FontWeight','bold');
    xticks([0 0.25 0.5 0.75 1]);
    
    %set(gca,'YColor','k');
    ylim([85 115])
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 14;
    ax.TickDir = 'out';
    ax.LineWidth = 2;


    % xlabel('Normalized session time', 'FontWeight', 'bold', 'FontSize', 12);
    % ylabel(sprintf('Normalized within-state %s', var_method), 'FontWeight', 'bold', 'FontSize', 12);
end
set(gcf, 'Units', 'normalized', 'Position', [0.666796875,0.152083333333333,0.238203125,0.772916666666667]);
    
fontsize(gcf,scale=1.5)
fontname('Arial')

% Save as high-quality PNG (300 dpi)
if save_output
    % fig = gcf;
    % fig.Units = 'pixels';
    % fig.Position = [1442 508 927 629];  % [left bottom width height]

    if strcmpi(response, 'TEP')
        print(gcf, fullfile(DATAOUT_VIS, sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s_slidingwindow.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits,  strjoin(bands_of_interest, ''))), '-dpng', '-r300');
    elseif strcmpi(response, 'MEP')
        exportgraphics(gcf, fullfile(DATAOUT_VIS, sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_slidingwindow.png', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(var_method), site, splits,  strjoin(bands_of_interest, ''))), 'Resolution', 600);
    end
end

%% REVISION CUMULATIVE EFFECT MEP AMPLITUDE INLCUDED

%% ====== Group-level sliding-window analysis ======
bands = bands_of_interest;
nBands = numel(bands);

% Define frequency bands in desired order + pretty names
allBandInfo = { ...
    'theta', 'theta';
    'alpha', 'mu-alpha';
    'beta', 'beta';
    'gamma', 'gamma'};

[~, idx] = ismember(bands_of_interest, allBandInfo(:,1));
idx = idx(idx > 0); % Remove zeros (bands not found)
bandInfo = allBandInfo(idx, :);
bands     = string(bandInfo(:,1));  % raw band names in desired order
plotNames = string(bandInfo(:,2));  % pretty labels

col_low = [59, 82, 139] / 255;    % Medium blue
col_high = [253, 70, 101] / 255;    % Vibrant red
zero_col = [0.3 0.3 0.3];          % gray line for zero

figure('Color','w','Name','Sliding-window');
%tiledlayout(1, nBands, "TileSpacing","loose", "Padding","compact")
xCommon = linspace(0,1,100);

for b = 1:nBands
    thisBand = bands{b};
    idx = strcmp({resultsTable.band}, thisBand);
    dat = resultsTable(idx);

    if isempty(dat)
        warning('No data for band: %s', thisBand);
        continue
    end

    % Normalized axis
    if iscell(dat(1).xNorm)
        x = dat(1).xNorm{1};
    else
        x = dat(1).xNorm;
    end

    nexttile;

    % ===== PLOT VARIABILITY (LEFT AXIS) =====
    VAR_high  = cat(1, dat(~cellfun(@isempty,{dat.([(plots_show), '_high_time_interp'])})).([(plots_show), '_high_time_interp']));
    VAR_low  = cat(1, dat(~cellfun(@isempty,{dat.([(plots_show), '_low_time_interp'])})).([(plots_show), '_low_time_interp']));
    VAR_high_mean = mean(VAR_high, 'omitnan');
    VAR_high_sem = std(VAR_high, 'omitnan') ./ sqrt(size(VAR_high,1));
    VAR_low_mean = mean(VAR_low, 'omitnan');
    VAR_low_sem = std(VAR_low, 'omitnan') ./ sqrt(size(VAR_low,1));  % Fixed typo: was VAR_high_low

    plot(xCommon, VAR_high_mean, '-', 'Color', col_high, 'LineWidth', 4);
    hold on;
    plot(xCommon, VAR_low_mean, '-', 'Color', col_low, 'LineWidth', 4);
    %fill([xCommon fliplr(xCommon)], [VAR_high_mean-VAR_high_sem fliplr(VAR_high_mean+VAR_high_sem)], col_high, 'FaceAlpha',0.15, 'EdgeColor','none');
    %fill([xCommon fliplr(xCommon)], [VAR_low_mean-VAR_low_sem fliplr(VAR_low_mean+VAR_low_sem)], col_low,  'FaceAlpha',0.15, 'EdgeColor','none');
    ylim([85 115])

    % ===== PLOT MEAN MEP (RIGHT AXIS) =====
    yyaxis right
    MEPMEAN_high = cat(1, dat(~cellfun(@isempty,{dat.MEAN_high_time_interp})).MEAN_high_time_interp);
    MEPMEAN_low = cat(1, dat(~cellfun(@isempty,{dat.MEAN_low_time_interp})).MEAN_low_time_interp);
    MEPMEAN_high_mean = mean(MEPMEAN_high, 'omitnan');
    MEPMEAN_low_mean = mean(MEPMEAN_low, 'omitnan');

    plot(xCommon, MEPMEAN_high_mean, '--', 'Color', col_high, 'LineWidth', 4, 'DisplayName', 'High/Random (Mean)');
    plot(xCommon, MEPMEAN_low_mean, '--', 'Color', col_low, 'LineWidth', 4, 'DisplayName', 'Low/Random (Mean)');

    % ===== PLOT PRESTIM POWER (RESCALED TO RIGHT AXIS) =====
    POWERMEAN_high = cat(1, dat(~cellfun(@isempty,{dat.MEAN_prestimVal_high_time_interp})).MEAN_prestimVal_high_time_interp);
    POWERMEAN_low = cat(1, dat(~cellfun(@isempty,{dat.MEAN_prestimVal_low_time_interp})).MEAN_prestimVal_low_time_interp);
    POWERMEAN_high_mean = mean(POWERMEAN_high, 'omitnan');
    POWERMEAN_low_mean = mean(POWERMEAN_low, 'omitnan');

    % Rescale prestim power to match MEP amplitude range
    current_ylim = ylim;
    low_all = min([POWERMEAN_high_mean POWERMEAN_low_mean]);
    high_all = max([POWERMEAN_high_mean POWERMEAN_low_mean]);
    POWERMEAN_high_mean_rescaled = rescale(POWERMEAN_high_mean, current_ylim(1), current_ylim(2), 'InputMin', low_all, 'InputMax', high_all);
    POWERMEAN_low_mean_rescaled = rescale(POWERMEAN_low_mean, current_ylim(1), current_ylim(2), 'InputMin', low_all, 'InputMax', high_all);

    plot(xCommon, POWERMEAN_high_mean_rescaled, ':', 'Color', col_high, 'LineWidth', 4, 'DisplayName', 'Prestim Power - High');
    plot(xCommon, POWERMEAN_low_mean_rescaled, ':', 'Color', col_low, 'LineWidth', 4, 'DisplayName', 'Prestim Power - Low');

    set(gca,'YColor','k');
    if strcmpi(response, 'MEP')
        ylim([800 1500]);
    elseif strcmpi(response, 'TEP')
        ylim([1, 3]);
    end

    yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');
    
    lgd = legend({'MEP variability - High', 'MEP variability - Low', 'MEP amplitude - High', ...
        'MEP amplitude - Low', 'Pre-stimulus M1-M1 iPLV - High', 'Pre-stimulus M1-M1 iPLV - Low'}, ...
        'Box', 'off', 'FontSize', 12);
    lgd.ItemTokenSize = [50, 18];  % increase 60 for even longer lines
    lgd.Orientation = 'vertical';
    lgd.Location = 'southoutside';
    %--- Style ---
    title(strrep(plotNames(b),'_','-'), 'Interpreter','none', 'FontWeight','bold');
    xticks([0 0.25 0.5 0.75 1]);
    
    %set(gca,'YColor','k');
    
    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 14;
    ax.TickDir = 'out';
    ax.LineWidth = 2;

    title(prestim_feature)

    xlabel('Normalized session time', 'FontWeight', 'bold', 'FontSize', 14);
    yyaxis left
    ylabel(sprintf('Normalized MEP variability (CV, %%)', var_method), 'FontWeight', 'bold', 'FontSize', 14);
    yyaxis right
    ylabel('Mean MEP amplitude (\muV)', 'FontSize', 14,  'FontWeight', 'bold');    
end
set(gcf, 'Units', 'normalized', 'Position', [0.57578125,0.152083333333333,0.32921875,0.772916666666667]);
    
fontsize(gcf,scale=1.6)
fontname('Arial')

if save_output
    exportgraphics(gcf, fullfile(DATAOUT_VIS, sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_slidingwindow_withMEPamplitude.png', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(var_method), site, splits,  strjoin(bands_of_interest, ''))), 'Resolution', 600);
end

