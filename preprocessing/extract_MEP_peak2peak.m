%% T2TVAR EXTRACT MEP AMPLITUDES
% By Timo van Hattem (timo.van-hattem@uni-tuebingen.de)
% Brain Networks & Plasticity Lab, University of Tuebingen
%
% Extract peak-to-peak motor-evoked potential (MEP) amplitudes from EMG
% recordings.

%% ----------- CLEAR WORKSPACE -----------

clear;          % Clear workspace variables
close all;      % Close all figure windows
clc;            % Clear command window

% Get full path of current script and move to its directory
filePath = matlab.desktop.editor.getActiveFilename;
mainpath = fileparts(filePath);   % Adjust to project root folder structure
cd(mainpath);

% Set paths to relevant toolboxes
TOOLBOXPATH = 'D:\2025-08 T2TVAR\Analysis\toolboxes';
if ~exist('eeglab', 'file')
    addpath(fullfile(TOOLBOXPATH, 'eeglab2024.2'));
end
eeglab nogui;

%% ----------- USER SETTINGS -----------

mep_win                 = [20 40];   % ms after TMS for peak-to-peak
chanMode                = 'best';   % options: 'ch1', 'ch2', 'best'
prominence              = 1;

outlier_rejection       = false;
non_response_rejection  = true;
coil_drift_rejection    = true;

save_output             = false;
generate_plots          = false;

%% ----------- SELECT DATA INPUT -----------

% % REFTEP
reftep = true; screen = false;
DATAIN_BASE  = '';
DATAOUT_BASE = '';
DATAOUT_VIS = '';
all_subjects = load('');

% % REFTEP REPEAT
% reftep = false; screen = false; reftep_repeat = true;
% DATAIN_BASE  = '';
% DATAOUT_BASE = '';
% DATAOUT_VIS = '';
% all_subjects = load('');

%% Extract MEP peak-to-peak
% EMG channel idx REFTEP TUEBINGEN: FDI = 127; APB = 128;
% EMG channel idx REFTEP AALTO: APB = 63; FDI = 64;
% EMG channel idx SCREEN3 & SCREEN5: APB = 65; FDI = 66;

alldata_MEP_master = [];
for i = 1:length(all_subjects)
    subID = all_subjects{i};
    subids = strrep(subID, 'REFTEP_', '');
    subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
    if contains(subID, 'screen5_')
        subID = ['SCREEN5_', subids];
        subids = num2str(str2num(subids) + 500);
    end
    sub_num = horzcat('sub-', subids);

    set_file = fullfile(DATAIN_BASE, subID, horzcat(sub_num, '_EMG_processed.set'));
    if ~isfile(set_file)
        fprintf('No EMG file for %s\n', subID);
        continue;
    end

    % Load EMG dataset
    EMG = pop_loadset(set_file);

    % Downsample to 1 kHz
    % EMG = pop_resample(EMG, 1000);

    % Use Miriams function
    % parameters = MEPamp(EMG, sub_num, 1);
    % MEP_final = parameters.final_meps;

    % Extract EMG data [channels x points x trials]
    data  = EMG.data;
    times = EMG.times;  

    % Get indices for MEP window
    [~, idx_start] = min(abs(times - mep_win(1)));
    [~, idx_end]   = min(abs(times - mep_win(2)));
    mep_win_data = data(:,idx_start:idx_end,:); % chans x time x trials

    % Compute peak-to-peak per trial and channel
    nTrials = size(mep_win_data,3);
    nChans  = size(mep_win_data,1);
    ptp_all = nan(nTrials,nChans);
    loc_max = nan(nTrials,nChans);    
    loc_min = nan(nTrials,nChans);   
    for tr = 1:nTrials
        for ch = 1:nChans
            sig = squeeze(mep_win_data(ch,:,tr));
            
            % find maximum peaks
            [pks_max, temp_loc_max] = findpeaks(sig, 'MinPeakProminence', prominence);
            % find minimum peaks
            [pks_min, temp_loc_min]  = findpeaks(-sig, 'MinPeakProminence', prominence);
            
            % if at least one peak was not found, enter NaN
            if isempty(pks_min) || isempty(pks_max) % if one peak was not found
                ptp_all(tr,ch) = NaN;
                loc_max(tr,ch) = 0;
                loc_min(tr,ch) = 0;
            else
                % get peak-to-peak maximum range 
                ptp_all(tr,ch) = max(pks_max) - min(-pks_min);
                % get location of the highest peaks
                loc_max(tr,ch) = temp_loc_max(find(pks_max == max(pks_max), 1, 'first'));
                loc_min(tr,ch) = temp_loc_min(find(pks_min == max(pks_min),1,'first'));
            end
        end
    end
    
    % Decide which channel to keep
    switch chanMode
        case 'ch1'
            chanIdx = 1;
        case 'ch2'
            chanIdx = 2;
        case 'best'
            mean_ptp = mean(ptp_all,1); % avg across trials
            [~, chanIdx] = max(mean_ptp);
        otherwise
            error('Unknown chanMode: %s', chanMode);
    end
    chanUsed = EMG.chanlocs(chanIdx).labels;
    chanUsed = regexprep(chanUsed, '.*APB.*', 'APBr', 'ignorecase');
    chanUsed = regexprep(chanUsed, '.*FDI.*', 'FDIr', 'ignorecase');

    % Extract PTP from chosen channel
    MEP_uV  = ptp_all(:,chanIdx);

    % Remove outliers (IQR)
    if outlier_rejection
        Q1 = prctile(MEP_uV, 25);
        Q3 = prctile(MEP_uV, 75);
        IQR = Q3 - Q1;
        lower_bound = Q1 - 1.5 * IQR;
        upper_bound = Q3 + 1.5 * IQR;
        MEP_uV(MEP_uV < lower_bound | MEP_uV > upper_bound) = NaN;
    end

    % Remove non-responses per MEP definition (MEPs < 50 microV)
    if non_response_rejection
       MEP_uV(MEP_uV < 50) = NaN;
    end
    N_meps_rejected = sum(isnan(MEP_uV));
    idx_meps_rejected = isnan(MEP_uV); 

    % Remove coil drift rejected trials (missing marker or coil deviations)
    if coil_drift_rejection
        if strcmpi(subID, 'REFTEP_018')
            rejected_coil_drift = false(nTrials,1); % All have NaNs because of missing trigger marker files
        else
            coil_drift_data = load(fullfile(DATAIN_BASE, subID, 'stimulation_times', [subID, '_coil_drift_data.mat']));

            nTrials = coil_drift_data.n_trials_after_preprocessing;
            %kept_idx = coil_drift_data.indices.kept + 1; % Python indices start at 0
            kept_idx = coil_drift_data.indices.kept_coil_only + 1; % Python indices start at 0 #NaN trials are still included, only trials rejected on recorded coil position

            rejected_coil_drift = true(nTrials,1);
            rejected_coil_drift(kept_idx) = false;

            MEP_uV(rejected_coil_drift) = NaN;
        end
    end

    % Get indices of kept trials after removals (MEP thresholds and coil deviations)
    trials_kept = ~isnan(MEP_uV);

    % Log transform
    MEP_log = log(MEP_uV);

    % Z-score (trialwise standardization) 
    MEP_log_z = (MEP_log - nanmean(MEP_log)) ./ nanstd(MEP_log);

    % Store into struct
    alldata_MEP_master(i).subID        = subID;
    alldata_MEP_master(i).chanMode     = chanMode;
    alldata_MEP_master(i).chanUsed     = chanUsed;
    alldata_MEP_master(i).raw_meps     = squeeze(mep_win_data(chanIdx,:,:))';
    alldata_MEP_master(i).outlier_rejection = outlier_rejection;
    alldata_MEP_master(i).non_response_rejection = non_response_rejection;
    alldata_MEP_master(i).N_meps_rejected     = N_meps_rejected;
    alldata_MEP_master(i).idx_meps_rejected     = idx_meps_rejected;
    alldata_MEP_master(i).coil_drift_rejection = coil_drift_rejection;
    if coil_drift_rejection
        alldata_MEP_master(i).N_coil_drift_rejection = sum(rejected_coil_drift);
    end
    alldata_MEP_master(i).total_trials_rejected = ~trials_kept;
    alldata_MEP_master(i).final_meps = squeeze(mep_win_data(chanIdx,:, trials_kept))';
    alldata_MEP_master(i).p2p_microvolts   = MEP_uV(trials_kept);     % raw in µV
    alldata_MEP_master(i).p2p_log_microvolts = MEP_log(trials_kept);     % log-transformed
    alldata_MEP_master(i).p2p_zscored_log     = MEP_log_z(trials_kept);   % z-scored log data
    alldata_MEP_master(i).loc_max     = loc_max(trials_kept, chanIdx);   % z-scored log data
    alldata_MEP_master(i).loc_min     = loc_min(trials_kept, chanIdx);   % z-scored log data

    if generate_plots
        trials_plot = figure('units','normalized','outerposition',[0.5 0 0.5 1]);

        % plot rejection histogram
        nexttile; hold on; grid on
        h = histogram('binedges', [0:1], 'bincounts', [sum(isnan(MEP_uV))], 'FaceColor','k');
        %write numbers on top
        y = h.BinCounts;
        x = [1];
        text(0.5, y+50, string(y))
        xticks(0.5)
        xticklabels({'rejected'})
        ylim([0 800])
        ylabel('number of trials rejected')
        title(['Rejection across subject ' num2str(subID)])

        % plot EMG response in peak-to-peak window
        nexttile([2 1]); hold on; grid on;
        x = times(idx_start:idx_end);
        y = squeeze(mep_win_data(chanIdx,:,:))';
        plot(x, y(~isnan(MEP_uV),:), 'Color', [0 0 0 0.2])
        try %#ok<TRYNC>
            plot(x, y(isnan(MEP_uV),:), 'Color', [1 0 1 0.2])
        end
        ylabel('EMG response [µV]')
        xlabel('t w.r.t. TMS pulse [ms]')
        title(['EMG response for subject ' num2str(subID)])

        if save_output
            filename = fullfile(DATAOUT_VIS, [subID, '_mep_trials_plot.png']);
            exportgraphics(trials_plot, filename);
        end
    end

    fprintf('Processed %s (%d trials), channel %d (%s mode)\n', ...
        subID, nTrials, chanIdx, chanMode);
end

% Save results
if save_output
    if reftep || screen
        save(fullfile(DATAOUT_BASE, 'alldata_MEP_master.mat'), 'alldata_MEP_master');
    else
        alldata_MEP_master_repeat = alldata_MEP_master;
        save(fullfile(DATAOUT_BASE, 'alldata_MEP_master_repeat.mat'), 'alldata_MEP_master_repeat');
    end
end

return

%% Plot final MEPs

load(fullfile(DATAOUT_BASE, 'alldata_MEP_master.mat'), 'alldata_MEP_master');

subID = 'REFTEP_027'; % INPUT

subids = strrep(subID, 'REFTEP_', '');
subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
if contains(subID, 'screen5_')
    subID = ['SCREEN5_', subids];
    subids = num2str(str2num(subids) + 500);
end
sub_num = horzcat('sub-', subids);

S = alldata_MEP_master(strcmpi({alldata_MEP_master.subID}, subID));
mep_data = S.final_meps;                 % trials x time
loc_max  = S.loc_max;                    % trials x 1 (or 1 x trials)
loc_min  = S.loc_min;                    % trials x 1 (or 1 x trials)

nTrials  = size(mep_data, 1);

% x-axis for the snippet
x = linspace(20,40,101);

% make sure indices are column vectors
loc_max = loc_max(:);
loc_min = loc_min(:);

% Figure
fig = figure('Name','MEP browser');
ax1 = axes(fig); grid(ax1,'on'); hold(ax1,'on');

trial_idx = 1;

% --- initial plot
cla(ax1); hold(ax1,'on');
y = mep_data(trial_idx,:);

plot(ax1, x, y, 'k', 'LineWidth', 1.5);

% plot detected peak locations (if valid)
if ~isnan(loc_max(trial_idx)) && loc_max(trial_idx) >= 1 && loc_max(trial_idx) <= numel(x)
    plot(ax1, x(loc_max(trial_idx)), y(loc_max(trial_idx)), 'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
end
if ~isnan(loc_min(trial_idx)) && loc_min(trial_idx) >= 1 && loc_min(trial_idx) <= numel(x)
    plot(ax1, x(loc_min(trial_idx)), y(loc_min(trial_idx)), 'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
end

ylabel(ax1,'EMG response [\muV]'); xlabel(ax1,'Time [ms]');
title(ax1, sprintf('%s | %s | Trial %d / %d', sub_num, S.chanUsed, trial_idx, nTrials));
drawnow;

% --- Advance only with SPACE
while ishandle(fig) && trial_idx < nTrials
    waitforbuttonpress;
    if ~ishandle(fig); break; end

    if strcmpi(get(fig,'CurrentKey'),'space')
        trial_idx = trial_idx + 1;

        cla(ax1); hold(ax1,'on');
        y = mep_data(trial_idx,:);

        plot(ax1, x, y, 'k', 'LineWidth', 1.5);

        if ~isnan(loc_max(trial_idx)) && loc_max(trial_idx) >= 1 && loc_max(trial_idx) <= numel(x)
            plot(ax1, x(loc_max(trial_idx)), y(loc_max(trial_idx)), 'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
        end
        if ~isnan(loc_min(trial_idx)) && loc_min(trial_idx) >= 1 && loc_min(trial_idx) <= numel(x)
            plot(ax1, x(loc_min(trial_idx)), y(loc_min(trial_idx)), 'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
        end

        ylabel(ax1,'EMG response [\muV]'); xlabel(ax1,'Time [ms]');
        title(ax1, sprintf('%s | %s | Trial %d / %d', sub_num, S.chanUsed, trial_idx, nTrials));

        drawnow;
    end
end