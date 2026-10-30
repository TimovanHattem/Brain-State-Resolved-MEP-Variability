clear;          % Clear workspace variables
close all;      % Close all figure windows
clc;            % Clear command window

% Get full path of current script and move to its directory
filePath = matlab.desktop.editor.getActiveFilename;
mainpath = fileparts(filePath);   % Adjust to project root folder structure
cd(mainpath);

% Set paths to relevant toolboxes
TOOLBOXPATH = 'D:\2025-08 T2TVAR\Analysis\toolboxes\';
addpath(fullfile(TOOLBOXPATH, 'eeglab2024.2'));
addpath(fullfile(TOOLBOXPATH, 'FastICA_25'));
eeglab nogui;
addpath(genpath(['D:\2025-08 T2TVAR\Analysis\toolboxes\phastimate']))

datapath_base = 'Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\converted'; %where to store the data
datapath_out = 'Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\preprocessed_alternative_phase';

save_output = true;

%% Estimate phase
subject_names = load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\reftep_active_subject_list.mat').reftep_active_subject_list;

for index = 1:29 % TUBINGEN ONLY
    reftep_subject = char(subject_names(index));
    subids = strrep(reftep_subject, 'REFTEP_', '');
    sub_num = horzcat('sub-', subids);

    %% Load data
    directory_path_in = char(fullfile(datapath_base,reftep_subject,"\")); %the directory of the EEG data to be loaded
    directory_path_out = char(fullfile(datapath_out,reftep_subject,"\")); %the directory of the EEG data to be loaded
    if ~exist(directory_path_out, 'dir')
        mkdir(directory_path_out);
    end
    filename_tmseeg = char(strcat(sub_num, '_task-tep_all_eeg.set')); %.set file name of the eeg data

    EEG = pop_loadset(filename_tmseeg, directory_path_in); %load the data
    
    channel_names_original = {EEG.chanlocs.labels}; %all channel names
    channel_names_original_lower = lower(channel_names_original); %all channel names lowered
    channels = channel_names_original(~contains(channel_names_original_lower,'input')); %select proper channels, this one is bad in tuebingen

    EMG_ch.FDI = 127;
    EMG_ch.APB = 128;

    % Load individual mu peak frequencies
    load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\individual_mu_peak_frequencies.mat');
    mu_cf_primary = mu_cf_primary(:);
    subjects = subjects(:);  
    mu_peak_ALL = table(subjects, mu_cf_primary, 'VariableNames', {'id', 'mu_peak'});

    %% Separate EMG and EEG channels
    EEG = pop_select(EEG, 'channel', channels);
    EMG = pop_select(EEG, 'channel', [EMG_ch.FDI EMG_ch.APB]);
    EEG = pop_select(EEG, 'rmchannel', [EMG_ch.FDI EMG_ch.APB] );
    
    %% Baseline correct EEG data
    EEG_epoched = pop_select(EEG, 'time', [-1.005, -0.005]);
    data = EEG_epoched.data;
    baselines = mean(data,2);

    data_rm = data - baselines;
    EEG_epoched.data = data_rm;

    clear data data_rm baselines

    %% Detrend
    t.detrend = [min(EEG_epoched.times), min(-5, max(EEG_epoched.times))];
    EEG_epoched = tesa_detrend(EEG_epoched, 'linear', t.detrend);
    EEG = EEG_epoched;

    clear EEG_epoched
    
    %% Find bad channels and exclude trials where the threshold was crossed for any channel in the laplacian
    hjorth.channel = 'C3'
    hjorth.electrodes = {'FCC5h', 'FCC3h', 'CCP5h', 'CCP3h'};
    % hjorth.electrodes = {'FC5', 'FC1', 'CP5', 'CP1'}
    electrode_labs = struct2cell(EEG.chanlocs);

    hjorth.channel_idx = find(strcmpi(hjorth.channel, electrode_labs(1,:)));
    for i = 1:length(hjorth.electrodes)
        hjorth.electrode_idx(i) = find(strcmpi(hjorth.electrodes(i), electrode_labs(1,:)));
    end

    % define time relevant for exclusion
    rejection.threshold_level = 150;

    % find relevant channels per trial
    threshold_marker = range(EEG.data(:,:,:),2)>rejection.threshold_level;
    % find trials that cross the threshold in any relevant channel for the
    % montage
    temp = squeeze(threshold_marker([hjorth.electrode_idx hjorth.channel_idx],:,:));

    % note down trials where the montage is not reliable
    rejection.threshold = sum(temp) >= 1;
    excluded_trials = find(rejection.threshold == 1);
    kept_trials = find(rejection.threshold == 0);
    fprintf('%.0f trials rejected by threshold. \n', sum(rejection.threshold))

    % reject selected trials
    EEG = pop_select(EEG, 'trial', find(~rejection.threshold));
    EMG = pop_select(EMG, 'trial', find(~rejection.threshold));

    %% Reject preinnervation EMG
    addpath('D:\2025-08 T2TVAR\Analysis\preprocessing_REFTEP_M1')
    EMG_temp = pop_select(EMG, 'time', [-0.1 0.1]);
    EMG_temp = pop_rmbase(EMG_temp, [-0.1 -0.01]*1000 ,[]);
    EMG_temp = pop_tesa_detrend(EMG_temp, 'linear',[-0.1 -0.01]*1000); %detrend the preact window
    epsilon = 1/EMG_temp.srate;
    EMG_pre = pop_select(EMG_temp, 'time', [-0.1 -0.01] + [0 epsilon]);
    
    bad_emg_trials = find_bad_emg_trials(EMG_pre, 50, 50 + [10 20 30 40 50]) %bad emg trials
    if ~isempty(bad_emg_trials)
        EEG = pop_select(EEG, 'notrial', bad_emg_trials);
        EMG = pop_select(EMG,'notrial', bad_emg_trials);
    end

    %% Preprocess EMG data
    EMG_epoched = pop_resample(EMG, 1000);

    mep_win                 = [20 40];   % ms after TMS for peak-to-peak
    chanMode                = 'best';   % options: 'ch1', 'ch2', 'best'
    prominence              = 1;

    outlier_rejection       = false;
    non_response_rejection  = true;

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

    % Get indices of kept trials after removals (MEP thresholds and coil deviations)
    trials_kept = ~isnan(MEP_uV);

    MEP_uV = MEP_uV(trials_kept);

    % Log transform
    MEP_log = log(MEP_uV);

    % Z-score (trialwise standardization) 
    MEP_log_z = (MEP_log - nanmean(MEP_log)) ./ nanstd(MEP_log);
    
    % reject selected trials
    EEG = pop_select(EEG, 'trial', find(trials_kept));
    EMG = pop_select(EMG, 'trial', find(trials_kept));

    %% Downsample data from 5000 to 1000 Hz
    EEG = pop_resample(EEG, 1000);

    %% Bandpass filter data 
    bp.frequencies = [9 13];
    EEG = pop_eegfiltnew(EEG, 'locutoff', bp.frequencies(1),'hicutoff', bp.frequencies(2));

    %% Create Laplacian filter
    EEG_lap = EEG;

    lap_filter = zeros(size(EEG.data,1),1);
    lap_filter(hjorth.channel_idx) = 1;
    lap_filter(hjorth.electrode_idx) = -1/4;

    temp = reshape(EEG_lap.data, 126, []);
    temp = temp.*lap_filter;
    temp = sum(temp, 1);
    temp = reshape(temp, length(EEG_lap.times), []);

    EEG_lap.data = temp;
    EEG_lap.chanlocs = 'C3-centered laplacian';
    EEG_lap.nbchan = 1;
    
    %% Phastimate 
    mu_row = find(strcmpi(mu_peak_ALL.id, reftep_subject));
    foi = table2array(mu_peak_ALL(mu_row, 2));
    if isnan(foi)
        foi = 10;
    end
    ds_fs = 1000;   
    D = designfilt('bandpassfir', 'FilterOrder', 80, 'StopbandFrequency1', foi-5 , 'PassbandFrequency1', foi-1, ...
        'PassbandFrequency2', foi+1, 'StopbandFrequency2', foi+5, 'SampleRate', ds_fs, 'DesignMethod', 'ls');

    phast.edge = 65;
    phast.ord = 30;
    phast.hilbertwindow = 128;
    phast.offset_correction = 4;

    phase = phastimate_adapted(EEG_lap.data, phast.edge, phast.ord, phast.hilbertwindow, phast.offset_correction);

    %% Generate table
    nPhaseBins = 8;
    phase_values = phase;
    NaN_idx = find(isnan(phase_values));
    if sum(isnan(phase_values))
        if sum(isnan(phase_values)) > 10
            return
        end
    end
    for i = 1:length(NaN_idx)
        random_phase_value = -pi + (2*pi)*rand;
        phase_values(NaN_idx(i)) = random_phase_value;
    end

    % Create centered bins
    binWidth = 2*pi / nPhaseBins;
    centers = -pi : binWidth : pi - binWidth;
    edges = centers - binWidth/2;
    edges(end+1) = edges(end) + binWidth;
    edges(1) = [];
    % edges = linspace(-pi, pi, nBins+1); % 9 edges = 8 bins NOT CENTERED

    % Assign phases to bins
    [phase_counts, phase_edges, phase_bin] = histcounts(phase_values, edges);
    phase_bin = phase_bin + 1; % WATCH OUT NaNs phase_values NOW GET TURNED INTO TROUGH

    % Get phase labels
    if nPhaseBins == 4
        labels = {"Trough","Rising","Peak","Falling"};
    elseif nPhaseBins == 8
        labels = {"Trough","Early Rising","Rising","Late Rising", ...
            "Peak","Early Falling","Falling","Late Falling"};
    end
    phase_labels = labels(phase_bin)';

    if length(phase) ~= length(MEP_uV)
        error('mismatch')
    else
       assign_tbl = table(repmat(1:length(phase),1)', phase_values', phase_bin',  phase_labels, MEP_uV, MEP_log, MEP_log_z, ...
        'VariableNames', {'Trial_index','Phase_C3H','Phase_Bin', 'Phase_Label', 'MEP_uV', 'MEP_log', 'MEP_log_z'});
    end

    if save_output
        output_filename = [sub_num, '_sensor_Phase_alternative_preprocessing_1_8phasebins_df_splits.mat'];
        output_path = fullfile(directory_path_out, output_filename);
        save(output_path, 'assign_tbl');
    end
end






