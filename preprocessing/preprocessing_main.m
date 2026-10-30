%% T2TVAR EEG/EMG PREPROCESSING
% By Timo van Hattem (timo.van-hattem@uni-tuebingen.de)
% Brain Networks & Plasticity Lab, University of Tuebingen
%
% Description: Preprocesses TMS-EEG and EMG recordings from the
% brain-state-resolved T2TVAR study. The pipeline includes separate
% preprocessing for the pre-stimulus and post-stimulus window, including
% artifact removal, bad-channel detection, ICA-based ocular artifact
% removal, EEG cleaning, trial rejection, and EMG preprocessing.
%
% Based on Ahola et al., 2025

%% ----------- CLEAR WORKSPACE -----------

clear;          % Clear workspace variables
close all;      % Close all figure windows
clc;            % Clear command window

% Get full path of current script and move to its directory
filePath = matlab.desktop.editor.getActiveFilename;
mainpath = fileparts(filePath);   % Adjust to project root folder structure
cd(mainpath);

% Set paths to relevant toolboxes
TOOLBOXPATH = '';
if ~exist('eeglab', 'file')
    addpath(fullfile(TOOLBOXPATH, 'eeglab2024.2'));
    addpath(fullfile(TOOLBOXPATH, 'FastICA_25'));
end
eeglab nogui

%% ----------- USER SETTINGS -----------

% EEG temporal windows
% High-pass filter applied before ICA and pre-stimulus bad-channel detection.
pre_highpass = 2; % Hz
% Baseline interval for post-stimulus EEG correction.
EEG_post_baseline_range = [-20 -10]; % ms
% Initial post-stimulus extraction window.
post_range = [-0.02 0.3]; % s
% Post-stimulus window used for bad-channel and trial detection.
post_range_true = [0.01 0.15]; % s
% Final post-stimulus EEG epoch.
post_range_final = [-0.02 0.15]; % s
% Initial and final pre-stimulus EEG epochs.
pre_range_init = [-1.25 -0.025]; % s
pre_range_final = [-1.045 -0.045]; % s

% Pulse-artifact removal and interpolation
% First pulse-artifact removal and interpolation.
pulse_artifact_window1 = [-5 6]; % ms
interpolation_window1 = [1 1]; % ms
% Second pulse-artifact removal and interpolation (wider window).
pulse_artifact_window2 = [-5 8]; % ms
interpolation_window2 = [1 1]; % ms

% Bad-channel detection
% Adjustment factors for rejection thresholds applied to frontal channels.
frontal_threshold_shift = [1.5 0.75];
% Peak-based rejection thresholds for pre-stimulation EEG.
pre_peak_threshes_large = [70 1 3];
pre_peak_threshes_small = [2 4 0];
% Z-score thresholds for bad-channel detection in pre- and post-stimulation EEG.
rejchan_threshold_pre = 2.5;
rejchan_threshold_post = 3;

% ICA-based ocular artifact removal
% Number of ICA components estimated after PCA compression.
n_ica_components = 35;
% ICLabel probability threshold for rejecting eye-related components.
ica_prob_threshold = 0.75;
% Figure dimensions for ICA component topographies.
ica_plot_window = [5 9];

% Pre-stimulation EEG filtering and resampling
% Butterworth filter order.
filter_order = 4;
% Band-pass and power-line noise filter ranges.
filter_range = [2 90]; % Hz
powerline_range = [48 52]; % Hz (band-stop)
% Target sampling frequency after resampling.
resample_to = 1000; % Hz

% EEG artifact suppression
% SOUND parameters: regularization parameter and number of iterations.
sound_lambda = 0.05;
sound_n_iters = 5;
% SSP-SIR percentile used for principal-component selection.
sspsir_percentile = 90;

% Trial rejection
% Global and local rejection thresholds for pre-stimulation EEG.
rejtrial_threshold_global_pre = 1.5;
rejtrial_thresholds_local_pre = [6 5];
nums_over_rejtrial_thresholds_local_pre = [0 1];
% Global and local rejection thresholds for post-stimulation EEG.
rejtrial_threshold_global_post = 2.5;
rejtrial_thresholds_local_post = [6 5];
nums_over_rejtrial_thresholds_local_post = [0 1];

% EMG preprocessing and trial rejection
% EMG epoch and pre-stimulation muscle activity assessment window.
time_window_emg = [-0.1 0.1]; % s
emg_preact_window = [-0.1 -0.01]; % s
% Primary threshold for rejecting trials with excessive pre-stimulation
% EMG activity, with additional thresholds for sensitivity assessment.
emg_preact_threshold_muv = 50; % µV
threshes_extra_emg = emg_preact_threshold_muv + [10 20 30 40 50];

set(0, 'DefaultFigureVisible', 'off');

%% ----------- SELECT DATA INPUT -----------

datapath_base = '';
datapath_out = '';
subject_names = load('');

save_output = false;

%% ----------- PREPROCESSING EEG/EMG -----------

for index = 1:length(subject_names)

    %% Get subject ID

    reftep_subject = char(subject_names(index))
    subids = strrep(reftep_subject, 'REFTEP_', '');
    sub_num = horzcat('sub-', subids);

    %% Define file paths and directories

    directory_path_in = char(fullfile(datapath_base,reftep_subject,"\")); %the directory of the EEG data to be loaded
    directory_path_out = char(fullfile(datapath_out,reftep_subject,"\")); %the directory of the EEG data to be loaded
    if ~exist(directory_path_out, 'dir')
        mkdir(directory_path_out);
    end
    filename_tmseeg = char(strcat(sub_num,'_task-tep_all_eeg.set')); %.set file name of the eeg data

    %% Load data, and divide EEG and EMG:

    % Divide the data into epochs
    EEG_and_EMG_epoched = pop_loadset(filename_tmseeg,directory_path_in); %load the data
    channel_names_original = {EEG_and_EMG_epoched.chanlocs.labels}; %all channel names
    channel_names_original_lower = lower(channel_names_original); %all channel names lowered

    block_identifiers_trials = load(fullfile(directory_path_in, horzcat(sub_num, '_block_identifiers.mat'))).block_identifiers_trials;

    emg_channels = channel_names_original(contains(channel_names_original_lower,'emg') | ...
        contains(channel_names_original_lower,'fdi') | contains(channel_names_original_lower,'apb')); %select emg channel names
    if length(emg_channels) ~=2
        disp("more or less than 2 emg channels detected")
    end
    EMG = pop_select(EEG_and_EMG_epoched, 'channel', emg_channels); %select only EMG channels
    for i=1:EMG.nbchan
        EMG.chanlocs(i).type = 'EMG'; %set EMG channel type
    end

    eeg_channels = channel_names_original(~contains(channel_names_original_lower,'emg') & ...
        ~contains(channel_names_original_lower,'fdi') & ~contains(channel_names_original_lower,'apb')); %select eeg channel names
    EEG = pop_select(EEG_and_EMG_epoched, 'channel', eeg_channels); %select eeg channels
    for i=1:EEG.nbchan
        EEG.chanlocs(i).type = 'EEG'; %set EEG channel type
    end
    clear EEG_and_EMG_epoched %not needed anymore, save data

    %% Split pre-stimulus and post-stimulus window

    EEG_pre = pop_select(EEG, 'time', pre_range_init);
    EEG_post = pop_select(EEG, 'time', post_range);
    post_raw_name = [sub_num,'_EEG_post_with_all_chans_epoched.set'];
    pop_saveset(EEG_post, 'filename',post_raw_name,'filepath',directory_path_out);
    pre_raw_name = [sub_num,'_EEG_pre_with_all_chans_epoched.set'];
    pop_saveset(EEG_pre, 'filename',pre_raw_name,'filepath',directory_path_out);
    EEG_pre_with_all_channels = EEG_pre;
    clear EEG

    %% Remove the TMS pulse artifact from post-stimulus window and interpolate missing data

    % Remove Pulse artefact
    EEG_post = pop_tesa_removedata(EEG_post, pulse_artifact_window1);
    % Interpolate the missing pulse data
    EEG_post = pop_tesa_interpdata(EEG_post, 'cubic', interpolation_window1);

    %% Detect bad channels from the pre-stimulus window and post-stimulus window

    EEG_post_true = pop_rmbase(EEG_post, EEG_post_baseline_range); %remove baseline and then select time
    EEG_post_true = pop_select(EEG_post_true, 'time', post_range_true); %post-stim eeg
    EEG_pre_filtered = pop_tesa_filtbutter(EEG_pre, pre_highpass, [], filter_order, 'highpass');
    EEG_pre_filtered = pop_select(EEG_pre_filtered, 'time', pre_range_final);
    bad_channels_pre = find_bad_channels(EEG_pre_filtered, rejchan_threshold_pre, frontal_threshold_shift, pre_peak_threshes_large, pre_peak_threshes_small);
    bad_channels_post = find_bad_channels(EEG_post_true, rejchan_threshold_post, [0 0], [], []);
    clear EEG_post_true EEG_pre_filtered %not used anymore, save memory

    % Combine bad channels detected from pre-stimulus window and post-stimulus window
    union_bad_elecs = union(bad_channels_pre, bad_channels_post); %combine bad channels detected from pre and post
    chan_names = {EEG_pre.chanlocs.labels}; %existing EEG channel names
    bad_chans = chan_names(union_bad_elecs) %bad channel names

    % Remove bad channels from both pre and post data and set datas to average reference
    EEG_pre = pop_select(EEG_pre, 'nochannel', bad_chans); %remove bad channels from the data
    EEG_post = pop_select(EEG_post, 'nochannel', bad_chans); %remove bad channels from the data
    EEG_pre = pop_reref(EEG_pre, []); %average reference the EEG data before ICA
    EEG_post = pop_reref(EEG_post, []); %average reference the EEG data before ICA

    %% Fit Independent Component Analysis (ICA) on pre-stimulus window

    % Highpass filter the data before ICA
    EEG_pre_filtered = pop_tesa_filtbutter(EEG_pre, pre_highpass, [], filter_order, 'highpass'); %highpass filter the data before ICA
    EEG_pre_filtered = pop_select(EEG_pre_filtered, 'time', pre_range_final);
    EEG_pre_filtered = pop_tesa_pcacompress(EEG_pre_filtered, 'compVal', n_ica_components, 'plot','off' );
    EEG_pre_filtered = pop_reref(EEG_pre_filtered, []); %average reference the EEG data before ICA
    EEG_pre_filtered = pop_tesa_fastica(EEG_pre_filtered, 'approach','symm','g','tanh','stabilization','off'); %run ICA
    % Run ICLabel
    EEG_pre_filtered = iclabel(EEG_pre_filtered); %run iclabel
    ic_probabilities = EEG_pre_filtered.etc.ic_classification.ICLabel.classifications; %n_channels x classes (or labels)
    possible_labels = EEG_pre_filtered.etc.ic_classification.ICLabel.classes; %possible IC labels

    % Find the IC probabilities of being an 'Eye' component
    index_of_label = find(strcmp(possible_labels, 'Eye')); %index of the 'Eye' label

    % Get the component indices that exceed the threshold
    component_indices_over_thresh = find(ic_probabilities(:,index_of_label) > ica_prob_threshold)
    ic_probabilities(:,index_of_label)

    %% Set the ICA properties from pre-stimulus window to the unfiltered pre- and post-stimulus data and remove the components

    % Set weights to EEG_pre
    EEG_pre.icaweights = EEG_pre_filtered.icaweights;
    EEG_pre.icasphere = EEG_pre_filtered.icasphere;
    EEG_pre.icawinv = EEG_pre_filtered.icawinv;
    EEG_pre.icachansind = EEG_pre_filtered.icachansind;

    % Set weights to EEG_post
    EEG_post.icaweights = EEG_pre_filtered.icaweights;
    EEG_post.icasphere = EEG_pre_filtered.icasphere;
    EEG_post.icawinv = EEG_pre_filtered.icawinv;
    EEG_post.icachansind = EEG_pre_filtered.icachansind;

    % Remove the selected components from EEG_pre and EEG_post
    EEG_pre = pop_subcomp(EEG_pre, component_indices_over_thresh, 0);
    EEG_post = pop_subcomp(EEG_post, component_indices_over_thresh, 0);

    %% Reconstruct bad channels in the pre-stimuls window with spherical interpolation

    EEG_pre = pop_interp(EEG_pre, EEG_pre_with_all_channels.chanlocs, 'spherical'); %interpolate bad channels in pre-stim

    %% Filter and downsample the pre-stimulus EEG

    EEG_pre = pop_tesa_filtbutter(EEG_pre, filter_range(1), filter_range(2), filter_order, 'bandpass');
    EEG_pre = pop_tesa_filtbutter(EEG_pre, powerline_range(1), powerline_range(2), filter_order, 'bandstop');
    EEG_pre = pop_resample(EEG_pre, resample_to);

    % Crop pre-stim window to its final time range and rereference
    EEG_pre = pop_select(EEG_pre, 'time', pre_range_final); %crop pre-stim to its final time range
    EEG_pre = pop_reref(EEG_pre,[]);

    %% Remove bad trials detected from pre-stimulus window from both pre-stimulus and post-stimulus EEG

    % Find bad trials
    bad_trials_pre = find_bad_trials(EEG_pre, rejtrial_threshold_global_pre, rejtrial_thresholds_local_pre, nums_over_rejtrial_thresholds_local_pre);
    % Remove bad trials from both windows
    if ~isempty(bad_trials_pre)
        EEG_pre = pop_select(EEG_pre,'notrial',bad_trials_pre);
        EEG_post = pop_select(EEG_post,'notrial',bad_trials_pre);
        EMG = pop_select(EMG,'notrial',bad_trials_pre);
        block_identifiers_trials(bad_trials_pre) = [];
    end

    %% SOUND on post-stimulus window: suppress noise and reconstruct bad channels

    EEG_post = pop_rmbase(EEG_post, EEG_post_baseline_range);
    % Run SOUND
    EEG_post = pop_tesa_sound(EEG_post, 'lambdaValue', sound_lambda, 'iter', sound_n_iters, ...
        'leadfieldInFile' ,[], 'leadfieldChansFile', [], 'replaceChans', char(fullfile(directory_path_out,post_raw_name)), 'multipleConds', []);

    %% SSP-SIR on post-stimulus window: suppress muscle artifacts

    % Run SSP-SIR
    EEG_post = pop_tesa_sspsir(EEG_post, 'artScale','automatic','PC', {'data',[sspsir_percentile]});

    %% Crop final post-stimulus window
    %take care of pulse artifact now on a widened window (it is highly
    %likely that the previous interpolation of the window interpolated
    %it with artifactual data)
    EEG_post = pop_tesa_removedata(EEG_post, pulse_artifact_window2);
    EEG_post = pop_tesa_interpdata(EEG_post, 'cubic', interpolation_window2);
    EEG_post = pop_rmbase(EEG_post, EEG_post_baseline_range);
    EEG_post = pop_reref(EEG_post,[]);


    %% Remove bad trials detected from post-stimulus window from both pre-stimulus and post-stimulus EEG

    EEG_post_in_post_range = pop_select(EEG_post, 'time', post_range_true);
    bad_trials_post = find_bad_trials(EEG_post_in_post_range, rejtrial_threshold_global_post, rejtrial_thresholds_local_post, nums_over_rejtrial_thresholds_local_post);

    if ~isempty(bad_trials_post)
        EEG_pre = pop_select(EEG_pre,'notrial',bad_trials_post);
        EEG_post = pop_select(EEG_post,'notrial',bad_trials_post);
        EMG = pop_select(EMG,'notrial',bad_trials_post);
        block_identifiers_trials(bad_trials_post) = [];
    end

    %% Crop final post-stimulus window
    EEG_post = pop_select(EEG_post, 'time', post_range_final);

    %% Process EMG data
    % Crop EMG to the analysis window, correct the baseline, and remove linear
    % trends from the pre-stimulation muscle activity interval.
    EMG = pop_select(EMG, 'time', time_window_emg);
    EMG = pop_rmbase(EMG, emg_preact_window*1000, []);
    EMG = pop_tesa_detrend(EMG, 'linear', emg_preact_window*1000);

    % Extract pre-stimulation EMG activity for trial rejection.
    % Add one sampling interval to include the final time point.
    epsilon = 1/EMG.srate;
    EMG_pre = pop_select(EMG, 'time', emg_preact_window + [0 epsilon]);

    % Identify trials with excessive pre-stimulation muscle activity.
    bad_emg_trials = find_bad_emg_trials(EMG_pre, emg_preact_threshold_muv, threshes_extra_emg);

    % Remove rejected trials from all datasets to maintain trial-wise
    % correspondence between EEG, EMG, and block identifiers.
    if ~isempty(bad_emg_trials)
        EEG_pre = pop_select(EEG_pre, 'notrial', bad_emg_trials);
        EEG_post = pop_select(EEG_post, 'notrial', bad_emg_trials);
        EMG = pop_select(EMG, 'notrial', bad_emg_trials);
        block_identifiers_trials(bad_emg_trials) = [];
    end

    %% Sanity checks

    EMG = eeg_checkset(EMG);
    EEG_pre = eeg_checkset(EEG_pre);
    EEG_post = eeg_checkset(EEG_post);

    fields_to_remove = {'event','eventlatency','eventtype','eventurevent'};
    trimmed_epochinfo_EEG_post = rmfield(EEG_post.epoch,fields_to_remove);
    trimmed_epochinfo_EEG_pre = rmfield(EEG_pre.epoch,fields_to_remove);

    if ~isequal(EMG.epoch, EEG_post.epoch) || ~isequal(trimmed_epochinfo_EEG_post, trimmed_epochinfo_EEG_pre)
        disp("structures have differing epoch information")
        return
    end

    if length(block_identifiers_trials) ~= EEG_pre.trials || ...
            length(block_identifiers_trials) ~= EEG_post.trials || ...
            length(block_identifiers_trials) ~= EMG.trials
        error('Mismatch: block_identifiers_trials length (%d) does not match number of trials (EEG_pre=%d, EEG_post=%d, EMG=%d)', ...
            length(block_identifiers_trials), EEG_pre.trials, EEG_post.trials, EMG.trials);
    end

    %% Save data

    if save_output
        pop_saveset(EMG, 'filename',[sub_num,'_EMG_processed.set'],'filepath',directory_path_out);
        pop_saveset(EEG_pre, 'filename',[sub_num,'_EEG_pre_processed.set'],'filepath',directory_path_out);
        pop_saveset(EEG_post, 'filename',[sub_num,'_EEG_post_processed.set'],'filepath',directory_path_out);
        save(fullfile(directory_path_out, [sub_num, '_block_identifiers_trials.mat']), 'block_identifiers_trials');
        urevent_identifiers_trials = [EEG_post.epoch.eventurevent];
        save(fullfile(directory_path_out, [sub_num, '_urevent_identifiers_trials.mat']), 'urevent_identifiers_trials');
    end
end