clear;          % Clear workspace variables
close all;      % Close all figure windows
clc;            % Clear command window

% Get full path of current script and move to its directory
filePath = matlab.desktop.editor.getActiveFilename;
mainpath = fileparts(filePath);   % Adjust to project root folder structure
cd(mainpath);

% Set paths to relevant toolboxes
TOOLBOXPATH = '';
addpath(fullfile(TOOLBOXPATH, 'eeglab2024.2'));
addpath(fullfile(TOOLBOXPATH, 'FastICA_25'));
eeglab nogui;

datapath_base = ''; %where to store the data
datapath_out = '';

%% Calculate PSD
subject_names = load('');

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
    filename_tmseeg = char(strcat(sub_num, '_task-rest_eeg.set')); %.set file name of the eeg data

    EEG = pop_loadset(filename_tmseeg,directory_path_in); %load the data
    channel_names_original = {EEG.chanlocs.labels}; %all channel names
    channel_names_original_lower = lower(channel_names_original); %all channel names lowered
    channels = channel_names_original(~contains(channel_names_original_lower,'input')); %select proper channels, this one is bad in tuebingen
   
    EMG_ch.FDI = 127;
    EMG_ch.APB = 128;

    %% Separate EMG and EEG channels

    EEG = pop_select(EEG, 'channel', channels);
    EEG = pop_select(EEG, 'rmchannel', [EMG_ch.FDI EMG_ch.APB] );

    %% Downsample data from 5000 to 1000 Hz
    EEG = pop_resample(EEG, 1000);

    %% C3-Hjorth Filter
    hjorth.channel = 'C3'
    hjorth.electrodes = {'FCC5h', 'FCC3h', 'CCP5h', 'CCP3h'};
    % hjorth.electrodes = {'FC5', 'FC1', 'CP5', 'CP1'}
    electrode_labs = struct2cell(EEG.chanlocs);

    hjorth.channel_idx = find(strcmpi(hjorth.channel, electrode_labs(1,:)));
    for i = 1:length(hjorth.electrodes)
        hjorth.electrode_idx(i) = find(strcmpi(hjorth.electrodes(i), electrode_labs(1,:)));
    end
    EEG_lap = EEG;

    lap_filter = zeros(size(EEG.data,1),1);
    lap_filter(hjorth.channel_idx) = 1;
    lap_filter(hjorth.electrode_idx) = -1/4;

    temp = reshape(EEG_lap.data, 126, []);       % 126 for Tue
    temp = temp.*lap_filter;
    temp = sum(temp, 1);
    temp = reshape(temp, length(EEG_lap.times), []);
    
    EEG_lap.data = temp;
    EEG_lap.chanlocs = 'C3-centered laplacian';
    EEG_lap.nbchan = 1;

    %% Calculate PSD
    Fs = EEG_lap.srate;
    f_range = [4 50];
    win_len  = round(4 * Fs);        % 500 ms
    noverlap = round(0.50 * win_len);
    window   = hanning(win_len);

    nFFT = win_len; 
    x = EEG_lap.data;   % channel 1, continuous

    % --- Welch PSD ---
    [Pxx, f] = pwelch(x, window, noverlap, nFFT, Fs);
    idx   = f >= f_range(1) & f <= f_range(2);
    psd   = Pxx(idx);
    freqs = f(idx);
    
    save(fullfile(directory_path_out, [sub_num, '_rsEEG_PSD.mat']), 'psd', 'freqs');
end







