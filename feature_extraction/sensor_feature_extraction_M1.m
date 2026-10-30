%% T2TVAR SPLIT TRIALS BASED ON PRE-STIMULUS FEATURES
% By Timo van Hattem (timo.van-hattem@uni-tuebingen.de)
% Brain Networks & Plasticity Lab, University of Tuebingen
%
% Description: Extract sensor-level EEG features from pre-stimulation
% recordings, including spectral power, connectivity, oscillatory phase,
% and phase-amplitude coupling (PAC), using conventional and Hjorth
% montages.

%% ----------- CLEAR WORKSPACE -----------

clear
clc
close all
tic

% Set paths to relevant toolboxes
TOOLBOXPATH = '';
addpath(fullfile(TOOLBOXPATH, 'eeglab2024.2'));
addpath(fullfile(TOOLBOXPATH, 'fieldtrip-20240603'));
addpath(genpath(''));
eeglab nogui
% Initialize fieldtrip
ft_defaults

% Define frequency bands
% freq_bands = {
%     'delta', 0.5, 4.0;
%     'theta', 4.0, 8.0;
%     'alpha', 8.0, 12.0;
%     'low_beta', 13.0, 18.0;
%     'high_beta', 19.0, 35.0;
%     'low_gamma', 36.0, 58.0;
%     'high_gamma', 58.0, 100.0
%     };
freq_bands = {
    'theta', 4.0, 8.0;
    'alpha', 8.0, 13.0;
    'beta', 13.0, 30.0;
    'gamma', 30.0, 80.0
    };

% Define the combinations for PAC calculation
PAC_combinations = {
    'theta', 'low_gamma';
    'theta', 'high_gamma';
    'alpha', 'low_gamma';
    'alpha', 'high_gamma';
    'low_beta', 'low_gamma';
    'low_beta', 'high_gamma';
    'high_beta', 'low_gamma';
    'high_beta', 'high_gamma'
    };

% Define phase of interest
phase_of_interest = 'alpha'; 

% --- Define channel groups (Hjorth montages)
channels_C3    = {'C3', 'FC1', 'FC5', 'CP1', 'CP5'};
channels_C4    = {'C4', 'FC2', 'FC6', 'CP2', 'CP6'};
channels_FCC3h = {'FCC3h', 'FFC1h', 'FFC5h', 'CCP1h', 'CCP5h'};
channels_FCC4h = {'FCC4h', 'FFC2h', 'FFC6h', 'CCP2h', 'CCP6h'};
% Define Hjorth weights
weights = [1, -1/4, -1/4, -1/4, -1/4];

%% Define input data

% % REFTEP
DATAIN = '';     
DATAOUT = '';
all_subjects = load('');
load('');
   
% % REFTEP REP
% DATAIN = '';     
% DATAOUT = '';
% all_subjects = load('');

%% Select feature extraction

psd = false;
connectivity = false;
hjorthpower = false;
hjorthphase = true;
hjorthconnectivity = false;
hjorthpac = false;

%% SENSOR FEATURE EXTRACTION

for i = 1:length(all_subjects)
    subID = all_subjects{i};
    subids = strrep(subID, 'REFTEP_', '');
    subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
    if contains(subID, 'screen5_')
        subID = ['SCREEN5_', subids];
        subids = num2str(str2num(subids) + 500);
    end
    sub_num = horzcat('sub-', subids);

    fprintf('\nLoading data: %s\n', subID)

    EEG = pop_loadset('filename', [sub_num,'_EEG_pre_processed.set'], 'filepath', fullfile(DATAIN, subID));

    % EEG.data: [channels x timepoints x trials]
    num_channels = size(EEG.data, 1);
    num_samples = size(EEG.data, 2);
    num_trials = size(EEG.data, 3);

    cfg = [];
    cfg.dataset = 'mydata.set';
    ft_data = [];
    ft_data.label = {EEG.chanlocs.labels};
    ft_data.fsample = EEG.srate;
    ft_data.trial = cell(1, num_trials);
    ft_data.time  = cell(1, num_trials);

    for t = 1:num_trials
        ft_data.trial{t} = EEG.data(:,:,t);          % [channel x timepoints]
        ft_data.time{t}  = EEG.times/1000;           % [1 x timepoints] in seconds
    end
    ft_data.trialinfo = (1:num_trials)';

    % Add electrode info (elec) if needed for plotting
    ft_data.elec = [];
    ft_data.elec.label = ft_data.label;
    ft_data.elec.chanpos = [[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]']; % 3D positions
    ft_data.elec.elecpos = ft_data.elec.chanpos;
    ft_data.elec.unit = 'mm';
    ft_data.elec.type = 'EEG';
    ft_data.elec.chantype = repmat({'EEG'}, length(ft_data.label), 1);
    ft_data.cfg = [];

    epochs = ft_data;
    fsample = epochs.fsample; % sampling frequency in Hz
    n_trials = length(ft_data.trial);
    n_samples = fsample + 1;

    data_base = ft_data;
    data_base.sampleinfo = [(1:n_trials)' * n_samples - n_samples + 1, (1:n_trials)' * n_samples];

    % Start a parallel pool (if not already started)
    if isempty(gcp('nocreate'))
        parpool; % Start with default settings
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PSD
    if psd
        fprintf('\nPower extraction: %s\n', subID)

        data = data_base;

        % Categorize channels based on Brain regions and Hemisphere
        brain_region_map = categorize_channels(data.label);
        brain_side_map = categorize_channels_rl(data.label);

        % Initialize cell arrays for brain regions and sides
        brain_regions = cell(1, length(data.label));
        brain_sides = cell(1, length(data.label));

        % Populate brain_regions and brain_sides based on the maps
        for i = 1:length(data.label)
            ch = data.label{i};
            brain_regions{i} = brain_region_map(ch);
            brain_sides{i} = brain_side_map(ch);
        end

        % Initialize results structures
        results = struct();

        % Loop through each frequency band
        for i = 1:size(freq_bands, 1)
            band = freq_bands{i, 1};
            low_freq = freq_bands{i, 2};
            high_freq = freq_bands{i, 3};

            % Set taper based on frequency band
            if strcmp(band, 'low_gamma') || strcmp(band, 'high_gamma')
                taper = 'dpss';
                tapsmofrq = 2; %for gamma
                step_size = 1; % Larger step for higher frequency bands
            else
                taper = 'hanning';
                tapsmofrq = []; % Not needed for hanning taper
                step_size = 1; % Smaller step for lower frequency bands
            end

            % Configuration for ft_freqanalysis
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.output = 'pow';
            cfg.keeptrials = 'yes';
            cfg.channel = data.label;
            cfg.pad = 'nextpow2';
            % Determine frequencies of interest for this band
            cfg.foi = low_freq:step_size:high_freq;
            cfg.taper = taper;
            if ~isempty(tapsmofrq)
                cfg.tapsmofrq = tapsmofrq;
            end

            % Perform frequency analysis
            temp_freq = ft_freqanalysis(cfg, data);

            % Store results
            results.(band) = temp_freq;
        end

        % Initialize results table
        n_trials = length(data.trial);
        n_channels = length(data.label);
        n_freq_bands = size(freq_bands, 1);
        n_conditions = 1;
        n_quant_status = 2; % 'Absolute' and 'Relative'
        n_rows = n_trials * n_channels * n_freq_bands * n_conditions * n_quant_status;

        % Precompute conditions array
        conditions_arr = ["TMS"];

        % Initialize the results table
        column_names = {'Trial_index', 'Condition', 'Channel', 'Brain_region', 'Brain_side', 'Freq_band', 'Quant_status', 'PSD'};
        results_table = table('Size', [n_rows, numel(column_names)], 'VariableTypes', {'double', 'categorical', 'string', 'categorical', 'categorical', 'categorical', 'categorical', 'double'}, 'VariableNames', column_names);

        % Parallel loop
        results_parts = cell(1, n_trials);
        parfor trialIdx = 1:n_trials
            temp_results = table('Size', [n_channels * n_freq_bands * n_conditions * n_quant_status, numel(column_names)], 'VariableTypes', {'double', 'categorical', 'string', 'categorical', 'categorical', 'categorical', 'categorical', 'double'}, 'VariableNames', column_names);
            temp_row = 1;

            for conditionIdx = 1:n_conditions
                condition = conditions_arr(conditionIdx);  % Use the precomputed array
                trial_index_adjusted = trialIdx;

                for channelIdx = 1:n_channels
                    channel_label = data.label(channelIdx);
                    brain_region = brain_regions{channelIdx};
                    brain_side = brain_sides{channelIdx};

                    for bandIdx = 1:n_freq_bands
                        band = freq_bands{bandIdx, 1};
                        low_freq = freq_bands{bandIdx, 2};
                        high_freq = freq_bands{bandIdx, 3};

                        freq_analysis = results.(band);
                        % Extract the PSD data and frequency values
                        band_data = freq_analysis.powspctrm(trialIdx, channelIdx, :);
                        freq_values = freq_analysis.freq;

                        % Find indices of the frequency range
                        freq_indices = find(freq_values >= low_freq & freq_values <= high_freq);

                        % Extract the relevant PSD data for bandpower calculation
                        relevant_psd_data = squeeze(band_data(:, freq_indices));
                        freq_range = freq_values(freq_indices);

                        % Calculate integrated power using summing approach
                        integrated_power = sum(relevant_psd_data);

                        % Convert to dB
                        abs_power_db = 10 * log10(integrated_power);

                        % Calculate total power across all bands for relative power
                        total_power = 0;
                        for b = 1:n_freq_bands
                            band_name = freq_bands{b, 1};
                            total_band_analysis = results.(band_name);
                            total_band_data = total_band_analysis.powspctrm(trialIdx, channelIdx, :);

                            % Use summing approach for total power calculation
                            total_power = total_power + sum(total_band_data(:));
                        end

                        % Calculate relative power
                        rel_power = integrated_power / total_power;

                        % Fill the temporary table
                        temp_results(temp_row, :) = {trial_index_adjusted, condition, channel_label, brain_region, brain_side, band, 'Absolute', abs_power_db};
                        temp_row = temp_row + 1;
                        temp_results(temp_row, :) = {trial_index_adjusted, condition, channel_label, brain_region, brain_side, band, 'Relative', rel_power};
                        temp_row = temp_row + 1;
                    end
                end
            end
            results_parts{trialIdx} = temp_results;
        end

        % Concatenate results from all parts
        results_table = vertcat(results_parts{:});

        % Rename the table to power_df
        power_df = results_table;

        % Convert categorical columns to cell arrays of character vectors
        power_df.Condition = cellstr(power_df.Condition);
        power_df.Brain_region = cellstr(power_df.Brain_region);
        power_df.Brain_side = cellstr(power_df.Brain_side);
        power_df.Freq_band = cellstr(power_df.Freq_band);
        power_df.Quant_status = cellstr(power_df.Quant_status);
        power_df.Channel = cellstr(power_df.Channel);

        % Sort the final table by 'Trial_index'
        power_df = sortrows(power_df, 'Trial_index');

        % Convert Table to Cell Array with Headers (easier to load in python)
        headers = power_df.Properties.VariableNames;  % Extract column names
        sensor_power_df = [headers; table2cell(power_df)];  % Combine headers and data

        % Save the power_df as a .mat file
        output_filename = sprintf('%s_sensor_Power_df.mat', sub_num);
        output_path = fullfile(DATAOUT, output_filename);
        save(output_path, 'sensor_power_df');

        % --- Clear large variables to free memory
        clear data brain_region_map ...
            results results_table results_parts power_df ...
            temp_freq band_data freq_values freq_indices relevant_psd_data ...
            freq_range integrated_power abs_power_db total_power rel_power ...
            headers
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HJORTH POWER
    if hjorthpower
        fprintf('\nHjorth power extraction: %s\n', subID)

        data = data_base;

        % --- Filter frequency bands to only include those of interest
        %freq_bands_of_interest = {'delta', 'theta', 'alpha', 'low_beta', 'high_beta', 'low_gamma', 'high_gamma'};
        freq_bands_of_interest = {'theta', 'alpha', 'beta', 'gamma'};
        %freq_bands_of_interest = {'alpha'};
        freq_band_names_all = freq_bands(:, 1);
        selected_idx = ismember(freq_band_names_all, freq_bands_of_interest);
        freq_bands_filtered = freq_bands(selected_idx, :);

        % Initialize results
        hjorth_power_results1 = [];
        n_trials = length(data.trial);  % number of trials
        fsample  = data.fsample;        % sampling rate
        labels   = data.label;          % channel names

        % --- Apply Hjorth spatial filters
        all_trials = cat(3, data.trial{:});  % channels x time x trials
        signal_C3    = spatialFilter(all_trials, labels, channels_C3, weights);
        signal_C4    = spatialFilter(all_trials, labels, channels_C4, weights);
        signal_FCC3h = spatialFilter(all_trials, labels, channels_FCC3h, weights);
        signal_FCC4h = spatialFilter(all_trials, labels, channels_FCC4h, weights);

        % --- Loop through each frequency band
        for bandIdx = 1:size(freq_bands_filtered, 1)
            bandName = freq_bands_filtered{bandIdx, 1};
            freq_range = [freq_bands_filtered{bandIdx, 2}, freq_bands_filtered{bandIdx, 3}];

            % % Use IAF
            % freq_range = [mu_cf_primary(i)-1 mu_cf_primary(i)+1]; 

            % Loop over trials
            for trial = 1:n_trials
                % Compute bandpower for each Hjorth channel
                power_C3    = bandpower(detrend(signal_C3(:, trial)), fsample, freq_range);
                power_C4    = bandpower(detrend(signal_C4(:, trial)), fsample, freq_range);
                power_FCC3h = bandpower(detrend(signal_FCC3h(:, trial)), fsample, freq_range);
                power_FCC4h = bandpower(detrend(signal_FCC4h(:, trial)), fsample, freq_range);

                power_C3_total    = bandpower(detrend(signal_C3(:, trial)), fsample, [4 45]);
                power_C3_relative = power_C3/power_C3_total;        

                % Store results as: {trial, condition, bandName, power values}
                hjorth_power_results1 = [hjorth_power_results1; ...
                    {trial, 'TMS', bandName, power_C3, power_C4, power_FCC3h, power_FCC4h, power_C3_total, power_C3_relative}];
            end
        end

        % --- Convert results to table
        hjorth_power_table1 = cell2table(hjorth_power_results1, ...
            'VariableNames', {'Trial_index', 'Condition', 'Freq_band', ...
            'Power_C3', 'Power_C4', 'Power_FCC3H', 'Power_FCC4H', 'Power_C3_total', 'Power_C3_relative'});

        % % --- Loop through trials first (to get spectra)
        hjorth_power_results2 = [];
        for trial = 1:n_trials
            % Get trial data
            trial_C3    = detrend(signal_C3(:, trial));

            % Compute full PSD for 1/f fitting (using pwelch or periodogram)
            [pxx_C3, f] = pwelch(trial_C3, hamming(fsample), fsample/2, 2^nextpow2(fsample*2), fsample);

            % --- 1/f CORRECTION for C3
            % Keep only 2–45 Hz for fitting range
            mask_full = (f >= 2) & (f <= 45);
            f_fit = f(mask_full);
            pxx_C3_fit = pxx_C3(mask_full);

            % Fit on "background" frequencies (avoid alpha 8-13, beta 13-30)
            mask_1f = ((f_fit >= 2 & f_fit < 8) | (f_fit > 30 & f_fit <= 45));

            % --- C3 correction
            x = log10(f_fit(mask_1f));
            y_C3 = log10(pxx_C3_fit(mask_1f));
            p_C3 = polyfit(x, y_C3, 1);
            yfit_log_C3 = polyval(p_C3, log10(f_fit));
            yfit_lin_C3 = 10.^yfit_log_C3;
            pxx_C3_corrected = pxx_C3_fit ./ yfit_lin_C3;

            % --- Now compute bandpower from CORRECTED spectra
            for bandIdx = 1:size(freq_bands_filtered, 1)
                bandName = freq_bands_filtered{bandIdx, 1};
                freq_range = [freq_bands_filtered{bandIdx, 2}, freq_bands_filtered{bandIdx, 3}];

                % Extract band from corrected spectrum
                band_mask = (f_fit >= freq_range(1)) & (f_fit <= freq_range(2));

                % Integrate power in band (corrected)
                power_C3_corrected = trapz(f_fit(band_mask), pxx_C3_corrected(band_mask));

                % Also compute UNCORRECTED for comparison
                power_C3_uncorrected = trapz(f_fit(band_mask), pxx_C3_fit(band_mask));

                % Total power (4-45 Hz, corrected)
                total_mask = (f_fit >= 4) & (f_fit <= 45);
                power_C3_total_corrected = trapz(f_fit(total_mask), pxx_C3_corrected(total_mask));
                power_C3_relative_corrected = power_C3_corrected / power_C3_total_corrected;

                % Store: trial, condition, band, powers (uncorrected), powers (corrected), total, relative, PSDs, frequencies, 1/f fits
                hjorth_power_results2 = [hjorth_power_results2; ...
                    {trial, 'TMS', bandName, power_C3_uncorrected, power_C3_corrected, ...
                    power_C3_total_corrected, power_C3_relative_corrected, f, pxx_C3, f_fit, pxx_C3_fit, pxx_C3_corrected, ...
                    yfit_lin_C3}];
            end
        end

        % --- Convert results to table
        hjorth_power_table2 = cell2table(hjorth_power_results2, ...
            'VariableNames', {'Trial_index', 'Condition', 'Freq_band', ...
            'Power_C3_uncorrected', 'Power_C3_corrected', 'Power_C3_total_corrected', 'Power_C3_relative_corrected', ...
            'f_full', 'PSD_C3_full', 'f_fit', 'PSD_C3_fit', 'PSD_C3_corrected', 'fit_1f_C3'});

        % --- Concatenate: keep first 3 columns from table1, add columns 4+ from table2
        hjorth_power_table = [hjorth_power_table1, hjorth_power_table2(:, 4:end)];

        % --- Convert Table to Cell Array with Headers (for Python loading)
        headers = hjorth_power_table.Properties.VariableNames;
        hjorth_power_table = [headers; table2cell(hjorth_power_table)];
       
        % --- Save results
        output_filename = sprintf('%s_sensor_Hjorth_Power_df.mat', sub_num);
        %output_filename = sprintf('%s_sensor_Hjorth_Power_IAF_df.mat', sub_num);
        output_path = fullfile(DATAOUT, output_filename);
        save(output_path, 'hjorth_power_table');

        % --- Clear large variables to free memory
        clear data all_trials signal_C3 signal_C4 signal_FCC3h signal_FCC4h ...
            power_C3 power_C4 power_FCC3h power_FCC4h hjorth_power_results power_C3_uncorrected power_C3_corrected ...
            power_C3_total_corrected power_C3_relative_corrected f pxx_C3 f_fit pxx_C3_fit pxx_C3_corrected ...
            yfit_lin_C3
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HJORTH CONNECTIVITY
    if hjorthconnectivity
        fprintf('\nHjorth connectivity extraction: %s\n', subID)

        data = data_base;

        % --- Filter frequency bands to only include those of interest
        %freq_bands_of_interest = {'delta', 'theta', 'alpha', 'low_beta', 'high_beta', 'low_gamma', 'high_gamma'};
        freq_bands_of_interest = {'theta', 'alpha', 'beta', 'gamma'};
        %freq_bands_of_interest = {'alpha'};
        freq_band_names_all = freq_bands(:, 1);
        selected_idx = ismember(freq_band_names_all, freq_bands_of_interest);
        freq_bands_filtered = freq_bands(selected_idx, :);

        % Apply Hjorth spatial filters
        labels   = data.label;
        all_trials = cat(3, data.trial{:});  % channels x time x trials
        signal_C3    = spatialFilter(all_trials, labels, channels_C3, weights);
        signal_C4    = spatialFilter(all_trials, labels, channels_C4, weights);
        signal_FCC3h = spatialFilter(all_trials, labels, channels_FCC3h, weights);
        signal_FCC4h = spatialFilter(all_trials, labels, channels_FCC4h, weights);

        signal_labels = {'C3', 'C4', 'FCC3h', 'FCC4h'};
        signal_combinations = nchoosek(1:length(signal_labels), 2); % All pairwise combinations of signals
        indxcmb = nchoosek(1:length(signal_labels), 2); % All pairwise combinations of signals
        fsample = epochs.fsample; % sampling frequency in Hz
        ARord = 30; % AR model order for padding
        numTrials = length(data.trial);
      
        % Preallocate arrays for connectivity results
        iPLV_results = cell(length(freq_bands_of_interest), 1);
        %PLI_results = cell(length(freq_bands_names), 1);
        wPLI_results = cell(length(freq_bands_of_interest), 1);
        Coh_results = cell(length(freq_bands_of_interest), 1);
        imCoh_results = cell(length(freq_bands_of_interest), 1);
        %lagCoh_results = cell(length(freq_bands_names), 1);
        %oPEC_results = cell(length(freq_bands_names), 1);

        % --- Loop over selected frequency bands
        for bandIdx = 1:size(freq_bands_filtered, 1)
            bandName = freq_bands_filtered{bandIdx, 1};
            bpfreq = [freq_bands_filtered{bandIdx, 2}, freq_bands_filtered{bandIdx, 3}];
            
            % % Use IAF
            % if strcmpi(bandName, 'alpha')
            %     bpfreq = [mu_cf_primary(i)-1 mu_cf_primary(i)+1];
            % end
             
            n = round((0.15*fsample)/2)*2; % filter order
            b = fir1(n, bpfreq/(fsample/2)); % filter coefficients
            padsize = round(2*fsample/mean(bpfreq)); % padding size

            % Initialize arrays for this band
            iPLV_band = zeros(size(indxcmb, 1), numTrials);
            %PLI_band = zeros(size(indxcmb, 1), numTrials);
            wPLI_band = zeros(size(indxcmb, 1), numTrials);
            Coh_band = zeros(size(indxcmb, 1), numTrials);
            imCoh_band = zeros(size(indxcmb, 1), numTrials);
            %lagCoh_band = zeros(size(indxcmb, 1), numTrials);
            %oPEC_band = zeros(size(indxcmb, 1), numTrials);

            % Loop over trials and signal combinations
            for trial = 1:numTrials
                % Extract trial data from each spatially filtered signal
                trial_data_C3 = signal_C3(:, trial);
                trial_data_C4 = signal_C4(:, trial);
                trial_data_FCC3h = signal_FCC3h(:, trial);
                trial_data_FCC4h = signal_FCC4h(:, trial);

                % Define the signals and their labels
                signals = {signal_C3, signal_C4, signal_FCC3h, signal_FCC4h};

                for combIdx = 1:size(signal_combinations, 1)
                    signal1 = signals{signal_combinations(combIdx, 1)}(:, trial); %has to be 1001 x 1
                    signal2 = signals{signal_combinations(combIdx, 2)}(:, trial);

                    % Preprocess the data
                    preprocessedData = preprocessData([signal1, signal2], b, ARord, padsize);

                    % Calculate connectivity measures
                    [iPLV, ~, wPLI, Coh, imCoh, ~, ~] = calculateConnectivity(preprocessedData, padsize, [1, 2]);

                    % Store results for this trial and band
                    iPLV_band(combIdx, trial) = iPLV;
                    %PLI_band(combIdx, trial) = PLI;
                    wPLI_band(combIdx, trial) = wPLI;
                    Coh_band(combIdx, trial) = Coh;
                    imCoh_band(combIdx, trial) = imCoh;
                    %lagCoh_band(combIdx, trial) = lagCoh;
                    %oPEC_band(combIdx, trial) = oPEC;
                end
            end
            % Store results for this frequency band
            iPLV_results{bandIdx} = iPLV_band;
            %PLI_results{bandIdx} = PLI_band;
            wPLI_results{bandIdx} = wPLI_band;
            Coh_results{bandIdx} = Coh_band;
            imCoh_results{bandIdx} = imCoh_band;
            %lagCoh_results{bandIdx} = lagCoh_band;
            %oPEC_results{bandIdx} = oPEC_band;
        end

        % Create source connectivity dataframe
        % Start a parallel pool (if not already started)
        if isempty(gcp('nocreate'))
            parpool; % Start with default settings
        end

        % Define signal pair labels
        signalPairLabels = arrayfun(@(x) sprintf('%s-%s', signal_labels{signal_combinations(x, 1)}, signal_labels{signal_combinations(x, 2)}), 1:size(signal_combinations, 1), 'UniformOutput', false);

        % Preallocate cell array for parallel results
        parallelResults = cell(length(freq_bands_filtered), 1);
        %parallelResults = cell(1, 1);

        % Loop through each frequency band in parallel
        parfor bandIdx = 1:length(freq_bands_filtered)
            bandTable = table();
            for trialIdx = 1:numTrials
                for signalPairIdx = 1:length(signalPairLabels)
                    % newRow = {trialIdx, conditions{trialIdx}, freq_bands_names{bandIdx}, signalPairLabels{signalPairIdx}, ...
                    %           iPLV_results{bandIdx}(signalPairIdx, trialIdx), PLI_results{bandIdx}(signalPairIdx, trialIdx), ...
                    %           wPLI_results{bandIdx}(signalPairIdx, trialIdx), Coh_results{bandIdx}(signalPairIdx, trialIdx), ...
                    %           imCoh_results{bandIdx}(signalPairIdx, trialIdx), lagCoh_results{bandIdx}(signalPairIdx, trialIdx), ...
                    %           oPEC_results{bandIdx}(signalPairIdx, trialIdx)};
                    newRow = {trialIdx, 'TMS', freq_bands_filtered{bandIdx}, signalPairLabels{signalPairIdx}, ...
                        iPLV_results{bandIdx}(signalPairIdx, trialIdx), ...
                        wPLI_results{bandIdx}(signalPairIdx, trialIdx), Coh_results{bandIdx}(signalPairIdx, trialIdx), ...
                        imCoh_results{bandIdx}(signalPairIdx, trialIdx)};
                    bandTable = [bandTable; newRow];
                end
            end
            parallelResults{bandIdx} = bandTable;
        end

        % Concatenate all tables
        hjorth_connectivity_df = vertcat(parallelResults{:});
        % hjorthconnectivity_df.Properties.VariableNames = {'Trial_index', 'Condition', 'Freq_band', 'Signal_comb', ...
        %                                                   'iPLV', 'PLI', 'wPLI', 'Coh', 'imCoh', 'lagCoh', 'oPEC'};
        hjorth_connectivity_df.Properties.VariableNames = {'Trial_index', 'Condition', 'Freq_band', 'Signal_comb', ...
            'iPLV', 'wPLI', 'Coh', 'imCoh'};

        % Convert Table to Cell Array with Headers (easier to load in python)
        headers = hjorth_connectivity_df.Properties.VariableNames;  % Extract column names
        hjorth_connectivity_df = [headers; table2cell(hjorth_connectivity_df)];  % Combine headers and data

        % --- Save
        output_filename = sprintf('%s_sensor_Hjorth_Connectivity_df.mat', sub_num);
        %output_filename = sprintf('%s_sensor_Hjorth_Connectivity_IAF_df.mat', sub_num);
        output_path = fullfile(DATAOUT, output_filename);
        save(output_path, 'hjorth_connectivity_df');

        % --- Clear large variables
        clear data all_trials signal_C3 signal_C4 signal_FCC3h signal_FCC4h ...
            iPLV_results wPLI_results Coh_results imCoh_results ...
            parallelResults bandTable preprocessedData ...
            iPLV_band wPLI_band Coh_band imCoh_band
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HJORTH PHASE
    if hjorthphase
        fprintf('\nHjorth mu phase extraction: %s\n', subID)

        data = data_base;

        % Define parameters for phase estimation
        arOrder = 30;
        edge = 65;
        hilbertwindow = 128;
        offset_correction = 44;
        sampling_rate = 1000;
        filter_order = round(sampling_rate/5);

        % Initialize filters
        freq_idx = find(strcmp(freq_bands(:,1), phase_of_interest));
        bpfreq  = [freq_bands{freq_idx,2}, freq_bands{freq_idx,3}];

        % % Use IAF
        % bpfreq = [mu_cf_primary(i)-1 mu_cf_primary(i)+1]; 

        D_bandpass = designfilt('bandpassfir', 'FilterOrder', filter_order, 'CutoffFrequency1', bpfreq(1) , 'CutoffFrequency2', bpfreq(2), 'SampleRate', sampling_rate, 'DesignMethod', 'window');
        %D_nyquist = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 0.5, 'SampleRate', sampling_rate);

        % --- Initialize result arrays        
        trial_indices = [];
        condition_list = {};
        phase_C3_list = [];
        phase_C4_list = [];
        phase_FCC3h_list = [];
        phase_FCC4h_list = [];

        % --- Loop over trials
        n_trials = length(data.trial);
        labels   = data.label;

        all_trials = cat(3, data.trial{:});  % channels x time x trials
        % Apply Hjorth spatial filters
        signal_C3    = spatialFilter(all_trials, labels, channels_C3, weights);
        signal_C4    = spatialFilter(all_trials, labels, channels_C4, weights);
        signal_FCC3h = spatialFilter(all_trials, labels, channels_FCC3h, weights);
        signal_FCC4h = spatialFilter(all_trials, labels, channels_FCC4h, weights);

        for trial = 1:n_trials
            % Extract trial data from each spatially filtered signal
            trial_data_C3 = signal_C3(:, trial);
            trial_data_C4 = signal_C4(:, trial);
            trial_data_FCC3h = signal_FCC3h(:, trial);
            trial_data_FCC4h = signal_FCC4h(:, trial);

            % Initialize lists to store phase for each signal
            phase_list = cell(1, 4);

            % Detrend, filter, and phase estimation for each signal
            for i = 1:4
                % Select the appropriate trial data
                switch i
                    case 1
                        current_signal = trial_data_C3;
                    case 2
                        current_signal = trial_data_C4;
                    case 3
                        current_signal = trial_data_FCC3h;
                    case 4
                        current_signal = trial_data_FCC4h;
                end

                % Detrend
                detrended_signal = detrend(current_signal);

                % Downsample to 500 Hz
                %downsampled_signal = detrended_signal(1:2:end);

                % Filter
                %filtered_signal = filter(D_nyquist, downsampled_signal);
                filtered_signal = filtfilt(D_bandpass, detrended_signal); %note that filtfilt uses reflection and sets the initial values

                % Phase estimation
                %[phase, ~] = phastimate(current_signal, D_bandpass, edge, arOrder, hilbertwindow, offset_correction);
                [phase, ~] = phastimate_adapted(filtered_signal, edge, arOrder, hilbertwindow, offset_correction);

                % Store the last phase value
                phase_list{i} = phase(end);
            end

            % Append trial results
            phase_C3_list    = [phase_C3_list; phase_list{1}];
            phase_C4_list    = [phase_C4_list; phase_list{2}];
            phase_FCC3h_list = [phase_FCC3h_list; phase_list{3}];
            phase_FCC4h_list = [phase_FCC4h_list; phase_list{4}];
            trial_indices    = [trial_indices; trial];
            condition_list   = [condition_list; {'TMS'}];  % all labeled TMS
        end
        %polarhistogram(phase_C3_list, 30)

        % --- Convert phase to sine/cosine representation
        Sin_Phase_C3H    = sin(phase_C3_list);
        Cos_Phase_C3H    = cos(phase_C3_list);
        Sin_Phase_C4H    = sin(phase_C4_list);
        Cos_Phase_C4H    = cos(phase_C4_list);
        Sin_Phase_FCC3H  = sin(phase_FCC3h_list);
        Cos_Phase_FCC3H  = cos(phase_FCC3h_list);
        Sin_Phase_FCC4H  = sin(phase_FCC4h_list);
        Cos_Phase_FCC4H  = cos(phase_FCC4h_list);

        % --- Build results table
        hjorth_phase_df = table(trial_indices, condition_list, ...
            phase_C3_list, phase_C4_list, phase_FCC3h_list, phase_FCC4h_list, ...
            Sin_Phase_C3H, Cos_Phase_C3H, ...
            Sin_Phase_C4H, Cos_Phase_C4H, ...
            Sin_Phase_FCC3H, Cos_Phase_FCC3H, ...
            Sin_Phase_FCC4H, Cos_Phase_FCC4H, ...
            'VariableNames', {'Trial_index','Condition', ...
            'Phase_C3H','Phase_C4H','Phase_FCC3H','Phase_FCC4H', ...
            'Sin_Phase_C3H','Cos_Phase_C3H', ...
            'Sin_Phase_C4H','Cos_Phase_C4H', ...
            'Sin_Phase_FCC3H','Cos_Phase_FCC3H', ...
            'Sin_Phase_FCC4H','Cos_Phase_FCC4H'});

        % --- Convert to cell array with headers (for Python loading)
        headers = hjorth_phase_df.Properties.VariableNames;
        hjorth_phase_df = [headers; table2cell(hjorth_phase_df)];

        % --- Save
        %output_filename = sprintf('%s_sensor_Hjorth_Phase_df.mat', sub_num);
        output_filename = sprintf('%s_sensor_Hjorth_Phase_IAF_df.mat', sub_num);
        output_path = fullfile(DATAOUT, output_filename);
        save(output_path, 'hjorth_phase_df');

        % --- Clear large variables
        clear trial_data signals sig_C3 sig_C4 sig_FCC3h sig_FCC4h ...
            phase_list detr_sig down_sig filt_sig phase ...
            Sin_Phase_C3H Cos_Phase_C3H Sin_Phase_C4H Cos_Phase_C4H ...
            Sin_Phase_FCC3H Cos_Phase_FCC3H Sin_Phase_FCC4H Cos_Phase_FCC4H
    end
end
time = toc

