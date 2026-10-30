%% T2TVAR BRAIN-STATE-RESOLVED VARIABILITY QUANTIFICATION
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

prestim_feature     = 'Phase'; %Power %Connectivity %Phase %PowerPhase %PowerConnectivity %Phase
bands_of_interest   = {'alpha'}; %{'theta','alpha','beta', 'gamma'};
if strcmpi(prestim_feature, 'Phase') || strcmpi(prestim_feature, 'PowerPhase') || strcmpi(prestim_feature, 'PowerConnectivity')
    bands_of_interest = {'alpha'};
end
splits              = '4phasebins'; %tails25 %tertiles %median %quintiles %octiles %deciles '4phasebins' '8phasebins'
nPhaseBins          = 4;
if strcmpi(prestim_feature, 'Phase') && nPhaseBins == 4
    splits          = '4phasebins';
elseif strcmpi(prestim_feature, 'Phase') && nPhaseBins == 8
    splits          = '8phasebins'; 
end
Phase_split = false;
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

individual_subjects = true;

save_output         = false;

generate_plots      = true;

plots_show          = 'VARrel'; %VARrel %VAR %MEAN_AMP
if strcmpi(plots_show, 'MEAN_AMP')
    MEP_feature =  'p2p_zscored_log';
    TEP_metric_type = 'z_log';
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

% % REVISION BRAIN STIMULATION RECONSTRUCTED CHANNELS POWER
% load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\reftep_C3_affected_subjects_list.mat');
% alldata_MEP_master = alldata_MEP_master(1:29);
% alldata_MEP_master = alldata_MEP_master(~ismember({alldata_MEP_master.subID}, C3_affected_ids));
% all_subjects = {alldata_MEP_master.subID};
% sub_range = 1:17;

% % REVISION BRAIN STIMULATION WITHOUT NON-RESPONSE REJECTION
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_without_nonresponse_rejection.mat'), 'alldata_MEP_master');
 
% % REVISION BRAIN STIMULATION WITH OUTLIER REJECTION
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_with_outlier_rejection.mat'), 'alldata_MEP_master');

% REVISION BRAIN STIMULATION APB FDI separately
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_APB.mat'), 'alldata_MEP_master');

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

    % --------- Variability quantification ---------
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')

        for fb = 1:numel(bands_of_interest)
            thisBand = bands_of_interest{fb};

            % Load splits table of band
            Tb = splitsTbl(strcmpi(splitsTbl.Freq_band, thisBand), :);
            if isempty(Tb)
                return
            end

            % Get indices for High/Low trials
            high_idx  = double(Tb.Trial_index(Tb.Index == 1));
            low_idx  = double(Tb.Trial_index(Tb.Index == 3));

            if fb == 1
                if ~strcmpi(var_method, 'TTV_NIU')
                    rand_idx = randperm(numel(all_data), numel(low_idx));
                elseif strcmpi(var_method, 'TTV_NIU')
                    rand_idx = randperm(size(all_data,1), numel(low_idx));
                end

                if do_bootstrap
                    n_rand = numel(low_idx);
                    rand_idx = cell(n_boot_iter,1);
                    if ~strcmpi(var_method, 'TTV_NIU')
                        n_total = numel(all_data);
                    else
                        n_total = size(all_data,1);
                    end
                    for i = 1:n_boot_iter
                        rand_idx{i} = randperm(n_total, n_rand);
                    end
                end
            end   

            if numel(high_idx) < nPerState_general || numel(low_idx) < nPerState_general
                disp(['Subject skipped. Not enough trials for ' subID])
                continue
            end

            % Compute mean per state
            MEAN_all = mean(all_data);
            MEAN_high = mean(all_data(high_idx));
            MEAN_low = mean(all_data(low_idx));            
            MEDIAN_all = median(all_data);
            MEDIAN_high = median(all_data(high_idx));
            MEDIAN_low = median(all_data(low_idx));

            if do_bootstrap
                n_iter = n_boot_iter;
                mean_rand_iter   = zeros(n_iter,1);
                median_rand_iter = zeros(n_iter,1);
                for i = 1:n_iter
                    idx = rand_idx{i,:};
                    mean_rand_iter(i)   = mean(all_data(idx));
                    median_rand_iter(i) = median(all_data(idx));
                end
                MEAN_rand   = mean(mean_rand_iter);
                MEDIAN_rand = mean(median_rand_iter);
            else
                MEAN_rand = mean(all_data(rand_idx));
                MEDIAN_rand = median(all_data(rand_idx));
            end

            % Compute variability per state
            switch upper(var_method)
                case 'SD' % Standard deviation (assumes log-microvolts)
                    VAR_all = std(all_data);
                    VAR_high = std(all_data(high_idx));
                    VAR_low = std(all_data(low_idx));

                    if do_bootstrap
                        n_iter = size(rand_idx,1);
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = rand_idx{i,:};
                            var_rand_iter(i)   = std(all_data(idx));
                        end
                        VAR_rand   = mean(var_rand_iter);
                    else
                        VAR_rand = std(all_data(rand_idx));
                    end
                    
                    % SD of mean-normalize data
                    if bin_mean_normalization
                        VAR_all = std((all_data / MEAN_all));
                        VAR_high   = std((all_data(high_idx) / MEAN_high));
                        VAR_low   = std((all_data(low_idx) / MEAN_low));
                        VAR_rand = std((all_data(rand_idx) / MEAN_rand));
                    end
                case 'CV' % Coefficient of variation (assumes raw microvolts)
                    VAR_all  = std(all_data) / mean(all_data);
                    VAR_high   = std(all_data(high_idx)) / mean(all_data(high_idx));
                    VAR_low   = std(all_data(low_idx)) / mean(all_data(low_idx));

                    if do_bootstrap
                        n_iter = size(rand_idx,1);
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = rand_idx{i,:};
                            var_rand_iter(i)   = std(all_data(idx)) / mean(all_data(idx));
                        end
                        VAR_rand   = mean(var_rand_iter);
                    else
                        VAR_rand = std(all_data(rand_idx)) / mean(all_data(rand_idx));
                    end
                    
                case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                    VAR_all  = 1.4826 * mad(all_data, 1);
                    VAR_high   = 1.4826 * mad(all_data(high_idx), 1);
                    VAR_low   = 1.4826 * mad(all_data(low_idx), 1);
                    VAR_rand = 1.4826 * mad(all_data(rand_idx), 1);
                case 'IQR' % Interquartile range (assumes log-microvolts)
                    VAR_all  = iqr(all_data);
                    VAR_high   = iqr(all_data(high_idx));
                    VAR_low   = iqr(all_data(low_idx));
                    VAR_rand = iqr(all_data(rand_idx));

                    % SD of mean-normalize data
                    if bin_mean_normalization
                        VAR_all  = iqr(all_data / MEAN_all);
                        VAR_high   = iqr(all_data(high_idx) / MEAN_high);
                        VAR_low   = iqr(all_data(low_idx) / MEAN_low);
                        VAR_rand = iqr(all_data(rand_idx) / MEAN_rand);
                    end
                case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)
                    VAR_all  = calculate_CQV(all_data);
                    VAR_high = calculate_CQV(all_data(high_idx));
                    VAR_low  = calculate_CQV(all_data(low_idx));
                    VAR_rand = calculate_CQV(all_data(rand_idx));
                case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                    VAR_all = calculate_TTV_Niu(all_data, nSurrogate);
                    VAR_high = calculate_TTV_Niu(all_data(high_idx,:), nSurrogate);
                    VAR_low = calculate_TTV_Niu(all_data(low_idx,:), nSurrogate);
                    VAR_rand = calculate_TTV_Niu(all_data(rand_idx,:), nSurrogate);
                otherwise
                    error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
            end

            % Compute relative variability
            VARrel_high = (VAR_high   / VAR_all) * 100;
            VARrel_low  = (VAR_low   / VAR_all) * 100;
            VARrel_rand = (VAR_rand / VAR_all) * 100;

            % Compute ratios relative to Random
            ratio_high = VAR_high / VAR_rand;
            ratio_low  = VAR_low / VAR_rand;
            logratio_high = log(ratio_high);
            logratio_low  = log(ratio_low);

            % % Compute deltas
            DeltaSpread_highrand    = VARrel_high - VARrel_rand;
            DeltaSpread_lowrand     = VARrel_low - VARrel_rand;
            DeltaSpread_highlow = VARrel_high  - VARrel_low;

            % Save results to table
            count = count + 1;
            resultsTable(count).subject = subID;
            resultsTable(count).site = site;
            resultsTable(count).type    = horzcat('Hjorth_', prestim_feature);
            resultsTable(count).band    = thisBand;
            resultsTable(count).MEAN_all       = MEAN_all;
            resultsTable(count).MEAN_high       = MEAN_high;
            resultsTable(count).MEAN_low        = MEAN_low;
            resultsTable(count).MEAN_rand       = MEAN_rand;
            resultsTable(count).VAR_method      = var_method;
            resultsTable(count).VAR_all         = VAR_all;
            resultsTable(count).VAR_high        = VAR_high;
            resultsTable(count).VAR_low         = VAR_low;
            resultsTable(count).VAR_rand        = VAR_rand;
            resultsTable(count).VARrel_high       = VARrel_high;
            resultsTable(count).VARrel_low        = VARrel_low;
            resultsTable(count).VARrel_rand       = VARrel_rand;
            resultsTable(count).ratio_high    = ratio_high;
            resultsTable(count).ratio_low     = ratio_low;
            resultsTable(count).logratio_high = logratio_high;
            resultsTable(count).logratio_low  = logratio_low;
            resultsTable(count).DeltaSpread_highrand  = DeltaSpread_highrand;
            resultsTable(count).DeltaSpread_lowrand   = DeltaSpread_lowrand;
            resultsTable(count).DeltaSpread_highlow = DeltaSpread_highlow;
        end

    elseif strcmpi(prestim_feature, 'Phase')
        Tb = splitsTbl; % ALPHA ONLY

        trialsPerBin = arrayfun(@(b) sum(Tb.Phase_Bin == b), unique(Tb.Phase_Bin));
        [minTrials, minTrials_idx] = min(trialsPerBin);
        if minTrials < nPerState_phase
            fprintf('%s does not have enough trials for bin %d (n = %d < %d)\n', subID, minTrials_idx, minTrials, nPerState_phase);
            disp(trialsPerBin')
            continue
        end

        VAR_phasebins = nan(1, nPhaseBins);
        MEAN_phasebins = nan(1, nPhaseBins);
        for b = 1:nPhaseBins

            % Get phase bins indices
            idx_thisBin = double(Tb.Trial_index(Tb.Phase_Bin == b));
            nTrials_thisBin = numel(idx_thisBin);

            if do_bootstrap
                n_rand = minTrials;
                idx_thisBin_thisSub = cell(n_boot_iter,1);
                for i = 1:n_boot_iter
                    rand_idx = randperm(nTrials_thisBin, n_rand);
                    idx_thisBin_thisSub{i} = idx_thisBin(rand_idx);
                end                
            else
                rand_idx = randperm(nTrials_thisBin, minTrials);
                idx_thisBin_thisSub = idx_thisBin(rand_idx);
                %idx_thisBin_thisSub = idx_thisBin;
            end

            % Compute mean per phase bin
            MEAN_all = mean(all_data);

            if do_bootstrap
                n_iter = n_boot_iter;
                mean_rand_iter   = zeros(n_iter,1);
                for i = 1:n_iter
                    idx = idx_thisBin_thisSub{i};
                    mean_rand_iter(i)   = mean(all_data(idx));
                end
                MEAN_phasebins(b)   = mean(mean_rand_iter);
            else
                MEAN_phasebins(b) = mean(all_data(idx_thisBin_thisSub));  % mean amplitude within bin
            end

            % Compute variability per phase bin
            switch upper(var_method)
                case 'SD' % Standard deviation (assumes log-microvolts)
                    VAR_all  = std(all_data);
                    
                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_thisBin_thisSub{i};
                            var_rand_iter(i)   = std(all_data(idx));
                        end
                        VAR_phasebins(b)  = mean(var_rand_iter);
                    else
                        VAR_phasebins(b)  = std(all_data(idx_thisBin_thisSub));   % spread within bin
                    end
                    
                    % SD of mean-normalize data
                    if bin_mean_normalization
                        VAR_all  = std(all_data / MEAN_all);
                        VAR_phasebins(b)  = std((all_data(idx_thisBin_thisSub) / MEAN_phasebins(b)));
                    end
                case 'CV' % Coefficient of variation (assumes raw microvolts)
                    VAR_all  = std(all_data) / mean(all_data);
                   
                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_thisBin_thisSub{i};
                            var_rand_iter(i)   = std(all_data(idx)) / MEAN_phasebins(b);
                        end
                        VAR_phasebins(b)  = mean(var_rand_iter);
                    else
                        VAR_phasebins(b)  = std(all_data(idx_thisBin_thisSub)) / MEAN_phasebins(b);
                    end

                case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                    VAR_all  = 1.4826 * mad(all_data, 1);
                    VAR_phasebins(b)  = 1.4826 * mad(all_data(idx_thisBin_thisSub), 1);
                case 'IQR' % Interquartile range (assumes log-microvolts)
                    VAR_all  = iqr(all_data);

                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_thisBin_thisSub{i};
                            var_rand_iter(i)   = iqr(all_data(idx));
                        end
                        VAR_phasebins(b)  = mean(var_rand_iter);
                    else
                        VAR_phasebins(b) = iqr(all_data(idx_thisBin_thisSub));
                    end
                    
                case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)
                    VAR_all  = calculate_CQV(all_data);

                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_thisBin_thisSub{i};
                            var_rand_iter(i)   = calculate_CQV(all_data(idx));
                        end
                        VAR_phasebins(b)  = mean(var_rand_iter);
                    else
                        VAR_phasebins(b) = calculate_CQV(all_data(idx_thisBin_thisSub));
                    end
                case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                    VAR_all = calculate_TTV_Niu(all_data, nSurrogate);
                    VAR_phasebins(b) = calculate_TTV_Niu(all_data(idx_thisBin_thisSub,:), nSurrogate);
                otherwise
                    error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
            end
        end

        % Add Random phases bin
        idx_allBins = double(splitsTbl.Trial_index(splitsTbl.Phase_Bin > 0));
        nTrials_allBins = numel(idx_allBins);

        if do_bootstrap
            n_rand = minTrials;
            idx_rand = cell(n_boot_iter,1);
            for i = 1:n_boot_iter
                rand_idx = randperm(nTrials_allBins, n_rand);
                idx_rand{i} = idx_allBins(rand_idx);
            end
        else
            rand_idx = randperm(nTrials_allBins, minTrials);
            idx_rand = idx_allBins(rand_idx);
        end

        if do_bootstrap
            n_iter = n_boot_iter;
            mean_rand_iter   = zeros(n_iter,1);
            for i = 1:n_iter
                idx = idx_rand{i};
                mean_rand_iter(i)   = mean(all_data(idx));
            end
            MEAN_rand   = mean(mean_rand_iter);
        else
            MEAN_rand = mean(all_data(idx_rand));
        end

        switch upper(var_method)
            case 'SD' % Standard deviation (assumes log-microvolts)
                
                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = std(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand  = std(all_data(idx_rand));
                end
                
                
                % SD of mean-normalize data
                if bin_mean_normalization
                    VAR_rand  = std((all_data(idx_rand) / MEAN_rand));
                end                
            case 'CV' % Coefficient of variation (assumes raw microvolts)

                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = std(all_data(idx)) / mean(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand = std(all_data(idx_rand)) / mean(all_data(idx_rand));
                end

            case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                VAR_rand = 1.4826 * mad(all_data(idx_rand), 1);
            case 'IQR' % Interquartile range (assumes log-microvolts)

                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = iqr(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand = iqr(all_data(idx_rand));
                end

            case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)

                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = calculate_CQV(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand = calculate_CQV(all_data(idx_rand));
                end

            case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                VAR_rand = calculate_TTV_Niu(all_data(idx_rand,:), nSurrogate);
            otherwise
                error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
        end

        % Compute relative variability
        VARrel_phasebins = (VAR_phasebins / VAR_all) * 100;
        VARrel_rand = (VAR_rand / VAR_all) * 100;

        % Save results to table
        count = count + 1;
        resultsTable(count).subject = subID;
        resultsTable(count).site = site;
        resultsTable(count).type    = horzcat('Hjorth_', prestim_feature);
        resultsTable(count).band    = 'alpha';
        resultsTable(count).MEAN_all       = MEAN_all;
        resultsTable(count).MEAN_phasebins = MEAN_phasebins;
        resultsTable(count).MEAN_rand = MEAN_rand;
        resultsTable(count).VAR_method      = var_method;
        resultsTable(count).VAR_all         = VAR_all;
        resultsTable(count).VAR_phasebins   = VAR_phasebins;
        resultsTable(count).VAR_rand   = VAR_rand;
        resultsTable(count).VARrel_phasebins       = VARrel_phasebins;
        resultsTable(count).VARrel_rand       = VARrel_rand;
        resultsTable(count).Phase_C3H       = Tb.Phase_C3H;

    elseif strcmpi(prestim_feature, 'PowerPhase') % ALPHA!!!

        splits_power = splits;
        splitsTbl_Power = load_splits_table(DATAIN_SPLIT, sub_num, 1, 'Power', splits);
        Tb_Power =  splitsTbl_Power(strcmpi(splitsTbl_Power.Freq_band, 'alpha'), :);

        if strcmpi(prestim_feature, 'PowerPhase') && nPhaseBins == 4
            splits          = '4phasebins';
        elseif strcmpi(prestim_feature, 'PowerPhase') && nPhaseBins == 8
            splits          = '8phasebins';
        end
        Tb_Phase = load_splits_table(DATAIN_SPLIT, sub_num, 1, 'Phase', splits);
        splits = splits_power; 

        % Check coherency between both features
        if isempty(Tb_Power) || isempty(Tb_Phase)
            warning('Splits not found for %s', subID);
            return
        end
        if height(Tb_Power) ~= height(Tb_Phase)
            warning('Mismatch in trials')
            return
        end
        assert(isequal(Tb_Power.Trial_index, Tb_Phase.Trial_index), 'Trial mismatch between Power and Phase');

        % Combine tables of both features
        phaseVars = Tb_Phase.Properties.VariableNames;
        phaseVars = setdiff(phaseVars, Tb_Power.Properties.VariableNames);
        Tb_PowerPhase = [Tb_Power Tb_Phase(:, phaseVars)];

        % Get indices for each feature combination
        powers = {'low', 'high'};
        phases = {'Trough', 'Rising', 'Peak', 'Falling'};
        powStr   = string(Tb_PowerPhase.Label);
        phaseStr = string(Tb_PowerPhase.Phase_Label);

        Tb_PowerPhase.Combo = strings(height(Tb_PowerPhase),1);
        idx_combos = struct();
        for iP = 1:numel(powers)
            for iPh = 1:numel(phases)
                comboName = sprintf('%s_%s', powers{iP}, phases{iPh});
                idx_thisCombo = strcmpi(powStr,   powers{iP}) & ...
                    strcmpi(phaseStr, phases{iPh});
                idx_combos.(comboName) = find(idx_thisCombo);
                Tb_PowerPhase.Combo(idx_thisCombo) = comboName;
            end
        end
        Tb_PowerPhase.Combo = categorical(Tb_PowerPhase.Combo);
        counts = countcats(Tb_PowerPhase.Combo);

        % Pseudo N-match all conditions
        comboNames = fieldnames(idx_combos);
        nCombos    = numel(comboNames);
        trialsPerCombo = cellfun(@(name) numel(idx_combos.(name)), comboNames);
        [minTrials, minTrials_idx] = min(trialsPerCombo);
        if minTrials < nPerState_combos
            fprintf('%s does not have enough trials for combo %s (n = %d < %d)\n', ...
                subID, comboNames{minTrials_idx}, minTrials, nPerState_combos);
            continue
        end

        % Compute variability per PowerPhase combo
        VAR_combos  = nan(1, nCombos);
        MEAN_combos = nan(1, nCombos);
        for c = 1:nCombos

            % Get indices for PowerPhase combo
            thisCombo = comboNames{c};
            idx_thisCombo = idx_combos.(thisCombo);   % trial indices for this combo
            nTrials_thisCombo = numel(idx_thisCombo);
            
            if do_bootstrap
                n_rand = minTrials;
                idx_thisCombo_thisSub = cell(n_boot_iter,1);
                for i = 1:n_boot_iter
                    idx = randperm(nTrials_thisCombo, n_rand);
                    idx_thisCombo_thisSub{i} = idx_thisCombo(idx);
                end
            else
                rand_idx = randperm(nTrials_thisCombo, minTrials);   % balance to minTrials
                idx_thisCombo_thisSub = idx_thisCombo(rand_idx);
            end

            % Compute mean for PowerPhase combo
            MEAN_all = mean(all_data);

            if do_bootstrap
                n_iter = n_boot_iter;
                mean_rand_iter   = zeros(n_iter,1);
                for i = 1:n_iter
                    idx = idx_thisCombo_thisSub{i};
                    mean_rand_iter(i)   = mean(all_data(idx));
                end
                MEAN_combos(c)   = mean(mean_rand_iter);
            else
                MEAN_combos(c) = mean(all_data(idx_thisCombo_thisSub));  % mean amplitude within bin
                % MEAN_combos(c) = mean(all_data(idx_thisCombo));
            end
           
            switch upper(var_method)
                case 'SD' % Standard deviation (assumes log-microvolts)
                    VAR_all  = std(all_data);
                    
                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_thisCombo_thisSub{i};
                            var_rand_iter(i)   = std(all_data(idx));
                        end
                        VAR_combos(c) = mean(var_rand_iter);
                    else
                        VAR_combos(c)  = std(all_data(idx_thisCombo_thisSub));   % spread within bin
                    end

                    
                    % SD of mean-normalize data
                    if bin_mean_normalization
                        VAR_all = std(all_data / MEAN_all);
                        VAR_combos(c)  = std((all_data(idx_thisCombo_thisSub) / MEAN_combos(c)));
                    end
                case 'CV' % Coefficient of variation (assumes raw microvolts)
                    VAR_all  = std(all_data) / mean(all_data);

                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_thisCombo_thisSub{i};
                            var_rand_iter(i)   = std(all_data(idx)) / mean(all_data(idx));
                        end
                        VAR_combos(c) = mean(var_rand_iter);
                    else
                        VAR_combos(c)  = std(all_data(idx_thisCombo_thisSub)) / mean(all_data(idx_thisCombo_thisSub));
                    end

                case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                    VAR_all  = 1.4826 * mad(all_data, 1);
                    VAR_combos(c)  = 1.4826 * mad(all_data(idx_thisCombo_thisSub), 1);
                case 'IQR' % Interquartile range (assumes log-microvolts)
                    VAR_all  = iqr(all_data);
                    VAR_combos(c) = iqr(all_data(idx_thisCombo_thisSub));
                case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)
                    VAR_all  = calculate_CQV(all_data);
                    VAR_combos(b) = calculate_CQV(all_data(idx_thisCombo_thisSub));
                case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                    VAR_all = calculate_TTV_Niu(all_data, nSurrogate);
                    VAR_combos(b) = calculate_TTV_Niu(all_data(idx_thisCombo_thisSub,:), nSurrogate);
                otherwise
                    error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
            end
        end

        % Add Random combos combo
        % idx_allCombos = [];
        % for c = 1:numel(comboNames)
        %     idx_allCombos = [idx_allCombos; idx_combos.(comboNames{c})];
        % end
        % nTrials_allCombos = numel(idx_allCombos);
        % rand_idx = randperm(nTrials_allCombos, minTrials);
        % idx_rand = idx_allCombos(rand_idx);
        
        if do_bootstrap
            n_rand = minTrials;
            idx_rand = cell(n_boot_iter,1);
            for i = 1:n_boot_iter
                rand_idx = randperm(numel(all_data), n_rand);
                idx_rand{i} = rand_idx;
            end
        else
            idx_rand = randperm(numel(all_data), minTrials);
        end

        if do_bootstrap
            n_iter = n_boot_iter;
            mean_rand_iter   = zeros(n_iter,1);
            for i = 1:n_iter
                idx = idx_rand{i};
                mean_rand_iter(i)   = mean(all_data(idx));
            end
            MEAN_rand   = mean(mean_rand_iter);
        else
            MEAN_rand = mean(all_data(idx_rand));
        end

        switch upper(var_method)
            case 'SD' % Standard deviation (assumes log-microvolts)
                
                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = std(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand  = std(all_data(idx_rand));
                end                
                
                % SD of mean-normalize data
                if bin_mean_normalization
                    VAR_rand  = std((all_data(idx_rand) / MEAN_rand));
                end
            case 'CV' % Coefficient of variation (assumes raw microvolts)
               
                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = std(all_data(idx)) / mean(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand = std(all_data(idx_rand)) / mean(all_data(idx_rand));
                end

            case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                VAR_rand = 1.4826 * mad(all_data(idx_rand), 1);
            case 'IQR' % Interquartile range (assumes log-microvolts)

                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = iqr(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand = iqr(all_data(idx_rand));
                end

            case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)

                if do_bootstrap
                    n_iter = n_boot_iter;
                    var_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_rand{i};
                        var_rand_iter(i)   = calculate_CQV(all_data(idx));
                    end
                    VAR_rand   = mean(var_rand_iter);
                else
                    VAR_rand = calculate_CQV(all_data(idx_rand));
                end

            case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                VAR_rand = calculate_TTV_Niu(all_data(idx_rand,:), nSurrogate);
            otherwise
                error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
        end

        % Compute relative variability
        VARrel_combos = (VAR_combos / VAR_all) * 100;
        VARrel_rand = (VAR_rand / VAR_all) * 100;

        % Save results to table
        count = count + 1;
        resultsTable(count).subject      = subID;
        resultsTable(count).site = site;
        resultsTable(count).type    = horzcat('Hjorth_', prestim_feature);
        resultsTable(count).band    = 'alpha';
        %resultsTable(count).ComboNames   = comboNames;
        resultsTable(count).MEAN_all     = MEAN_all;
        resultsTable(count).MEAN_combos  = MEAN_combos;
        resultsTable(count).MEAN_rand  = MEAN_rand;
        resultsTable(count).VAR_method      = var_method;
        resultsTable(count).VAR_all      = VAR_all;
        resultsTable(count).VAR_combos   = VAR_combos;
        resultsTable(count).VAR_rand   = VAR_rand;
        resultsTable(count).VARrel_combos = VARrel_combos;
        resultsTable(count).VARrel_rand = VARrel_rand;

    elseif strcmpi(prestim_feature, 'PowerConnectivity')

        splitsTbl_Power = load_splits_table(DATAIN_SPLIT, sub_num, 1, 'Power', splits);
        splitsTbl_Connectivity = load_splits_table(DATAIN_SPLIT, sub_num, 1, 'Connectivity', splits);

        for fb = 1:numel(bands_of_interest)
            thisBand = bands_of_interest{fb};

            % Load splits table of band
            Tb_Power = splitsTbl_Power(strcmpi(splitsTbl_Power.Freq_band, thisBand), :);
            Tb_Connectivity = splitsTbl_Connectivity(strcmpi(splitsTbl_Connectivity.Freq_band, thisBand), :);

            % Check coherency between both features
            if isempty(Tb_Power) || isempty(Tb_Connectivity)
                warning('No splits found for %s', subID);
                return
            end
            if height(Tb_Power) ~= height(Tb_Connectivity)
                warning('mismatch in trials')
                return
            end
            assert(isequal(Tb_Power.Trial_index, Tb_Connectivity.Trial_index), 'Trial mismatch between Power and Phase');

            % Combine tables of both features
            powerVars = setdiff(Tb_Power.Properties.VariableNames, 'Trial_index');
            connVars = setdiff(Tb_Connectivity.Properties.VariableNames, 'Trial_index');
            overlapVars = intersect(powerVars, connVars);
            for i = 1:length(overlapVars)
                Tb_Power.Properties.VariableNames{overlapVars{i}} = [overlapVars{i} '_Pow'];
                Tb_Connectivity.Properties.VariableNames{overlapVars{i}} = [overlapVars{i} '_Conn'];
            end
            Tb_PowerConnectivity = join(Tb_Power, Tb_Connectivity, 'Keys', 'Trial_index');

            % Get indices for each feature combination
            powers = {'low', 'high'};
            connectivities = {'low', 'high'};
            powStr   = string(Tb_PowerConnectivity.Label_Pow);
            connStr = string(Tb_PowerConnectivity.Label_Conn);

            Tb_PowerConnectivity.Combo = strings(height(Tb_PowerConnectivity),1);
            idx_combos = struct();
            for iP = 1:numel(powers)
                for iPh = 1:numel(connectivities)
                    comboName = sprintf('%s_%s', powers{iP}, connectivities{iPh});
                    idx = strcmpi(powStr,   powers{iP}) & ...
                        strcmpi(connStr, connectivities{iPh});
                    idx_combos.(comboName) = find(idx);
                    Tb_PowerConnectivity.Combo(idx) = comboName;
                end
            end
            Tb_PowerConnectivity.Combo = categorical(Tb_PowerConnectivity.Combo);
            counts = countcats(Tb_PowerConnectivity.Combo);

            % Pseudo N-match all conditions
            comboNames = fieldnames(idx_combos);
            nCombos    = numel(comboNames);
            trialsPerCombo = cellfun(@(name) numel(idx_combos.(name)), comboNames);
            [minTrials, minTrials_idx] = min(trialsPerCombo);
            if minTrials < nPerState_combos
                fprintf('%s does not have enough trials for combo %s (n = %d < %d)\n', ...
                    subID, comboNames{minTrials_idx}, minTrials, nPerState_combos);
                continue
            end

            % Compute variability per PowerPhase combo
            VAR_combos  = nan(1, nCombos);
            MEAN_combos = nan(1, nCombos);
            for c = 1:nCombos

                % Get indices for PowerPhase combo
                thisCombo = comboNames{c};
                idx_thisCombo = idx_combos.(thisCombo);   % trial indices for this combo
                nTrials_thisCombo = numel(idx_thisCombo);

                if do_bootstrap
                    n_rand = minTrials;
                    idx_thisCombo_thisSub = cell(n_boot_iter,1);
                    for i = 1:n_boot_iter
                        idx = randperm(nTrials_thisCombo, n_rand);
                        idx_thisCombo_thisSub{i} = idx_thisCombo(idx);
                    end
                else
                    rand_idx = randperm(nTrials_thisCombo, minTrials);   % balance to minTrials
                    idx_thisCombo_thisSub = idx_thisCombo(rand_idx);
                end

                % Compute mean for PowerPhase combo
                MEAN_all = mean(all_data);

                if do_bootstrap
                    n_iter = n_boot_iter;
                    mean_rand_iter   = zeros(n_iter,1);
                    for i = 1:n_iter
                        idx = idx_thisCombo_thisSub{i};
                        mean_rand_iter(i)   = mean(all_data(idx));
                    end
                    MEAN_combos(c)   = mean(mean_rand_iter);
                else
                    MEAN_combos(c) = mean(all_data(idx_thisCombo_thisSub));  % mean amplitude within bin
                    % MEAN_combos(c) = mean(all_data(idx_thisCombo));
                end

                switch upper(var_method)
                    case 'SD' % Standard deviation (assumes log-microvolts)
                        VAR_all  = std(all_data);

                        if do_bootstrap
                            n_iter = n_boot_iter;
                            var_rand_iter   = zeros(n_iter,1);
                            for i = 1:n_iter
                                idx = idx_thisCombo_thisSub{i};
                                var_rand_iter(i)   = std(all_data(idx));
                            end
                            VAR_combos(c) = mean(var_rand_iter);
                        else
                            VAR_combos(c)  = std(all_data(idx_thisCombo_thisSub));   % spread within bin
                        end


                        % SD of mean-normalize data
                        if bin_mean_normalization
                            VAR_all = std(all_data / MEAN_all);
                            VAR_combos(c)  = std((all_data(idx_thisCombo_thisSub) / MEAN_combos(c)));
                        end
                    case 'CV' % Coefficient of variation (assumes raw microvolts)
                        VAR_all  = std(all_data) / mean(all_data);

                        if do_bootstrap
                            n_iter = n_boot_iter;
                            var_rand_iter   = zeros(n_iter,1);
                            for i = 1:n_iter
                                idx = idx_thisCombo_thisSub{i};
                                var_rand_iter(i)   = std(all_data(idx)) / mean(all_data(idx));
                            end
                            VAR_combos(c) = mean(var_rand_iter);
                        else
                            VAR_combos(c)  = std(all_data(idx_thisCombo_thisSub)) / mean(all_data(idx_thisCombo_thisSub));
                        end

                    case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                        VAR_all  = 1.4826 * mad(all_data, 1);
                        VAR_combos(c)  = 1.4826 * mad(all_data(idx_thisCombo_thisSub), 1);
                    case 'IQR' % Interquartile range (assumes log-microvolts)
                        VAR_all  = iqr(all_data);
                        VAR_combos(c) = iqr(all_data(idx_thisCombo_thisSub));
                    case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)
                        VAR_all  = calculate_CQV(all_data);
                        VAR_combos(b) = calculate_CQV(all_data(idx_thisCombo_thisSub));
                    case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                        VAR_all = calculate_TTV_Niu(all_data, nSurrogate);
                        VAR_combos(b) = calculate_TTV_Niu(all_data(idx_thisCombo_thisSub,:), nSurrogate);
                    otherwise
                        error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
                end
            end

            % Add Random combos combo
            % idx_allCombos = [];
            % for c = 1:numel(comboNames)
            %     idx_allCombos = [idx_allCombos; idx_combos.(comboNames{c})];
            % end
            % nTrials_allCombos = numel(idx_allCombos);
            % rand_idx = randperm(nTrials_allCombos, minTrials);
            % idx_rand = idx_allCombos(rand_idx);

            if do_bootstrap
                n_rand = minTrials;
                idx_rand = cell(n_boot_iter,1);
                for i = 1:n_boot_iter
                    rand_idx = randperm(numel(all_data), n_rand);
                    idx_rand{i} = rand_idx;
                end
            else
                idx_rand = randperm(numel(all_data), minTrials);
            end

            if do_bootstrap
                n_iter = n_boot_iter;
                mean_rand_iter   = zeros(n_iter,1);
                for i = 1:n_iter
                    idx = idx_rand{i};
                    mean_rand_iter(i)   = mean(all_data(idx));
                end
                MEAN_rand   = mean(mean_rand_iter);
            else
                MEAN_rand = mean(all_data(idx_rand));
            end

            switch upper(var_method)
                case 'SD' % Standard deviation (assumes log-microvolts)

                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_rand{i};
                            var_rand_iter(i)   = std(all_data(idx));
                        end
                        VAR_rand   = mean(var_rand_iter);
                    else
                        VAR_rand  = std(all_data(idx_rand));
                    end

                    % SD of mean-normalize data
                    if bin_mean_normalization
                        VAR_rand  = std((all_data(idx_rand) / MEAN_rand));
                    end
                case 'CV' % Coefficient of variation (assumes raw microvolts)

                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_rand{i};
                            var_rand_iter(i)   = std(all_data(idx)) / mean(all_data(idx));
                        end
                        VAR_rand   = mean(var_rand_iter);
                    else
                        VAR_rand = std(all_data(idx_rand)) / mean(all_data(idx_rand));
                    end

                case 'MAD' % Mean absolute deviation (assumes log-microvolts)
                    VAR_rand = 1.4826 * mad(all_data(idx_rand), 1);
                case 'IQR' % Interquartile range (assumes log-microvolts)

                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_rand{i};
                            var_rand_iter(i)   = iqr(all_data(idx));
                        end
                        VAR_rand   = mean(var_rand_iter);
                    else
                        VAR_rand = iqr(all_data(idx_rand));
                    end

                case 'CQV' % Coefficient of quartile variation: (Q3-Q1)/(Q3+Q1)

                    if do_bootstrap
                        n_iter = n_boot_iter;
                        var_rand_iter   = zeros(n_iter,1);
                        for i = 1:n_iter
                            idx = idx_rand{i};
                            var_rand_iter(i)   = calculate_CQV(all_data(idx));
                        end
                        VAR_rand   = mean(var_rand_iter);
                    else
                        VAR_rand = calculate_CQV(all_data(idx_rand));
                    end

                case 'TTV_NIU' % Cross-correlation matrices (Niu et al., 2025)
                    VAR_rand = calculate_TTV_Niu(all_data(idx_rand,:), nSurrogate);
                otherwise
                    error('Unknown var_method: %s. Use SD, CV, MAD, or IQR.', var_method)
            end

            % Compute relative variability
            VARrel_combos = (VAR_combos / VAR_all) * 100;
            VARrel_rand = (VAR_rand / VAR_all) * 100;

            % Save results to table
            count = count + 1;
            resultsTable(count).subject      = subID;
            resultsTable(count).site = site;
            resultsTable(count).type    = horzcat('Hjorth_', prestim_feature);
            resultsTable(count).band    = thisBand;
            %resultsTable(count).ComboNames   = comboNames;
            resultsTable(count).MEAN_all     = MEAN_all;
            resultsTable(count).MEAN_combos  = MEAN_combos;
            resultsTable(count).MEAN_rand  = MEAN_rand;
            resultsTable(count).VAR_method      = var_method;
            resultsTable(count).VAR_all      = VAR_all;
            resultsTable(count).VAR_combos   = VAR_combos;
            resultsTable(count).VAR_rand   = VAR_rand;
            resultsTable(count).VARrel_combos = VARrel_combos;
            resultsTable(count).VARrel_rand = VARrel_rand;
        end
    end
end
T_all = struct2table(resultsTable);

% Save output table
if save_output & ~strcmpi(plots_show, 'MEAN_AMP')
    if strcmpi(response, 'TEP')
        outfile = sprintf('T_VAR_allblocks_%s_%s_%s_%s_%s_%s_%s_%s_%s.csv', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''));
        outpath = DATAIN_EP;
    elseif strcmpi(response, 'MEP')
        outfile = sprintf('T_VAR_allblocks_%s_%s_%s_%s_%s_%s_%s.csv', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(var_method), site, splits, strjoin(bands_of_interest, ''));
        outpath = DATAIN_EP;
    end
    writetable(T_all, fullfile(outpath, outfile));
end

%% ----------- GROUP-LEVEL SUMMARY -----------
format bank
if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity') || strcmpi(prestim_feature, 'PowerConnectivity')

    % Variables to average across subjects Variables to average across
    % subjects
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')
        vars2avg = {'MEAN_all', 'MEAN_high', 'MEAN_low', 'MEAN_rand', ...
            'VAR_all', 'VAR_high', 'VAR_low', 'VAR_rand', ...
            'VARrel_high', 'VARrel_low', 'VARrel_rand', ...
            'ratio_high', 'ratio_low', 'logratio_high', 'logratio_low'};
    else
        vars2avg = {'MEAN_all', 'MEAN_combos', 'MEAN_rand', ...
            'VAR_all', 'VAR_combos', 'VAR_rand', ...
            'VARrel_combos', 'VARrel_rand'};
    end

    % Group by frequency band
    G = grpstats(T_all, "band", {@mean,@std,@numel}, 'DataVars', vars2avg);

    % Compute SEM for each
    for v = vars2avg
        var = v{1};
        meanName = ['mean_' var];
        stdName  = ['std_' var];
        nName    = ['numel_' var];
        semName  = ['sem_' var];

        if ismember(stdName, G.Properties.VariableNames) && ismember(nName, G.Properties.VariableNames)
            G.(semName) = G.(stdName) ./ sqrt(G.(nName));
        else
            warning('Missing columns for %s', var);
        end
    end

elseif strcmpi(prestim_feature, 'Phase')

    MEAN_phasebins_mat = cat(1, T_all.MEAN_phasebins);
    VAR_phasebins_mat  = cat(1, T_all.VAR_phasebins);
    VARrel_phasebins_mat  = cat(1, T_all.VARrel_phasebins);

    % Compute group averages
    mean_MEAN_all = mean(T_all.MEAN_all);
    sem_MEAN_all = std(T_all.MEAN_all,  'omitnan') ./ sqrt(size(T_all.MEAN_all,1));
    mean_MEAN_phasebins = mean(MEAN_phasebins_mat, 'omitnan');
    sem_MEAN_phasebins  = std(MEAN_phasebins_mat,  'omitnan') ./ sqrt(size(MEAN_phasebins_mat,1));
    mean_MEAN_rand = mean(T_all.MEAN_rand, 'omitnan');
    sem_MEAN_rand  = std(T_all.MEAN_rand,  'omitnan') ./ sqrt(size(T_all.MEAN_rand,1));
    mean_VAR_all = mean(T_all.VAR_all);
    sem_VAR_all = std(T_all.VAR_all,  'omitnan') ./ sqrt(size(T_all.VAR_all,1));
    mean_VAR_phasebins = mean(VAR_phasebins_mat, 'omitnan');
    sem_VAR_phasebins  = std(VAR_phasebins_mat,  'omitnan') ./ sqrt(size(VAR_phasebins_mat,1));
    mean_VAR_rand = mean(T_all.VAR_rand, 'omitnan');
    sem_VAR_rand  = std(T_all.VAR_rand,  'omitnan') ./ sqrt(size(T_all.VAR_rand,1));
    mean_VARrel_phasebins = mean(VARrel_phasebins_mat, 'omitnan');
    sem_VARrel_phasebins  = std(VARrel_phasebins_mat,  'omitnan') ./ sqrt(size(VARrel_phasebins_mat,1));
    mean_VARrel_rand = mean(T_all.VARrel_rand, 'omitnan');
    sem_VARrel_rand  = std(T_all.VARrel_rand,  'omitnan') ./ sqrt(size(T_all.VARrel_rand,1));


    G = table(mean_MEAN_all, sem_MEAN_all, ...
        mean_MEAN_phasebins, sem_MEAN_phasebins, ...
        mean_MEAN_rand, sem_MEAN_rand, ...
        mean_VAR_all, sem_VAR_all, ...
        mean_VAR_phasebins, sem_VAR_phasebins, ...
        mean_VAR_rand, sem_VAR_rand, ...
        mean_VARrel_phasebins, sem_VARrel_phasebins, ...
        mean_VARrel_rand, sem_VARrel_rand, ...
        'VariableNames', {'mean_MEAN_all', 'sem_MEAN_all', ...
        'mean_MEAN_phasebins', 'sem_MEAN_phasebins', ...
        'mean_MEAN_rand', 'sem_MEAN_rand', ...
        'mean_VAR_all', 'sem_VAR_all', ...
        'mean_VAR_phasebins', 'sem_VAR_phasebins', ...
        'mean_VAR_rand', 'sem_VAR_rand', ...
        'mean_VARrel_phasebins', 'sem_VARrel_phasebins', ...
        'mean_VARrel_rand', 'sem_VARrel_rand'});

elseif strcmpi(prestim_feature, 'PowerPhase')

    MEAN_combos_mat = cat(1, T_all.MEAN_combos);
    VAR_combos_mat  = cat(1, T_all.VAR_combos);
    VARrel_combos_mat  = cat(1, T_all.VARrel_combos);

    % Compute group averages
    mean_MEAN_all = mean(T_all.MEAN_all);
    sem_MEAN_all = std(T_all.MEAN_all,  'omitnan') ./ sqrt(size(T_all.MEAN_all,1));
    mean_MEAN_combos = mean(MEAN_combos_mat, 'omitnan');
    sem_MEAN_combos  = std(MEAN_combos_mat,  'omitnan') ./ sqrt(size(MEAN_combos_mat,1));
    mean_MEAN_rand = mean(T_all.MEAN_rand, 'omitnan');
    sem_MEAN_rand  = std(T_all.MEAN_rand,  'omitnan') ./ sqrt(size(T_all.MEAN_rand,1));
    mean_VAR_all = mean(T_all.VAR_all);
    sem_VAR_all = std(T_all.VAR_all,  'omitnan') ./ sqrt(size(T_all.VAR_all,1));
    mean_VAR_combos = mean(VAR_combos_mat, 'omitnan');
    sem_VAR_combos  = std(VAR_combos_mat,  'omitnan') ./ sqrt(size(VAR_combos_mat,1));
    mean_VAR_rand = mean(T_all.VAR_rand, 'omitnan');
    sem_VAR_rand  = std(T_all.VAR_rand,  'omitnan') ./ sqrt(size(T_all.VAR_rand,1));
    mean_VARrel_combos = mean(VARrel_combos_mat, 'omitnan');
    sem_VARrel_combos  = std(VARrel_combos_mat,  'omitnan') ./ sqrt(size(VARrel_combos_mat,1));
    mean_VARrel_rand = mean(T_all.VARrel_rand, 'omitnan');
    sem_VARrel_rand  = std(T_all.VARrel_rand,  'omitnan') ./ sqrt(size(T_all.VARrel_rand,1));


    G = table(mean_MEAN_all, sem_MEAN_all, ...
        mean_MEAN_combos, sem_MEAN_combos, ...
        mean_MEAN_rand, sem_MEAN_rand, ...
        mean_VAR_all, sem_VAR_all, ...
        mean_VAR_combos, sem_VAR_combos, ...
        mean_VAR_rand, sem_VAR_rand, ...
        mean_VARrel_combos, sem_VARrel_combos, ...
        mean_VARrel_rand, sem_VARrel_rand, ...
        'VariableNames', {'mean_MEAN_all', 'sem_MEAN_all', ...
        'mean_MEAN_combos', 'sem_MEAN_combos', ...
        'mean_MEAN_rand', 'sem_MEAN_rand', ...
        'mean_VAR_all', 'sem_VAR_all', ...
        'mean_VAR_combos', 'sem_VAR_combos', ...
        'mean_VAR_rand', 'sem_VAR_rand', ...
        'mean_VARrel_combos', 'sem_VARrel_combos', ...
        'mean_VARrel_rand', 'sem_VARrel_rand'});
end

% % Save group-level summary
% if save_output
%     if strcmpi(response, 'TEP')
%         outfile = sprintf('T_VAR_allblocks_grouplevel_%s_%s_%s_%s_%s_%s_%s_%s_%s.csv', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''));
%         outpath = DATAIN_EP;
%     elseif strcmpi(response, 'MEP')
%         outfile = sprintf('T_VAR_allblocks_grouplevel_%s_%s_%s_%s_%s_%s_%s.csv', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(var_method), site, splits, strjoin(bands_of_interest, ''));
%         outpath = DATAIN_EP;
%     end
%     writetable(G, fullfile(outpath, outfile));
% end

%% ----------- PLOTTING -----------

if generate_plots
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')

        % Define frequency bands in desired order + pretty names
        allBandInfo = { ...
            'theta', 'Theta';
            'alpha', 'Mu';
            'beta', 'Beta';
            'gamma', 'Gamma'};

        [~, idx] = ismember(bands_of_interest, allBandInfo(:,1));
        idx = idx(idx > 0); % Remove zeros (bands not found)
        bandInfo = allBandInfo(idx, :);
        bands     = string(bandInfo(:,1));  % raw band names in desired order
        plotNames = string(bandInfo(:,2));  % pretty labels
        [~, order] = ismember(bands, string(G.band));
        G = G(order,:);
        nBands = height(G);

        % Define plot colours
        col_high = [253, 70, 101] / 255;    % Vibrant red
        col_low = [59, 82, 139] / 255;    % Medium blue
        col_rand = [120, 120, 130] / 255;  % Cool gray with slight blue tint

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (strcmpi(plots_show,'VARrel') || strcmpi(plots_show, 'VAR')) & length(bands_of_interest) > 1
            figure('Name','Relative dispersion across bands (all-trials)','Color','w'); hold on;

            % Number of bands and conditions
            nBands = length(bands_of_interest);
            conditions = {'low', 'high'};
            nConditions = length(conditions);

            % Define colors for conditions
            condition_colors = [col_low; col_high];

            % Prepare data structure using cell arrays
            band_data = cell(nBands + 1, 1); % +1 for random

            % Add connecting boxplot line
            summary_low  = nan(nBands,1);
            summary_high = nan(nBands,1);
            x_low  = nan(nBands,1);
            x_high = nan(nBands,1);

            % Extract data for each band
            scatter_offset = 0;
            band_step = 1;      % Fixed spacing between band centers
            cond_sep  = 0.25;   % Separation between low/high within a band
            jitter_amount = 0.02;

            for b = 1:nBands
                band_name = bands_of_interest{b};
                band_table = T_all(strcmp(T_all.band, band_name), :);

                data_low  = band_table.([plots_show, '_low']);
                data_high = band_table.([plots_show, '_high']);

                % positions (clean) - use simple integer spacing
                band_center = b;
                positions   = band_center + [-cond_sep/2, +cond_sep/2];

                % Mean and SEM
                mu_high  = mean(data_high, 'omitnan');
                sem_high = std(data_high, 'omitnan') / sqrt(sum(~isnan(data_high)));
                mu_low   = mean(data_low,  'omitnan');
                sem_low  = std(data_low,  'omitnan') / sqrt(sum(~isnan(data_low)));

                errorbar(positions(1), mu_low, sem_low, 'o', ...
                    'Color', col_low, 'MarkerFaceColor', col_low, 'MarkerEdgeColor', col_low, ...
                    'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

                errorbar(positions(2), mu_high, sem_high, 'o', ...
                    'Color', col_high, 'MarkerFaceColor', col_high, 'MarkerEdgeColor', col_high, ...
                    'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

                % Store data in cell array
                band_data{b}.name = band_name;
                band_data{b}.center = band_center;
                band_data{b}.positions = positions;
                band_data{b}.data_low = data_low;
                band_data{b}.data_high = data_high;

                summary_low(b)  = mean(data_low,  'omitnan');
                summary_high(b) = mean(data_high, 'omitnan');
                x_low(b)  = positions(1);
                x_high(b) = positions(2);
            end

            % % === Add random (single-box) ===
            % data_random = unique(T_all.([plots_show, '_rand']));    % Random (same for all bands)
            % mu_rand  = mean(data_random, 'omitnan');
            % sem_rand = std(data_random, 'omitnan') / sqrt(sum(~isnan(data_random)));
            % 
            % % Place random at the next integer position
            % pos_rand = nBands + 1;
            % 
            % % Plot error bar
            % errorbar(pos_rand, mu_rand, sem_rand, ...
            %     'o', ...
            %     'Color', col_rand, ...
            %     'MarkerFaceColor', col_rand, ...
            %     'MarkerEdgeColor', col_rand, ...
            %     'MarkerSize', 10, ...
            %     'LineWidth', 4, ...
            %     'CapSize', 25);
            % hold on;
            % 
            % % Scatter points with jitter for random
            % scatter(pos_rand + scatter_offset + jitter_amount*(rand(size(data_random))-0.5), ...
            %     data_random, 25, ...
            %     'MarkerFaceColor', col_rand, 'MarkerEdgeColor', 'k', ...
            %     'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
            % 
            % % Store for tick placement and separators
            % band_data{nBands+1}.name = 'Random';
            % band_data{nBands+1}.positions = pos_rand;

            % Add individual data points with jitter
            for b = 1:nBands
                positions = band_data{b}.positions + scatter_offset;
                % Scatter points for each condition
                scatter(positions(1) + jitter_amount*(rand(size(band_data{b}.data_low))-0.5), ...
                    band_data{b}.data_low, 25, ...
                    'MarkerFaceColor', col_low, 'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
                scatter(positions(2) + jitter_amount*(rand(size(band_data{b}.data_high))-0.5), ...
                    band_data{b}.data_high, 25, ...
                    'MarkerFaceColor', col_high, 'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
            end

            % === Add connection line between means ===
            % plot(x_low, summary_low, ...
            %     '-', 'Color', [col_low], ...
            %     'MarkerFaceColor',  col_low, 'LineWidth', 2);
            % plot(x_high, summary_high, ...
            %     '-', 'Color', [col_high], ...
            %     'MarkerFaceColor',  col_high, 'LineWidth', 2);

            % === Formatting ===
            tick_positions = 1:(nBands+1); % Simple integer positions
            band_labels = [plotNames(:); {'Random'}];
            set(gca, 'XTick', tick_positions, 'XTickLabel', band_labels, ...
                'FontSize', 12, 'LineWidth', 1);
            ylim([50 150])

            % Add separator lines between bands
            for b = 1:nBands-1
                line_pos = b + 0.5;   % midpoint between band b and b+1
                plot([line_pos line_pos], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
            end

            % Create custom legend
            h_legend = zeros(2,1);
            h_legend(1) = plot(nan, nan, 's', 'Color', col_high, 'MarkerFaceColor', col_high, 'MarkerSize', 12);
            h_legend(2) = plot(nan, nan, 's', 'Color', col_low, 'MarkerFaceColor', col_low, 'MarkerSize', 12);
            legend(h_legend, {'High', 'Low'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 14);

            set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
            ax = gca;
            set(ax, 'FontSize', 14);
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 2;
            fontsize(gcf, scale=1.5)
            fontname('Arial')
            xlim([0.5 4.5])

            yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

            % xlim_range = xlim;
            % % plot(xlim_range, [mu_rand mu_rand], ':', ...
            % %     'Color', [0.7 0.7 0.7], ...
            % %     'LineWidth', 1 ,'HandleVisibility', 'off');
            %  plot(xlim_range, [100 100], ':', ...
            %     'Color', [0.7 0.7 0.7], ...
            %     'LineWidth', 1 ,'HandleVisibility', 'off');

            %xlabel(sprintf('%s', prestim_feature), 'FontSize', 16, 'FontWeight','bold');
            %ylabel(sprintf('Normalized within-state %s', var_method), 'FontSize', 16, 'FontWeight','bold');
            %xlab = ax.XLabel;
            %xlab.Position(2) = xlab.Position(2) - 1.1;
            %ylab = ax.YLabel;
            %ylab.Position(1) = ylab.Position(1) - 0.1;

            % === Add significance
            [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, strjoin(bands_of_interest, ''), reftep);
            sigBands_fields = fields(sigBands);

            if length(bands_of_interest) ~= 4
                allBands = {'theta','alpha','beta','gamma'};
                keepIdx = find(ismember(allBands, bands_of_interest));

                fn = fieldnames(sigBands);
                for i = 1:numel(fn)
                    key = fn{i};
                    v = sigBands.(key);
                    s = starText.(key);

                    if isempty(v)
                        sigBands.(key) = [];
                        starText.(key) = {};
                        continue
                    end

                    mask = ismember(v, keepIdx);
                    sigBands.(key) = v(mask);
                    starText.(key) = s(mask);
                end

                map = zeros(1, numel(allBands));
                map(keepIdx) = 1:numel(keepIdx);
            end

            for j = 1:size(sigBands_fields,1) % lowhigh lowrand highrand
                thisPair = sigBands_fields{j};
                sign_thisPair = sigBands.(thisPair);
                starText_thisPair = starText.(thisPair);

                if length(sign_thisPair) ~= length(starText_thisPair)
                    disp('Mismatch significance')
                    return
                end
                if isempty(sign_thisPair)
                    continue
                end

                for i = 1:length(sign_thisPair)
                    sign_now = sign_thisPair(i);
                    starText_now = starText_thisPair{i};

                    if numel(bands_of_interest) ~= 4
                        sign_now = map(sign_now);
                        if sign_now == 0
                            continue
                        end
                    end

                    %x_rand = pos_rand;
                    if strcmpi(thisPair, 'lowhigh')
                        x_first = x_low(sign_now);
                        x_second =  x_high(sign_now);
                    elseif strcmpi(thisPair, 'lowrand')
                        x_first = x_low(sign_now);
                        x_second = x_rand;
                        yOffset = yOffset - 2;
                    else
                        x_first = x_high(sign_now);
                        x_second = x_rand;
                        yOffset = yOffset - 2;
                    end

                    % Draw significance
                    ylim_range = ylim;
                    tickLength = 0.005 * diff(ylim_range);
                    y_top = ylim_range(2) - yOffset;
                    y_bottom = y_top - tickLength;

                    line([x_first, x_first, x_second, x_second], ...
                        [y_bottom, y_top, y_top, y_bottom], ...
                        'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
                    text(((x_first + x_second) / 2), ylim_range(2) - yOffset + 0.5, starText_now, ...
                        'HorizontalAlignment', 'center', 'FontSize', 28, 'FontWeight', 'bold', ...
                        'HandleVisibility', 'off');
                end
            end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''))), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
                end
            end
        end

        if (strcmpi(plots_show,'VARrel') || strcmpi(plots_show, 'VAR')) && length(bands_of_interest) == 1
            figure('Name','Relative dispersion across bands (all-trials)','Color','w'); hold on;
            jitter_amount = 0.02;
           
            % Single band case
            band_name = bands_of_interest{1};  % Get the single band
            band_table = T_all(strcmp(T_all.band, band_name), :);

            % Extract data
            data_low  = band_table.([plots_show, '_low']);
            data_high = band_table.([plots_show, '_high']);
            data_random = unique(T_all.([plots_show, '_rand']));  % Random (same for all bands)

            % Positions - three separate x positions for Low, High, Random
            positions = [1, 2, 3];  % Each condition gets its own x-tick

            % Calculate mean and SEM
            mu_low   = mean(data_low,  'omitnan');
            sem_low  = std(data_low,  'omitnan') / sqrt(sum(~isnan(data_low)));

            mu_high  = mean(data_high, 'omitnan');
            sem_high = std(data_high, 'omitnan') / sqrt(sum(~isnan(data_high)));

            mu_rand  = mean(data_random, 'omitnan');
            sem_rand = std(data_random, 'omitnan') / sqrt(sum(~isnan(data_random)));

            % Plot error bars
            errorbar(positions(1), mu_low, sem_low, 'o', ...
                'Color', col_low, 'MarkerFaceColor', col_low, 'MarkerEdgeColor', col_low, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            hold on;
            errorbar(positions(2), mu_high, sem_high, 'o', ...
                'Color', col_high, 'MarkerFaceColor', col_high, 'MarkerEdgeColor', col_high, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);
            hold on;

            errorbar(positions(3), mu_rand, sem_rand, 'o', ...
                'Color', col_rand, 'MarkerFaceColor', col_rand, 'MarkerEdgeColor', col_rand, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            hold on;

            % Scatter points with jitter
            scatter(positions(1) + jitter_amount*(rand(size(data_low))-0.5), ...
                data_low, 25, ...
                'MarkerFaceColor', col_low, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
            scatter(positions(2) + jitter_amount*(rand(size(data_high))-0.5), ...
                data_high, 25, ...
                'MarkerFaceColor', col_high, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
            scatter(positions(3) + jitter_amount*(rand(size(data_random))-0.5), ...
                data_random, 25, ...
                'MarkerFaceColor', col_rand, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);

            % === Formatting ===
            set(gca, 'XTick', positions, 'XTickLabel', {'Low', 'High', 'Random'}, ...
                'FontSize', 12, 'LineWidth', 1);
            ylim([70 130])


            yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
    
            box on;
            set(gcf, 'Units', 'normalized', 'Position', [0.666796875,0.152083333333333,0.238203125,0.772916666666667]);
            ax = gca;
            ax.Box = 'off';
            ax.TickDir = 'out';
            set(ax, 'FontSize', 14);
            ax.LineWidth = 2;
            fontsize(gcf,scale=1.5)
            fontname('Arial')
            
            xlim([0, 4])
            % xlim_range = xlim;
            % % plot(xlim_range, [mu_rand mu_rand], ':', ...
            % %     'Color', [0.7 0.7 0.7], ...
            % %     'LineWidth', 0.00000001 ,'HandleVisibility', 'off');
            % plot(xlim_range, [100 100], ':', ...
            %     'Color', [0.7 0.7 0.7], ...
            %     'LineWidth', 0.00000001 ,'HandleVisibility', 'off');

            % ax.Position = [0.18 0.13 0.72 0.77];
            % xlim_range = xlim;
            % % plot(xlim_range, [mu_rand mu_rand], ':', ...
            % %     'Color', [0.7 0.7 0.7], ...
            % %     'LineWidth', 1 ,'HandleVisibility', 'off');
            % plot(xlim_range, [100 100], ':', ...
            %     'Color', [0.7 0.7 0.7], ...
            %     'LineWidth', 0.00000001 ,'HandleVisibility', 'off');


            % xlabel(sprintf('%s bin', prestim_feature), 'FontSize', 16, 'FontWeight','bold');
            % ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 16, 'FontWeight','bold');
            % xlab = ax.XLabel;
            % xlab.Position(2) = xlab.Position(2) - 1.1;
            % ylab = ax.YLabel;
            % ylab.Position(1) = ylab.Position(1) - 0.1;

            % === Add significance
            [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, strjoin(bands_of_interest, ''), reftep);
            sigBands_fields = fields(sigBands);

            if length(bands_of_interest) ~= 4
                allBands = {'theta','alpha','beta','gamma'};
                keepIdx = find(ismember(allBands, bands_of_interest));

                fn = fieldnames(sigBands);
                for i = 1:numel(fn)
                    key = fn{i};
                    v = sigBands.(key);
                    s = starText.(key);

                    if isempty(v)
                        sigBands.(key) = [];
                        starText.(key) = {};
                        continue
                    end

                    mask = ismember(v, keepIdx);
                    sigBands.(key) = v(mask);
                    starText.(key) = s(mask);
                end

                map = zeros(1, numel(allBands));
                map(keepIdx) = 1:numel(keepIdx);
            end

            for j = 1:size(sigBands_fields,1) % lowhigh lowrand highrand
                thisPair = sigBands_fields{j};
                sign_thisPair = sigBands.(thisPair);
                starText_thisPair = starText.(thisPair);

                if length(sign_thisPair) ~= length(starText_thisPair)
                    disp('Mismatch significance')
                    return
                end
                if isempty(sign_thisPair)
                    continue
                end

                for i = 1:length(sign_thisPair)
                    sign_now = sign_thisPair(i);
                    starText_now = starText_thisPair{i};

                    if strcmpi(thisPair, 'lowhigh')
                        x_first = positions(1);
                        x_second =  positions(2);
                    elseif strcmpi(thisPair, 'lowrand')
                        x_first = positions(1);
                        x_second = positions(3);
                        yOffset = yOffset - 2;
                    else
                        x_first = positions(2);
                        x_second = positions(3);
                        yOffset = yOffset - 2;
                    end

                    ylim_range = ylim;
                    tickLength = 0.005 * diff(ylim_range);
                    y_top = ylim_range(2) - yOffset;
                    y_bottom = y_top - tickLength;

                    line([x_first, x_first, x_second, x_second], ...
                        [y_bottom, y_top, y_top, y_bottom], ...
                        'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
                    text(((x_first + x_second) / 2), ylim_range(2) - yOffset + 0.5, starText_now, ...
                        'HorizontalAlignment', 'center', 'FontSize', 28, 'FontWeight', 'bold', ...
                        'HandleVisibility', 'off');
                end
            end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''))), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
                end
            end
        end        

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FINAL PLOTS MEANS
        if strcmpi(plots_show, 'MEAN_AMP') & length(bands_of_interest) > 1
            figure('Name','Means across bands (all-trials)','Color','w'); hold on;

            % Number of bands and conditions
            nBands = length(bands_of_interest);
            conditions = {'low', 'high'};
            nConditions = length(conditions);

            % Define colors for conditions
            condition_colors = [col_low; col_high];

            % Prepare data structure using cell arrays
            band_data = cell(nBands + 1, 1); % +1 for random

            % Add connecting boxplot line
            summary_low  = nan(nBands,1);
            summary_high = nan(nBands,1);
            x_low  = nan(nBands,1);
            x_high = nan(nBands,1);

            % Extract data for each band
            scatter_offset = 0;
            band_step = 1;      % Fixed spacing between band centers
            cond_sep  = 0.25;   % Separation between low/high within a band
            jitter_amount = 0.02;

            for b = 1:nBands
                band_name = bands_of_interest{b};
                band_table = T_all(strcmp(T_all.band, band_name), :);

                data_low  = band_table.(['MEAN_low']);
                data_high = band_table.(['MEAN_high']);

                % positions (clean) - use simple integer spacing
                band_center = b;
                positions   = band_center + [-cond_sep/2, +cond_sep/2];

                % Mean and SEM
                mu_high  = mean(data_high, 'omitnan');
                sem_high = std(data_high, 'omitnan') / sqrt(sum(~isnan(data_high)));
                mu_low   = mean(data_low,  'omitnan');
                sem_low  = std(data_low,  'omitnan') / sqrt(sum(~isnan(data_low)));

                errorbar(positions(1), mu_low, sem_low, 'o', ...
                    'Color', col_low, 'MarkerFaceColor', col_low, 'MarkerEdgeColor', col_low, ...
                    'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

                errorbar(positions(2), mu_high, sem_high, 'o', ...
                    'Color', col_high, 'MarkerFaceColor', col_high, 'MarkerEdgeColor', col_high, ...
                    'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

                % Store data in cell array
                band_data{b}.name = band_name;
                band_data{b}.center = band_center;
                band_data{b}.positions = positions;
                band_data{b}.data_low = data_low;
                band_data{b}.data_high = data_high;

                summary_low(b)  = mean(data_low,  'omitnan');
                summary_high(b) = mean(data_high, 'omitnan');
                x_low(b)  = positions(1);
                x_high(b) = positions(2);
            end

            % === Add random (single-box) ===
            data_random = unique(T_all.(['MEAN_rand']));    % Random (same for all bands)
            mu_rand  = mean(data_random, 'omitnan');
            sem_rand = std(data_random, 'omitnan') / sqrt(sum(~isnan(data_random)));

            % Place random at the next integer position
            pos_rand = nBands + 1;

            % Plot error bar
            errorbar(pos_rand, mu_rand, sem_rand, ...
                'o', ...
                'Color', col_rand, ...
                'MarkerFaceColor', col_rand, ...
                'MarkerEdgeColor', col_rand, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);
            hold on;

            % Scatter points with jitter for random
            scatter(pos_rand + scatter_offset + jitter_amount*(rand(size(data_random))-0.5), ...
                data_random, 25, ...
                'MarkerFaceColor', col_rand, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);

            % Store for tick placement and separators
            band_data{nBands+1}.name = 'Random';
            band_data{nBands+1}.positions = pos_rand;

            % Add individual data points with jitter
            for b = 1:nBands
                positions = band_data{b}.positions + scatter_offset;
                % Scatter points for each condition
                scatter(positions(1) + jitter_amount*(rand(size(band_data{b}.data_low))-0.5), ...
                    band_data{b}.data_low, 25, ...
                    'MarkerFaceColor', col_low, 'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
                scatter(positions(2) + jitter_amount*(rand(size(band_data{b}.data_high))-0.5), ...
                    band_data{b}.data_high, 25, ...
                    'MarkerFaceColor', col_high, 'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
            end

            % === Add connection line between means ===
            plot(x_low, summary_low, ...
                '-', 'Color', [col_low 0.5], ...
                'MarkerFaceColor',  col_low, 'LineWidth', 1.2);
            plot(x_high, summary_high, ...
                '-', 'Color', [col_high 0.5], ...
                'MarkerFaceColor',  col_high, 'LineWidth', 1.2);

            % === Formatting ===
            tick_positions = 1:(nBands+1); % Simple integer positions
            band_labels = [plotNames(:); {'random'}];
            set(gca, 'XTick', tick_positions, 'XTickLabel', band_labels, ...
                'FontSize', 12, 'LineWidth', 1);
            ylim([-0.5 0.5])

            % Add separator lines between bands
            for b = 1:nBands
                line_pos = b + 0.5;   % midpoint between band b and b+1
                plot([line_pos line_pos], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            end

            % Create custom legend
            h_legend = zeros(2,1);
            h_legend(1) = plot(nan, nan, 's', 'Color', col_high, 'MarkerFaceColor', col_high, 'MarkerSize', 12);
            h_legend(2) = plot(nan, nan, 's', 'Color', col_low, 'MarkerFaceColor', col_low, 'MarkerSize', 12);
            legend(h_legend, {'High (Q4)', 'Low (Q1)'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12);

            box on;
            set(gcf, 'Units', 'normalized', 'Position', [0.515234375,0.152083333333333,0.389765625,0.772916666666667]);
            ax = gca;
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 1.4;
            fontsize(gcf,scale=1.3)
            fontname('Arial')

            xlim_range = xlim;
            % plot(xlim_range, [mu_rand mu_rand], ':', ...
            %     'Color', [0.7 0.7 0.7], ...
            %     'LineWidth', 1 ,'HandleVisibility', 'off');
            plot(xlim_range, [100 100], ':', ...
                'Color', [0.7 0.7 0.7], ...
                'LineWidth', 1 ,'HandleVisibility', 'off');

            xlabel(sprintf('%s bin', prestim_feature), 'FontSize', 16, 'FontWeight','bold');
            ylabel('MEP amplitude (log-transformed, z-scored)', 'FontSize', 16, 'FontWeight','bold');
            xlab = ax.XLabel;
            xlab.Position(2) = xlab.Position(2) - 0.02;
            ylab = ax.YLabel;
            ylab.Position(1) = ylab.Position(1) - 0.1;
            
            % === Add significance
            [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, 'alphabeta', reftep);
            sigBands_fields = fields(sigBands);

            if length(bands_of_interest) ~= 4
                allBands = {'theta','alpha','beta','gamma'};
                keepIdx = find(ismember(allBands, bands_of_interest));

                fn = fieldnames(sigBands);
                for i = 1:numel(fn)
                    key = fn{i};
                    v = sigBands.(key);
                    s = starText.(key);

                    if isempty(v)
                        sigBands.(key) = [];
                        starText.(key) = {};
                        continue
                    end

                    mask = ismember(v, keepIdx);
                    sigBands.(key) = v(mask);
                    starText.(key) = s(mask);
                end

                map = zeros(1, numel(allBands));
                map(keepIdx) = 1:numel(keepIdx);
            end

            for j = 1:size(sigBands_fields,1) % lowhigh lowrand highrand
                thisPair = sigBands_fields{j};
                sign_thisPair = sigBands.(thisPair);
                starText_thisPair = starText.(thisPair);

                if length(sign_thisPair) ~= length(starText_thisPair)
                    disp('Mismatch significance')
                    return
                end
                if isempty(sign_thisPair)
                    continue
                end

                for i = 1:length(sign_thisPair)
                    sign_now = sign_thisPair(i);
                    starText_now = starText_thisPair{i};

                    if numel(bands_of_interest) ~= 4
                        sign_now = map(sign_now);
                        if sign_now == 0
                            continue
                        end
                    end

                    x_rand = pos_rand;
                    if strcmpi(thisPair, 'lowhigh')
                        x_first = x_low(sign_now);
                        x_second =  x_high(sign_now);
                    elseif strcmpi(thisPair, 'lowrand')
                        x_first = x_low(sign_now);
                        x_second = x_rand;
                        yOffset = yOffset - 0.025;
                    else
                        x_first = x_high(sign_now);
                        x_second = x_rand;
                        yOffset = yOffset - 0.025;
                    end

                    % Draw significance
                    ylim_range = ylim;
                    tickLength = 0.005 * diff(ylim_range);
                    y_top = ylim_range(2) - yOffset;
                    y_bottom = y_top - tickLength;

                    line([x_first, x_first, x_second, x_second], ...
                        [y_bottom, y_top, y_top, y_bottom], ...
                        'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    text(((x_first + x_second) / 2), ylim_range(2) - yOffset + 0.005, starText_now, ...
                        'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', ...
                        'HandleVisibility', 'off');
                end
            end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''))), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), '-dpng', '-r300');
                end
            end
        end

        if strcmpi(plots_show,'MEAN_AMP') & length(bands_of_interest) == 1
            figure('Name','Relative dispersion across bands (all-trials)','Color','w'); hold on;
            
            jitter_amount = 0.02;
            % Single band case
            band_name = bands_of_interest{1};  % Get the single band
            band_table = T_all(strcmp(T_all.band, band_name), :);

            % Extract data
            data_low  = band_table.(['MEAN', '_low']);
            data_high = band_table.(['MEAN', '_high']);
            data_random = unique(T_all.(['MEAN', '_rand']));  % Random (same for all bands)

            % Positions - three separate x positions for Low, High, Random
            positions = [1, 2, 3];  % Each condition gets its own x-tick

            % Calculate mean and SEM
            mu_low   = mean(data_low,  'omitnan');
            sem_low  = std(data_low,  'omitnan') / sqrt(sum(~isnan(data_low)));

            mu_high  = mean(data_high, 'omitnan');
            sem_high = std(data_high, 'omitnan') / sqrt(sum(~isnan(data_high)));

            mu_rand  = mean(data_random, 'omitnan');
            sem_rand = std(data_random, 'omitnan') / sqrt(sum(~isnan(data_random)));

            % Plot error bars
            errorbar(positions(1), mu_low, sem_low, 'o', ...
                'Color', col_low, 'MarkerFaceColor', col_low, 'MarkerEdgeColor', col_low, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            hold on;
            errorbar(positions(2), mu_high, sem_high, 'o', ...
                'Color', col_high, 'MarkerFaceColor', col_high, 'MarkerEdgeColor', col_high, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);
            hold on;

            errorbar(positions(3), mu_rand, sem_rand, 'o', ...
                'Color', col_rand, 'MarkerFaceColor', col_rand, 'MarkerEdgeColor', col_rand, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            hold on;

            % Scatter points with jitter
            scatter(positions(1) + jitter_amount*(rand(size(data_low))-0.5), ...
                data_low, 25, ...
                'MarkerFaceColor', col_low, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
            scatter(positions(2) + jitter_amount*(rand(size(data_high))-0.5), ...
                data_high, 25, ...
                'MarkerFaceColor', col_high, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);
            scatter(positions(3) + jitter_amount*(rand(size(data_random))-0.5), ...
                data_random, 25, ...
                'MarkerFaceColor', col_rand, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);

            % === Formatting ===
            set(gca, 'XTick', positions, 'XTickLabel', {'Low', 'High', 'Random'}, ...
                'FontSize', 12, 'LineWidth', 1);
            ylim([-0.5 0.5])

            box on;
            set(gcf, 'Units', 'normalized', 'Position', [0.678125,0.152083333333333,0.226875,0.772916666666667]);
            ax = gca;
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 1.4;
            fontsize(gcf,scale=1.3)
            fontname('Arial')

            xlim_range = xlim;
            % plot(xlim_range, [mu_rand mu_rand], ':', ...
            %     'Color', [0.7 0.7 0.7], ...
            %     'LineWidth', 0.000000000001 ,'HandleVisibility', 'off');
            plot(xlim_range, [100 100], ':', ...
                'Color', [0.7 0.7 0.7], ...
                'LineWidth', 0.000000000001 ,'HandleVisibility', 'off');

            ax.Position = [0.18 0.13 0.72 0.77];
            xlim_range = xlim;
            % plot(xlim_range, [mu_rand mu_rand], ':', ...
            %     'Color', [0.7 0.7 0.7], ...
            %     'LineWidth', 1 ,'HandleVisibility', 'off');
            plot(xlim_range, [100 100], ':', ...
                'Color', [0.7 0.7 0.7], ...
                'LineWidth', 1 ,'HandleVisibility', 'off');

            xlabel(sprintf('%s bin', prestim_feature), 'FontSize', 16, 'FontWeight','bold');
            ylabel('MEP amplitude (log-transformed, z-scored)', 'FontSize', 16, 'FontWeight','bold');
            xlab = ax.XLabel;
            xlab = ax.XLabel;
            xlab.Position(2) = xlab.Position(2) - 0.02;
            ylab = ax.YLabel;
            ylab.Position(1) = ylab.Position(1) - 0.1;

            % === Add significance
            [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, strjoin(bands_of_interest, ''), reftep);
            sigBands_fields = fields(sigBands);

            if length(bands_of_interest) ~= 4
                allBands = {'theta','alpha','beta','gamma'};
                keepIdx = find(ismember(allBands, bands_of_interest));

                fn = fieldnames(sigBands);
                for i = 1:numel(fn)
                    key = fn{i};
                    v = sigBands.(key);
                    s = starText.(key);

                    if isempty(v)
                        sigBands.(key) = [];
                        starText.(key) = {};
                        continue
                    end

                    mask = ismember(v, keepIdx);
                    sigBands.(key) = v(mask);
                    starText.(key) = s(mask);
                end

                map = zeros(1, numel(allBands));
                map(keepIdx) = 1:numel(keepIdx);
            end

            for j = 1:size(sigBands_fields,1) % lowhigh lowrand highrand
                thisPair = sigBands_fields{j};
                sign_thisPair = sigBands.(thisPair);
                starText_thisPair = starText.(thisPair);

                if length(sign_thisPair) ~= length(starText_thisPair)
                    disp('Mismatch significance')
                    return
                end
                if isempty(sign_thisPair)
                    continue
                end

                for i = 1:length(sign_thisPair)
                    sign_now = sign_thisPair(i);
                    starText_now = starText_thisPair{i};

                    if strcmpi(thisPair, 'lowhigh')
                        x_first = positions(1);
                        x_second =  positions(2);
                    elseif strcmpi(thisPair, 'lowrand')
                        x_first = positions(1);
                        x_second = positions(3);
                        yOffset = yOffset - 2;
                    else
                        x_first = positions(2);
                        x_second = positions(3);
                        yOffset = yOffset - 2;
                    end

                    ylim_range = ylim;
                    tickLength = 0.005 * diff(ylim_range);
                    y_top = ylim_range(2) - yOffset;
                    y_bottom = y_top - tickLength;

                    line([x_first, x_first, x_second, x_second], ...
                        [y_bottom, y_top, y_top, y_bottom], ...
                        'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    text(((x_first + x_second) / 2), ylim_range(2) - yOffset + 0.005, starText_now, ...
                        'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', ...
                        'HandleVisibility', 'off');
                end
            end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''))), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), '-dpng', '-r300');
                end
            end
        end

    elseif strcmpi(prestim_feature, 'Phase')

        if nPhaseBins == 4
            nTotalBins = 5;         % including random bin
            phaseLabels = {'Trough','Rising','Peak','Falling','Random'};
        elseif nPhaseBins == 8
            nTotalBins = 9;         % including random bin
            phaseLabels = {'Trough','Early rising','Rising','Late rising', ...
                'Peak','Early falling','Falling','Late falling', ...
                'Random'};
        end

        % Define colours
        col_phasebins  = [150, 125, 160] / 255;
        col_rand = [120, 120, 130] / 255;  % Cool gray with slight blue tint

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if strcmpi(plots_show,'VARrel') || strcmpi(plots_show, 'VAR')

            figure('Name','Relative dispersion across phase bins','Color','w'); hold on;

            % === Plot phase bins boxplots ===
            position_offset = -0.25;
            band_data = cell(nPhaseBins+1,1);
            summary_phase = nan(nPhaseBins,1);
            x_phase = nan(nPhaseBins,1);
            jitter_amount = 0.02;
            for b = 1:nPhaseBins

                % Extract data for this phase bin
                data_phase = T_all.VARrel_phasebins(:,b);

                pos_thisBin = position_offset + 1;

                mu   = mean(data_phase,  'omitnan');
                sem  = std(data_phase,  'omitnan') / sqrt(sum(~isnan(data_phase)));

                errorbar(pos_thisBin, mu, sem, 'o', ...
                    'Color', col_phasebins, 'MarkerFaceColor', col_phasebins, 'MarkerEdgeColor', col_phasebins, ...
                    'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

                % Scatter dots
                scatter(pos_thisBin + jitter_amount*(rand(size(data_phase))-0.5), ...
                    data_phase, 25, ...
                    'MarkerFaceColor', col_phasebins, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, ...
                    'MarkerFaceAlpha', 0.6);

                % Store mean for connection line
                summary_phase(b) = mean(data_phase,'omitnan');
                x_phase(b) = pos_thisBin;

                band_data{b}.positions = pos_thisBin;
                position_offset = position_offset + 1;
            end

            % === Add random phase bin boxplot ===
            data_random = T_all.VARrel_rand;
            mu_rand   = mean(data_random,  'omitnan');
            sem_rand  = std(data_random,  'omitnan') / sqrt(sum(~isnan(data_random)));

            pos_rand = position_offset + 1;
            errorbar(pos_rand, mu_rand, sem_rand, 'o', ...
                'Color', col_rand, 'MarkerFaceColor', col_rand, 'MarkerEdgeColor', col_rand, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            % Add scatter for random
            scatter(pos_rand + jitter_amount*(rand(size(data_random))-0.5), ...
                data_random, 25, ...
                'MarkerFaceColor', col_rand, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, ...
                'MarkerFaceAlpha', 0.6);

            % === Add connection line between means ===
            plot(x_phase, summary_phase, ...
                '-', 'Color', col_phasebins, ...
                'MarkerFaceColor', col_phasebins, ...
                'LineWidth', 4);

            % 100% reference line
            yline(100, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

            % ===  Formatting ===
            xt = cellfun(@(x) x.positions, band_data(1:nPhaseBins), 'UniformOutput', true);
            xt = xt(:)';   % force row
            if nPhaseBins == 4
                set(gca, ...
                    'XTick', [xt, pos_rand], ...
                    'XTickLabel', phaseLabels, ...
                    'XTickLabelRotation', 0, ...
                    'FontSize', 12, ...
                    'LineWidth', 1.2);
            elseif nPhaseBins == 8
                set(gca, ...
                    'XTick', [xt, pos_rand], ...
                    'XTickLabel', phaseLabels, ...
                    'XTickLabelRotation',45, ...
                    'FontSize', 8, ...
                    'LineWidth', 1.2);
            end

            % Set x-axis limits with padding
            xlim([min(xt)-0.5, pos_rand+0.5]);  % <-- ADD THIS LINE
            xlim([0.25 5.25])
            
            if nPhaseBins == 4
                ylim([85 115])
            else
                ylim([70 130])
                xlim([0 pos_rand+0.5])
            end
            line_pos = pos_rand - 0.5;  % halfway in the gap
            plot([line_pos line_pos], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
            %plot(xlim, [mu_rand mu_rand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            %plot(xlim, [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

            box off;

            set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
            ax = gca;
            set(ax, 'FontSize', 14);
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 2;
            fontsize(gcf, scale=1.5)
            fontname('Arial')

            % xlabel('Phase', 'FontSize', 16, 'FontWeight','bold');
            % ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 16, 'FontWeight','bold');
            % title('Relative dispersion across phase bins');
            % xlab = ax.XLabel;
            % xlab.Position(2) = xlab.Position(2) - 1;
            % ylab = ax.YLabel;
            % ylab.Position(1) = ylab.Position(1) - 0.1;           

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, 'alpha')), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmpi(plots_show, 'MEAN_AMP')
              figure('Name','Mean amplitude across phase bins','Color','w'); hold on;

            % === Plot phase bins boxplots ===
            position_offset = -0.25;
            band_data = cell(nPhaseBins+1,1);
            summary_phase = nan(nPhaseBins,1);
            x_phase = nan(nPhaseBins,1);
            jitter_amount = 0.02;
            for b = 1:nPhaseBins

                % Extract data for this phase bin
                data_phase = T_all.MEAN_phasebins(:,b);

                pos_thisBin = position_offset + 1;

                mu   = mean(data_phase,  'omitnan');
                sem  = std(data_phase,  'omitnan') / sqrt(sum(~isnan(data_phase)));

                errorbar(pos_thisBin, mu, sem, 'o', ...
                    'Color', col_phasebins, 'MarkerFaceColor', col_phasebins, 'MarkerEdgeColor', col_phasebins, ...
                    'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

                % Scatter dots
                scatter(pos_thisBin + jitter_amount*(rand(size(data_phase))-0.5), ...
                    data_phase, 25, ...
                    'MarkerFaceColor', col_phasebins, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, ...
                    'MarkerFaceAlpha', 0.6);

                % Store mean for connection line
                summary_phase(b) = mean(data_phase,'omitnan');
                x_phase(b) = pos_thisBin;

                band_data{b}.positions = pos_thisBin;
                position_offset = position_offset + 1;
            end

            % === Add random phase bin boxplot ===
            data_random = T_all.MEAN_rand;
            mu_rand   = mean(data_random,  'omitnan');
            sem_rand  = std(data_random,  'omitnan') / sqrt(sum(~isnan(data_random)));

            pos_rand = position_offset + 1;
            errorbar(pos_rand, mu_rand, sem_rand, 'o', ...
                'Color', col_rand, 'MarkerFaceColor', col_rand, 'MarkerEdgeColor', col_rand, ...
                'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25);

            % Add scatter for random
            scatter(pos_rand + jitter_amount*(rand(size(data_random))-0.5), ...
                data_random, 25, ...
                'MarkerFaceColor', col_rand, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, ...
                'MarkerFaceAlpha', 0.6);

            % === Add connection line between means ===
            plot(x_phase, summary_phase, ...
                '-', 'Color', col_phasebins, ...
                'MarkerFaceColor', col_phasebins, ...
                'LineWidth', 4);

            % 100% reference line
            yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

            % ===  Formatting ===
            xt = cellfun(@(x) x.positions, band_data(1:nPhaseBins), 'UniformOutput', true);
            xt = xt(:)';   % force row
            if nPhaseBins == 4
                set(gca, ...
                    'XTick', [xt, pos_rand], ...
                    'XTickLabel', phaseLabels, ...
                    'XTickLabelRotation', 0, ...
                    'FontSize', 12, ...
                    'LineWidth', 1.2);
            elseif nPhaseBins == 8
                set(gca, ...
                    'XTick', [xt, pos_rand], ...
                    'XTickLabel', phaseLabels, ...
                    'XTickLabelRotation', 0, ...
                    'FontSize', 10, ...
                    'LineWidth', 1.2);
            end

            % Set x-axis limits with padding
            xlim([min(xt)-0.5, pos_rand+0.5]);  % <-- ADD THIS LINE
            xlim([0.25 5.25])

            line_pos = pos_rand - 0.5;  % halfway in the gap
            plot([line_pos line_pos], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
            %plot(xlim, [mu_rand mu_rand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            %plot(xlim, [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

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

            % xlabel('Phase', 'FontSize', 16, 'FontWeight','bold');
            % ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 16, 'FontWeight','bold');
            % title('Relative dispersion across phase bins');
            % xlab = ax.XLabel;
            % xlab.Position(2) = xlab.Position(2) - 1;
            % ylab = ax.YLabel;
            % ylab.Position(1) = ylab.Position(1) - 0.1;     

            % === Add significance
            if nPhaseBins == 4
                [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, 'alpha', reftep);

                x_phase = [x_phase; pos_rand];
                for n_sig = 1:size(sigBands, 1)
                    idx_pair1 = sigBands(n_sig, 1);
                    idx_pair2 = sigBands(n_sig, 2);

                    % Draw significance
                    ylim_range = ylim;
                    tickLength = 0.005 * diff(ylim_range);
                    y_top = ylim_range(2) - yOffset;
                    y_bottom = y_top - tickLength;

                    x_first = x_phase(idx_pair1);
                    x_second = x_phase(idx_pair2);
                    starText_now = starText(1);

                    line([x_first, x_first, x_second, x_second], ...
                        [y_bottom, y_top, y_top, y_bottom], ...
                        'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
                    text(((x_first + x_second) / 2), ylim_range(2) - yOffset + 0.0017, starText_now, ...
                        'HorizontalAlignment', 'center', 'FontSize', 28, 'FontWeight', 'bold', ...
                        'HandleVisibility', 'off');

                    yOffset = yOffset - 0.009;

                end
            end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, 'alpha')), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Phase_split
            insig_phase_subjects = load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\reftep_insig_phase_subject_list.mat').reftep_insig_phase_subject_list;
            sig_phase_subjects = load('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\extracted_ep_info\reftep_sig_phase_subject_list.mat').reftep_sig_phase_subject_list;

            insig_phase_subject_idx = find(ismember(all_subjects, insig_phase_subjects));
            sig_phase_subject_idx   = find(ismember(all_subjects, sig_phase_subjects));

            figure; hold on;

            groups      = {sig_phase_subject_idx, insig_phase_subject_idx};
            group_names = {'Significant phase effect', 'Insignificant phase effect'};
            group_colors = {col_phasebins, col_phasebins*0.2};
            group_offset = [-0.2, 0.2];   % left/right offset within each bin
            jitter_amount = 0.05;

            for g = 1:2
                subj_idx  = groups{g};
                col_group = group_colors{g};

                summary_phase = nan(nPhaseBins, 1);
                x_phase       = nan(nPhaseBins, 1);

                % === Phase bins ===
                for b = 1:nPhaseBins
                    data_phase  = T_all.VARrel_phasebins(subj_idx, b);
                    pos_thisBin = b + group_offset(g);

                    mu  = mean(data_phase, 'omitnan');
                    sem = std(data_phase,  'omitnan') / sqrt(sum(~isnan(data_phase)));

                    errorbar(pos_thisBin, mu, sem, 'o', ...
                        'Color', col_group, 'MarkerFaceColor', col_group, 'MarkerEdgeColor', col_group, ...
                        'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25,'HandleVisibility', 'off');

                    scatter(pos_thisBin + jitter_amount*(rand(size(data_phase))-0.5), ...
                        data_phase, 25, ...
                        'MarkerFaceColor', col_group, 'MarkerEdgeColor', 'k', ...
                        'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');

                    summary_phase(b) = mu;
                    x_phase(b)       = pos_thisBin;
                end

                % === Connection line ===
                plot(x_phase, summary_phase, '-', ...
                    'Color', col_group, 'LineWidth', 4, 'DisplayName', group_names{g});

                % === Random bin ===
                data_random = T_all.VARrel_rand(subj_idx);
                mu_rand  = mean(data_random, 'omitnan');
                sem_rand = std(data_random,  'omitnan') / sqrt(sum(~isnan(data_random)));
                pos_rand = nPhaseBins + 1 + group_offset(g);

                errorbar(pos_rand, mu_rand, sem_rand, 'o', ...
                    'Color', col_rand, 'MarkerFaceColor', col_rand, 'MarkerEdgeColor', col_rand, ...
                    'MarkerSize', 10, 'LineWidth', 4, 'CapSize', 25, 'HandleVisibility', 'off');

                scatter(pos_rand + jitter_amount*(rand(size(data_random))-0.5), ...
                    data_random, 25, ...
                    'MarkerFaceColor', col_rand, 'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
            end

            % === Legend ===
            legend(group_names, 'Location', 'northwest', 'Box', 'off');

            % === Reference line ===
            yline(100, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom', 'HandleVisibility', 'off');

            % ===  Formatting ===
            xt = cellfun(@(x) x.positions, band_data(1:nPhaseBins), 'UniformOutput', true);
            xt = xt(:)';   % force row
            if nPhaseBins == 4
                set(gca, ...
                    'XTick', [xt, pos_rand], ...
                    'XTickLabel', phaseLabels, ...
                    'XTickLabelRotation', 0, ...
                    'FontSize', 12, ...
                    'LineWidth', 1.2);
            elseif nPhaseBins == 8
                set(gca, ...
                    'XTick', [xt, pos_rand], ...
                    'XTickLabel', phaseLabels, ...
                    'XTickLabelRotation',45, ...
                    'FontSize', 8, ...
                    'LineWidth', 1.2);
            end

            % Set x-axis limits with padding
            xlim([min(xt)-0.5, pos_rand+0.5]);  % <-- ADD THIS LINE
            xlim([0.25 5.75])

            if nPhaseBins == 4
                ylim([85 115])
            else
                ylim([70 130])
                xlim([0 pos_rand+0.5])
            end
            line_pos = pos_rand - 0.70;  % halfway in the gap
            plot([line_pos line_pos], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'HandleVisibility', 'off');
            %plot(xlim, [mu_rand mu_rand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            %plot(xlim, [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

            box off;

            set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
            ax = gca;
            set(ax, 'FontSize', 14);
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 2;
            fontsize(gcf, scale=1.5)
            fontname('Arial')
        end

        if individual_subjects
            if strcmpi(plots_show,'VARrel') || strcmpi(plots_show, 'VAR')

                nSubjects = 29;
                colors = parula(nSubjects);  % or lines(nSubjects), jet(nSubjects), turbo(nSubjects), etc.


                figure('Name','Relative dispersion across phase bins','Color','w'); hold on;

                % === Plot phase bins boxplots ===
                position_offset = -0.25;
                band_data = cell(nPhaseBins+1,1);
                summary_phase = nan(nPhaseBins,1);
                x_phase = nan(nPhaseBins,1);
                jitter_amount = 0.02;
                for blaba = 1:29
                    mu = T_all.VARrel_phasebins(blaba,:);
                    data_random = T_all.VARrel_rand(blaba);
                    databla = [mu data_random];

                    x = [1, 2, 3, 4, 5] + position_offset;

                    plot(x, databla, '-o', ...
                        'Color', colors(blaba,:), ...
                        'MarkerFaceColor', colors(blaba,:), ...
                        'MarkerEdgeColor', colors(blaba,:), ...
                        'MarkerSize', 10, ...
                        'LineWidth', 3);
                end

                mu_hat = mean(T_all.VARrel_phasebins, 1, 'omitnan');
                plot(x(1:4), mu_hat, '-o', ...
                    'Color', 'k', ...
                    'LineWidth', 4, 'MarkerFaceColor', 'k', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 10);

                % 100% reference line
                yline(100, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

                % ===  Formatting ===
                %xt = cellfun(@(x) x.positions, band_data(1:nPhaseBins), 'UniformOutput', true);
                %xt = xt(:)';   % force row
                if nPhaseBins == 4
                    set(gca, ...
                        'XTick', [xt, pos_rand], ...
                        'XTickLabel', phaseLabels, ...
                        'XTickLabelRotation', 0, ...
                        'FontSize', 12, ...
                        'LineWidth', 1.2);
                elseif nPhaseBins == 8
                    set(gca, ...
                        'XTick', [xt, pos_rand], ...
                        'XTickLabel', phaseLabels, ...
                        'XTickLabelRotation',45, ...
                        'FontSize', 8, ...
                        'LineWidth', 1.2);
                end

                % Set x-axis limits with padding
                xlim([min(xt)-0.5, pos_rand+0.5]);  % <-- ADD THIS LINE
                xlim([0.25 5.25])

                if nPhaseBins == 4
                    ylim([80 120])
                else
                    ylim([70 130])
                    xlim([0 pos_rand+0.5])
                end
                line_pos = pos_rand - 0.5;  % halfway in the gap
                plot([line_pos line_pos], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
                %plot(xlim, [mu_rand mu_rand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
                %plot(xlim, [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

                box off;

                set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
                ax = gca;
                set(ax, 'FontSize', 14);
                ax.Box = 'off';
                ax.TickDir = 'out';
                ax.LineWidth = 2;
                fontsize(gcf, scale=1.5)
                fontname('Arial')

                % xlabel('Phase', 'FontSize', 16, 'FontWeight','bold');
                % ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 16, 'FontWeight','bold');
                % title('Relative dispersion across phase bins');
                % xlab = ax.XLabel;
                % xlab.Position(2) = xlab.Position(2) - 1;
                % ylab = ax.YLabel;
                % ylab.Position(1) = ylab.Position(1) - 0.1;

                % Save as high-quality PNG (300 dpi)
                if save_output
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_subjects.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmpi(plots_show, 'MEAN_AMP')
               
                nSubjects = 29;
                colors = parula(nSubjects);  % or lines(nSubjects), jet(nSubjects), turbo(nSubjects), etc.
              
                figure('Name','Relative dispersion across phase bins','Color','w'); hold on;

                % === Plot phase bins boxplots ===
                position_offset = -0.25;
                band_data = cell(nPhaseBins+1,1);
                summary_phase = nan(nPhaseBins,1);
                x_phase = nan(nPhaseBins,1);
                jitter_amount = 0.02;
                for blaba = 1:29
                    mu = T_all.MEAN_phasebins(blaba,:);
                    data_random = T_all.MEAN_rand(blaba);
                    databla = [mu data_random];

                    x = [1, 2, 3, 4, 5] + position_offset;

                    plot(x, databla, '-o', ...
                        'Color', colors(blaba,:), ...
                        'MarkerFaceColor', colors(blaba,:), ...
                        'MarkerEdgeColor', colors(blaba,:), ...
                        'MarkerSize', 10, ...
                        'LineWidth', 3);
                end

                mu_hat = mean(T_all.MEAN_phasebins, 1, 'omitnan');
                %mu_hat = [mu_hat mean(T_all.MEAN_rand, 1, 'omitnan')]
                plot(x(1:4), mu_hat, '-o', ...
                    'Color', 'k', ...
                    'LineWidth', 4, 'MarkerFaceColor', 'k', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 10);

                % 100% reference line
                yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

                % ===  Formatting ===
                %xt = cellfun(@(x) x.positions, band_data(1:nPhaseBins), 'UniformOutput', true);
                %xt = xt(:)';   % force row
                if nPhaseBins == 4
                    set(gca, ...
                        'XTick', [xt, pos_rand], ...
                        'XTickLabel', phaseLabels, ...
                        'XTickLabelRotation', 0, ...
                        'FontSize', 12, ...
                        'LineWidth', 1.2);
                elseif nPhaseBins == 8
                    set(gca, ...
                        'XTick', [xt, pos_rand], ...
                        'XTickLabel', phaseLabels, ...
                        'XTickLabelRotation',45, ...
                        'FontSize', 8, ...
                        'LineWidth', 1.2);
                end

                % Set x-axis limits with padding
                xlim([min(xt)-0.5, pos_rand+0.5]);  % <-- ADD THIS LINE
                xlim([0.25 5.25])

                ylim([-0.4 0.4])
                line_pos = pos_rand - 0.5;  % halfway in the gap
                plot([line_pos line_pos], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
                %plot(xlim, [mu_rand mu_rand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
                %plot(xlim, [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

                box off;

                set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
                ax = gca;
                set(ax, 'FontSize', 14);
                ax.Box = 'off';
                ax.TickDir = 'out';
                ax.LineWidth = 2;
                fontsize(gcf, scale=1.5)
                fontname('Arial')

                % xlabel('Phase', 'FontSize', 16, 'FontWeight','bold');
                % ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 16, 'FontWeight','bold');
                % title('Relative dispersion across phase bins');
                % xlab = ax.XLabel;
                % xlab.Position(2) = xlab.Position(2) - 1;
                % ylab = ax.YLabel;
                % ylab.Position(1) = ylab.Position(1) - 0.1;

                % Save as high-quality PNG (300 dpi)
                if save_output
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s_subjects.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end
        end

    elseif strcmpi(prestim_feature, 'PowerPhase')

        comboLabels = {'Trough','Rising','Peak','Falling', ...
            'Trough','Rising','Peak','Falling', ...
            'random'};
        % Define plot colours
        col_combos = [150, 125, 160] / 255;
        col_rand   = [120, 120, 130] / 255;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmpi(plots_show,'VARrel') || strcmpi(plots_show, 'VAR')

            figure('Name','Relative dispersion across Power × Phase combos','Color','w');
            hold on;

            % x positions for the 4 phase bins only
            xpos = 1:4;

            % Assuming the order in G.mean_VARrel_combos is:
            % [Low-Trough, Low-Rising, Low-Peak, Low-Falling, ...
            %  High-Trough, High-Rising, High-Peak, High-Falling]

            low_mean  = G.mean_VARrel_combos(1:4);
            low_sem   = G.sem_VARrel_combos(1:4);

            high_mean = G.mean_VARrel_combos(5:8);
            high_sem  = G.sem_VARrel_combos(5:8);

            % choose different colors for low and high
            col_high = [253, 70, 101] / 255;    % Vibrant red
            col_low = [59, 82, 139] / 255;    % Medium blue

            % Low power
            h1 = errorbar(xpos, low_mean, low_sem, '-o', ...
                'Color', col_low, ...
                'MarkerFaceColor', col_low, ...
                'MarkerEdgeColor', col_low, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

            % High power
            h2 = errorbar(xpos, high_mean, high_sem, '-o', ...
                'Color', col_high, ...
                'MarkerFaceColor', col_high, ...
                'MarkerEdgeColor', col_high, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

               % 100% reference line
            yline(100, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer','bottom');

            ylim([85 115]);
            xlim([0.25 4.75]);

            set(gca, ...
                'XTick', xpos, ...
                'XTickLabel', {'Trough','Rising','Peak','Falling'}, ...
                'FontSize', 12, ...
                'LineWidth', 1.4);

        
            %xlabel('Phase', 'FontSize', 14, 'FontWeight','bold');
            %ylabel(sprintf('Normalized within-state %s (%%)', var_method), ...
            %    'FontSize', 14, 'FontWeight','bold');

            legend([h2 h1], {'High mu power','Low mu power','Random'}, ...
                'Location', 'northwest', ...
                'Box', 'off', 'FontSize', 14);

            set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
            ax = gca;
            set(ax, 'FontSize', 14);
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 2;
            fontsize(gcf, scale=1.5)
            fontname('Arial')

      




            % figure('Name','Relative dispersion across Power × Phase combos','Color','w'); hold on;
            % 
            % combo_step = 0.25;
            % group_gap  = 0.35;
            % start_x    = 0.25;
            % 
            % xpos = zeros(1, 9);
            % xpos(1:4) = start_x + (0:3) * combo_step;
            % xpos(5:8) = xpos(4) + group_gap + (1:4) * combo_step;
            % xpos(9)   = xpos(8) + group_gap + combo_step;
            % 
            % % Mean ± SEM for combos
            % errorbar(xpos(1:8), ...
            %     G.mean_VARrel_combos(1:8), ...
            %     G.sem_VARrel_combos(1:8), ...
            %     'o', ...
            %     'Color', col_combos, ...
            %     'MarkerFaceColor', col_combos, ...
            %     'MarkerEdgeColor', col_combos, ...
            %     'MarkerSize', 10, ...
            %     'LineWidth', 4, ...
            %     'CapSize', 25);
            % 
            % % Random
            % errorbar(xpos(9), ...
            %     G.mean_VARrel_rand, ...
            %     G.sem_VARrel_rand, ...
            %     'o', ...
            %     'Color', col_rand, ...
            %     'MarkerFaceColor', col_rand, ...
            %     'MarkerEdgeColor', col_rand, ...
            %     'MarkerSize', 10, ...
            %     'LineWidth', 4, ...
            %     'CapSize', 25);
            % 
            % % === Add connection line between means ===
            % plot(xpos(1:8), G.mean_VARrel_combos(1:8), ...
            %     '-', 'Color', col_combos, ...
            %     'MarkerFaceColor', col_combos, ...
            %     'LineWidth', 1.4);
            % 
            % 
            % set(gca, ...
            %     'XTick', xpos, ...
            %     'XTickLabel', comboLabels, ...
            %     'XTickLabelRotation', 0, ...
            %     'FontSize', 10, ...
            %     'LineWidth', 1.2);
            % 
            % % % Add custom text labels
            % % xticks_range = xticks;
            % % y_offset = ax.YLim(1) - 0.065* diff(ax.YLim); % Adjust vertical position
            % % text(((xticks_range(1)+xticks_range(4))/2), y_offset, 'Low Power', ...
            % %     'HorizontalAlignment', 'center', ...
            % %     'Rotation', 0, ...
            % %     'FontSize', 14);
            % % text(((xticks_range(5)+xticks_range(8))/2), y_offset, 'High Power', ...
            % %     'HorizontalAlignment', 'center', ...
            % %     'Rotation', 0, ...
            % %     'FontSize', 14);
            % 
            % yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % 
            % ylim([85 115]);
            % xlim([0 xpos(9) + 0.25]); % Added padding on right side
            % 
            % set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
            % ax = gca;
            % set(ax, 'FontSize', 12);
            % ax.Box = 'off';
            % ax.TickDir = 'out';
            % ax.LineWidth = 1.6;
            % fontsize(gcf, scale=1.4)
            % fontname('Arial')           
            % 
            % %xlabel('Power × Phase combination', 'FontSize', 14, 'FontWeight','bold');
            % %ylabel(sprintf('Normalized within-state %s (%%)', var_method), 'FontSize', 14, 'FontWeight','bold');
            % %xlab = ax.XLabel;
            % %xlab.Position(2) = xlab.Position(2) - 2.5;
            % %ylab = ax.YLabel;
            % %ylab.Position(1) = ylab.Position(1) - 0.1;
            % 
            % % Add vertical separators
            % v1 = mean([xpos(4), xpos(5)]);   % Low ↔ High
            % v2 = mean([xpos(8), xpos(9)]);   % High ↔ Random
            % plot([v1 v1], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % plot([v2 v2], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % %plot(xlim, [G.mean_VARrel_rand G.mean_VARrel_rand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % %plot(xlim, [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % 
            % fontsize(gcf,scale=1.3)
            % fontname('Arial')

            % % === Add significance
            % [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, 'alpha', reftep);
            % 
            % for thisPair = 1:size(sigBands.pair,2)
            %     % Trough Low - High
            %     % Rising Low - High
            %     % Peak Low - High
            %     % Falling Low - High
            %     idx_pair1 = sigBands.pair(thisPair);
            %     idx_pair2 = sigBands.pair(thisPair) + 4;
            % 
            % 
            %     % Draw significance
            %     ylim_range = ylim;
            %     tickLength = 0.005 * diff(ylim_range);
            %     y_top = ylim_range(2) - yOffset;
            %     y_bottom = y_top - tickLength;
            % 
            %     x_first = xpos(idx_pair1);
            %     x_second = xpos(idx_pair2);
            %     starText_now = starText.asterix(thisPair);
            % 
            %     line([x_first, x_first, x_second, x_second], ...
            %         [y_bottom, y_top, y_top, y_bottom], ...
            %         'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            %     text(((x_first + x_second) / 2), ylim_range(2) - yOffset + 0.5, starText_now, ...
            %         'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', ...
            %         'HandleVisibility', 'off');
            % 
            %     yOffset = yOffset - 2;
            % end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, 'alpha')), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % -------- Plot 2: Mean amplitudes across PowerPhase combos --------
        if strcmpi(plots_show, 'MEAN_AMP')
           
            figure('Name','Relative dispersion across Power × Phase combos','Color','w');
            hold on;

            % x positions for the 4 phase bins only
            xpos = 1:4;

            % Assuming the order in G.mean_VARrel_combos is:
            % [Low-Trough, Low-Rising, Low-Peak, Low-Falling, ...
            %  High-Trough, High-Rising, High-Peak, High-Falling]

            low_mean  = G.mean_MEAN_combos(1:4);
            low_sem   = G.sem_MEAN_combos(1:4);

            high_mean = G.mean_MEAN_combos(5:8);
            high_sem  = G.sem_MEAN_combos(5:8);

            % choose different colors for low and high
            col_high = [253, 70, 101] / 255;    % Vibrant red
            col_low = [59, 82, 139] / 255;    % Medium blue

            % Low power
            h1 = errorbar(xpos, low_mean, low_sem, '-o', ...
                'Color', col_low, ...
                'MarkerFaceColor', col_low, ...
                'MarkerEdgeColor', col_low, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

            % High power
            h2 = errorbar(xpos, high_mean, high_sem, '-o', ...
                'Color', col_high, ...
                'MarkerFaceColor', col_high, ...
                'MarkerEdgeColor', col_high, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

               % 100% reference line
            yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

            ylim([-0.2 0.2]);
            xlim([0.25 4.75]);

            set(gca, ...
                'XTick', xpos, ...
                'XTickLabel', {'Trough','Rising','Peak','Falling'}, ...
                'FontSize', 12, ...
                'LineWidth', 1.4);

        
            %xlabel('Phase', 'FontSize', 14, 'FontWeight','bold');
            %ylabel(sprintf('Normalized within-state %s (%%)', var_method), ...
            %    'FontSize', 14, 'FontWeight','bold');

            legend([h2 h1], {'High mu power','Low mu power','Random'}, ...
                'Location', 'northeast', ...
                'Box', 'off', 'FontSize', 14);

            set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
            ax = gca;
            set(ax, 'FontSize', 14);
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 2;
            fontsize(gcf, scale=1.5)
            fontname('Arial')

            % figure('Name','Mean amplitude across Power × Phase combos','Color','w'); hold on;

            % combo_step = 0.25;
            % group_gap  = 0.35;
            % start_x    = 0.25;
            % 
            % xpos = zeros(1, 9);
            % xpos(1:4) = start_x + (0:3) * combo_step;
            % xpos(5:8) = xpos(4) + group_gap + (1:4) * combo_step;
            % xpos(9)   = xpos(8) + group_gap + combo_step;
            % 
            % % Mean ± SEM for combos
            % errorbar(xpos(1:8), ...
            %     G.mean_MEAN_combos(1:8), ...
            %     G.sem_MEAN_combos(1:8), ...
            %     'o', ...
            %     'Color', col_combos, ...
            %     'MarkerFaceColor', col_combos, ...
            %     'MarkerEdgeColor', col_combos, ...
            %     'MarkerSize', 10, ...
            %     'LineWidth', 4, ...
            %     'CapSize', 25);
            % 
            % % Random
            % errorbar(xpos(9), ...
            %     G.mean_MEAN_rand, ...
            %     G.sem_MEAN_rand, ...
            %     'o', ...
            %     'Color', col_rand, ...
            %     'MarkerFaceColor', col_rand, ...
            %     'MarkerEdgeColor', col_rand, ...
            %     'MarkerSize', 10, ...
            %     'LineWidth', 4, ...
            %     'CapSize', 25);
            % 
            % % === Add connection line between means ===
            % plot(xpos(1:8), G.mean_MEAN_combos(1:8), ...
            %     '-', 'Color', col_combos, ...
            %     'MarkerFaceColor', col_combos, ...
            %     'LineWidth', 1.4);
            % 
            % ylim([-0.5 0.5]);
            % xlim([0 xpos(9) + 0.25]); % Added padding on right side
            % 
            % set(gca, ...
            %     'XTick', xpos, ...
            %     'XTickLabel', comboLabels, ...
            %     'XTickLabelRotation', 0, ...
            %     'FontSize', 10, ...
            %     'LineWidth', 1.2);
            % box on;
            % set(gcf, 'Units', 'normalized', 'Position', [0.515234375,0.152083333333333,0.389765625,0.772916666666667]);
            % ax = gca;
            % ax.Box = 'off';
            % ax.TickDir = 'out';
            % ax.LineWidth = 1.4;
            % 
            % % Add custom text labels
            % xticks_range = xticks;
            % y_offset = ax.YLim(1) - 0.065* diff(ax.YLim); % Adjust vertical position
            % text(((xticks_range(1)+xticks_range(4))/2), y_offset, 'Low Power', ...
            %     'HorizontalAlignment', 'center', ...
            %     'Rotation', 0, ...
            %     'FontSize', 14);
            % text(((xticks_range(5)+xticks_range(8))/2), y_offset, 'High Power', ...
            %     'HorizontalAlignment', 'center', ...
            %     'Rotation', 0, ...
            %     'FontSize', 14);
            % 
            % xlabel('Power × Phase combination', 'FontSize', 16, 'FontWeight','bold');
            % ylabel('MEP amplitude (log-transformed, z-scored', 'FontSize', 16, 'FontWeight','bold');
            % 
            % % Adjust label positions
            % xlab = ax.XLabel;
            % xlab.Position(2) = xlab.Position(2) - 2.5;
            % ylab = ax.YLabel;
            % ylab.Position(1) = ylab.Position(1) - 0.1;
            % 
            % % Add vertical separators
            % v1 = mean([xpos(4), xpos(5)]);   % Low ↔ High
            % v2 = mean([xpos(8), xpos(9)]);   % High ↔ Random
            % plot([v1 v1], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % plot([v2 v2], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % % plot(xlim, [G.mean_MEAN_rand G.mean_MEAN_rand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % plot(xlim, [0 0], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % 
            % fontsize(gcf,scale=1.3)
            % fontname('Arial')

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, 'alpha')), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end
        end

    elseif strcmpi(prestim_feature, 'PowerConnectivity')

        bands      = string(G.band);   % e.g., "theta","alpha","beta","gamma"
        nBands     = height(G);
        nCombos    = size(G.mean_MEAN_combos, 2);
        nTotalBins = nCombos + 1;
        comboLabels = {'Low FC','High FC'};

        % Define plot colours
        col_combos = [150, 125, 160] / 255;
        col_rand   = [120, 120, 130] / 255;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -------- Plot 1: Relative dispersion across Power x Connectivity combos --------
        if strcmpi(plots_show,'VARrel') || strcmpi(plots_show, 'VAR')

            figure('Name','Relative dispersion across Power × Phase combos','Color','w');
            hold on;

            % x positions for the 4 phase bins only
            xpos = 1:2;

            % Assuming the order in G.mean_VARrel_combos is:
            % [Low-Trough, Low-Rising, Low-Peak, Low-Falling, ...
            %  High-Trough, High-Rising, High-Peak, High-Falling]

            low_mean  = G.mean_VARrel_combos(1:2);
            low_sem   = G.sem_VARrel_combos(1:2);

            high_mean = G.mean_VARrel_combos(3:4);
            high_sem  = G.sem_VARrel_combos(3:4);

            % choose different colors for low and high
            col_high = [253, 70, 101] / 255;    % Vibrant red
            col_low = [59, 82, 139] / 255;    % Medium blue

            % Low power
            h1 = errorbar(xpos, low_mean, low_sem, '-o', ...
                'Color', col_low, ...
                'MarkerFaceColor', col_low, ...
                'MarkerEdgeColor', col_low, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

            % High power
            h2 = errorbar(xpos, high_mean, high_sem, '-o', ...
                'Color', col_high, ...
                'MarkerFaceColor', col_high, ...
                'MarkerEdgeColor', col_high, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

               % 100% reference line
            yline(100, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

            ylim([85 115]);
            xlim([0.5 2.50]);

            set(gca, ...
                'XTick', xpos, ...
                'XTickLabel', {'Low','High'}, ...
                'FontSize', 12, ...
                'LineWidth', 1.4);

        
            %xlabel('Phase', 'FontSize', 14, 'FontWeight','bold');
            %ylabel(sprintf('Normalized within-state %s (%%)', var_method), ...
            %    'FontSize', 14, 'FontWeight','bold');

            legend([h2 h1], {'High mu power','Low mu power','Random'}, ...
                'Location', 'northwest', ...
                'Box', 'off', 'FontSize', 14);

            set(gcf, 'Units', 'normalized', 'Position', [0.666796875,0.152083333333333,0.238203125,0.772916666666667]);
            ax = gca;
            set(ax, 'FontSize', 14);
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 2;
            fontsize(gcf, scale=1.5)
            fontname('Arial')  

            % % Define x-positions
            % combo_step = 0.25;
            % group_gap  = 0.35;
            % start_x    = 0.25;
            % xpos = zeros(1, nTotalBins);
            % 
            % % Position combos (adjust indices based on your combo structure)
            % nCombosPerGroup = nCombos / 2; % Assuming two groups (Low/High)
            % xpos(1:nCombosPerGroup) = start_x + (0:(nCombosPerGroup-1)) * combo_step;
            % xpos((nCombosPerGroup+1):nCombos) = xpos(nCombosPerGroup) + group_gap + (1:nCombosPerGroup) * combo_step;
            % xpos(nTotalBins) = xpos(nCombos) + group_gap + combo_step; % Random position
            % 
            % for iB = 1:nBands
            %     figure('Name','Relative dispersion across Power × Connectivity combos','Color','w'); hold on;
            % 
            %     % Get data for this band
            %     yComb = G.mean_VARrel_combos(iB, 1:nCombos);
            %     eComb = G.sem_VARrel_combos(iB, 1:nCombos);
            %     yRand = G.mean_VARrel_rand(iB);
            %     eRand = G.sem_VARrel_rand(iB);
            % 
            %     % Mean ± SEM for combos
            %     errorbar(xpos(1:nCombos), yComb, eComb, ...
            %         'o', ...
            %         'Color', col_combos, ...
            %         'MarkerFaceColor', col_combos, ...
            %         'MarkerEdgeColor', col_combos, ...
            %         'MarkerSize', 10, ...
            %         'LineWidth', 4, ...
            %         'CapSize', 25);
            % 
            %     % Random
            %     errorbar(xpos(nTotalBins), yRand, eRand, ...
            %         'o', ...
            %         'Color', col_rand, ...
            %         'MarkerFaceColor', col_rand, ...
            %         'MarkerEdgeColor', col_rand, ...
            %         'MarkerSize', 10, ...
            %         'LineWidth', 4, ...
            %         'CapSize', 25);
            % 
            %     % Add connection line between combo means
            %     plot(xpos(1:nCombos), yComb, '-', ...
            %         'Color', col_combos, ...
            %         'LineWidth', 1.4);
            % 
            %     % Axis settings
            %     ylim([80 120]);
            %     xlim([0 xpos(nTotalBins) + 0.25]);
            %     set(gca, ...
            %         'XTick', xpos, ...
            %         'XTickLabel', comboLabels, ...
            %         'XTickLabelRotation', 0, ...
            %         'FontSize', 10, ...
            %         'LineWidth', 1.2, ...
            %         'Box', 'off', ...
            %         'TickDir', 'out');
            %     set(gcf, 'Units', 'normalized', ...
            %         'Position', [0.515234375, 0.152083333333333, 0.389765625, 0.772916666666667]);
            % 
            %     ax = gca;
            %     ax.LineWidth = 1.4;
            %     % Add custom group labels (adjust based on your grouping)
            %     xticks_range = xticks;
            %     y_offset = ax.YLim(1) - 0.065 * diff(ax.YLim);
            %     text(mean(xticks_range(1:nCombosPerGroup)), y_offset, 'Low Power', ...
            %         'HorizontalAlignment', 'center', ...
            %         'Rotation', 0, ...
            %         'FontSize', 14);
            %     text(mean(xticks_range((nCombosPerGroup+1):nCombos)), y_offset, 'High Power', ...
            %         'HorizontalAlignment', 'center', ...
            %         'Rotation', 0, ...
            %         'FontSize', 14);
            % 
            %     % Axis labels
            %     xlabel('Power × Connectivity combination', 'FontSize', 16, 'FontWeight','bold');
            %     ylabel(sprintf('Normalized within-state %s', var_method)', 'FontSize', 16, 'FontWeight','bold');
            %     %title(sprintf('%s', bands(iB)), 'FontSize', 16);
            % 
            %     % Adjust label positions
            %     xlab = ax.XLabel;
            %     xlab.Position(2) = xlab.Position(2) - 2.5;
            %     ylab = ax.YLabel;
            %     ylab.Position(1) = ylab.Position(1) - 0.1;
            % 
            %     % Add vertical separators
            %     v1 = mean([xpos(nCombosPerGroup), xpos(nCombosPerGroup+1)]);   % Low ↔ High
            %     v2 = mean([xpos(nCombos), xpos(nTotalBins)]);                   % High ↔ Random
            %     plot([v1 v1], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            %     plot([v2 v2], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            %     %plot(xlim, [yRand yRand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            %     plot(xlim, [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % 
            %     fontsize(gcf, scale=1.3);
            %     fontname('Arial');
            % end

            % % === Add significance
            % [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, 'alpha', reftep);
            % 
            % for thisPair = 1:size(sigBands.pair,2)
            %     % Low Low - High Low
            %     % Low Low - Low High
            %     % Low Low - High High
            %     % High Low - Low High
            %     % High Low - High High
            %     % Low High - High High
            %     switch ~isempty(sigBands)
            %         case sigBands.pair(thisPair) == 1
            %             idx_pair1 = 1;
            %             idx_pair2 = 3;
            %         case sigBands.pair(thisPair) == 2
            %             idx_pair1 = 1;
            %             idx_pair2 = 2;
            %         case sigBands.pair(thisPair) == 3
            %             idx_pair1 = 1;
            %             idx_pair2 = 4;
            %         case sigBands.pair(thisPair) == 4
            %             idx_pair1 = 3;
            %             idx_pair2 = 2;
            %         case sigBands.pair(thisPair) == 5
            %             idx_pair1 = 3;
            %             idx_pair2 = 4;
            %         case sigBands.pair(thisPair) == 6
            %             idx_pair1 = 2;
            %             idx_pair2 = 4;
            %     end
            % 
            %     % Draw significance
            %     ylim_range = ylim;
            %     tickLength = 0.005 * diff(ylim_range);
            %     y_top = ylim_range(2) - yOffset;
            %     y_bottom = y_top - tickLength;
            % 
            %     x_first = xpos(idx_pair1);
            %     x_second = xpos(idx_pair2);
            %     starText_now = starText.asterix(thisPair);
            % 
            %     line([x_first, x_first, x_second, x_second], ...
            %         [y_bottom, y_top, y_top, y_bottom], ...
            %         'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            %     text(((x_first + x_second) / 2), ylim_range(2) - yOffset + 0.5, starText_now, ...
            %         'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', ...
            %         'HandleVisibility', 'off');
            % 
            %     yOffset = yOffset - 2;
            % end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, 'alpha')), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end
        end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % -------- Plot 2: Mean amplitudes across Power x Connectivity combos --------
        if strcmpi(plots_show, 'MEAN_AMP')
            
             figure('Name','Relative dispersion across Power × Phase combos','Color','w');
            hold on;

            % x positions for the 4 phase bins only
            xpos = 1:2;

            % Assuming the order in G.mean_VARrel_combos is:
            % [Low-Trough, Low-Rising, Low-Peak, Low-Falling, ...
            %  High-Trough, High-Rising, High-Peak, High-Falling]

            low_mean  = G.mean_MEAN_combos(1:2);
            low_sem   = G.sem_MEAN_combos(1:2);

            high_mean = G.mean_MEAN_combos(3:4);
            high_sem  = G.sem_MEAN_combos(3:4);

            % choose different colors for low and high
            col_high = [253, 70, 101] / 255;    % Vibrant red
            col_low = [59, 82, 139] / 255;    % Medium blue

            % Low power
            h1 = errorbar(xpos, low_mean, low_sem, '-o', ...
                'Color', col_low, ...
                'MarkerFaceColor', col_low, ...
                'MarkerEdgeColor', col_low, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

            % High power
            h2 = errorbar(xpos, high_mean, high_sem, '-o', ...
                'Color', col_high, ...
                'MarkerFaceColor', col_high, ...
                'MarkerEdgeColor', col_high, ...
                'MarkerSize', 10, ...
                'LineWidth', 4, ...
                'CapSize', 25);

               % 100% reference line
            yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

            ylim([-0.2 0.2]);
            xlim([0.5 2.50]);

            set(gca, ...
                'XTick', xpos, ...
                'XTickLabel', {'Low','High'}, ...
                'FontSize', 12, ...
                'LineWidth', 1.4);

        
            %xlabel('Phase', 'FontSize', 14, 'FontWeight','bold');
            %ylabel(sprintf('Normalized within-state %s (%%)', var_method), ...
            %    'FontSize', 14, 'FontWeight','bold');

            legend([h2 h1], {'High mu power','Low mu power','Random'}, ...
                'Location', 'northeast', ...
                'Box', 'off', 'FontSize', 14);

     
            set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
            ax = gca;
            set(ax, 'FontSize', 14);
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.LineWidth = 2;
            fontsize(gcf, scale=1.5)
            fontname('Arial')
            
            % % Define x-positions
            % combo_step = 0.25;
            % group_gap  = 0.35;
            % start_x    = 0.25;
            % xpos = zeros(1, nTotalBins);
            % 
            % % Position combos (adjust indices based on your combo structure)
            % nCombosPerGroup = nCombos / 2; % Assuming two groups (Low/High)
            % xpos(1:nCombosPerGroup) = start_x + (0:(nCombosPerGroup-1)) * combo_step;
            % xpos((nCombosPerGroup+1):nCombos) = xpos(nCombosPerGroup) + group_gap + (1:nCombosPerGroup) * combo_step;
            % xpos(nTotalBins) = xpos(nCombos) + group_gap + combo_step; % Random position
            % 
            % for iB = 1:nBands
            %     figure('Name','Mean amplitude across Power × Connectivity combos','Color','w'); hold on;
            % 
            %     % Get data for this band
            %     yComb = G.mean_MEAN_combos(iB, 1:nCombos);
            %     eComb = G.sem_MEAN_combos(iB, 1:nCombos);
            %     yRand = G.mean_MEAN_rand(iB);
            %     eRand = G.sem_MEAN_rand(iB);
            % 
            %     % Mean ± SEM for combos
            %     errorbar(xpos(1:nCombos), yComb, eComb, ...
            %         'o', ...
            %         'Color', col_combos, ...
            %         'MarkerFaceColor', col_combos, ...
            %         'MarkerEdgeColor', col_combos, ...
            %         'MarkerSize', 10, ...
            %         'LineWidth', 4, ...
            %         'CapSize', 25);
            % 
            %     % Random
            %     errorbar(xpos(nTotalBins), yRand, eRand, ...
            %         'o', ...
            %         'Color', col_rand, ...
            %         'MarkerFaceColor', col_rand, ...
            %         'MarkerEdgeColor', col_rand, ...
            %         'MarkerSize', 10, ...
            %         'LineWidth', 4, ...
            %         'CapSize', 25);
            % 
            %     % Add connection line between combo means
            %     plot(xpos(1:nCombos), yComb, '-', ...
            %         'Color', col_combos, ...
            %         'LineWidth', 1.4);
            % 
            %     % Axis settings
            %     ylim([-0.5 0.5]);
            %     xlim([0 xpos(nTotalBins) + 0.25]);
            %     set(gca, ...
            %         'XTick', xpos, ...
            %         'XTickLabel', comboLabels, ...
            %         'XTickLabelRotation', 0, ...
            %         'FontSize', 10, ...
            %         'LineWidth', 1.2, ...
            %         'Box', 'off', ...
            %         'TickDir', 'out');
            %     set(gcf, 'Units', 'normalized', ...
            %         'Position', [0.515234375, 0.152083333333333, 0.389765625, 0.772916666666667]);
            % 
            %     ax = gca;
            %     ax.LineWidth = 1.4;
            %     % Add custom group labels (adjust based on your grouping)
            %     xticks_range = xticks;
            %     y_offset = ax.YLim(1) - 0.065 * diff(ax.YLim);
            %     text(mean(xticks_range(1:nCombosPerGroup)), y_offset, 'Low Power', ...
            %         'HorizontalAlignment', 'center', ...
            %         'Rotation', 0, ...
            %         'FontSize', 14);
            %     text(mean(xticks_range((nCombosPerGroup+1):nCombos)), y_offset, 'High Power', ...
            %         'HorizontalAlignment', 'center', ...
            %         'Rotation', 0, ...
            %         'FontSize', 14);
            % 
            %     % Axis labels
            %     xlabel('Power × Connectivity combination', 'FontSize', 14, 'FontWeight','bold');
            %     ylabel('MEP amplitude (log-transformed, z-scored)', 'FontSize', 14, 'FontWeight','bold');
            %     %title(sprintf('%s', bands(iB)), 'FontSize', 16);
            % 
            %     % Adjust label positions
            %     xlab = ax.XLabel;
            %     xlab.Position(2) = xlab.Position(2) - 2.5;
            %     ylab = ax.YLabel;
            %     ylab.Position(1) = ylab.Position(1) - 0.1;
            % 
            %     % Add vertical separators
            %     v1 = mean([xpos(nCombosPerGroup), xpos(nCombosPerGroup+1)]);   % Low ↔ High
            %     v2 = mean([xpos(nCombos), xpos(nTotalBins)]);                   % High ↔ Random
            %     plot([v1 v1], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            %     plot([v2 v2], ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % 
            %     % Horizontal reference line at random mean
            %     plot(xlim, [yRand yRand], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            % 
            %     fontsize(gcf, scale=1.3);
            %     fontname('Arial');
            % end

            % Save as high-quality PNG (300 dpi)
            if save_output
                if strcmpi(response, 'TEP')
                    print(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, 'alpha')), '-dpng', '-r300');
                elseif strcmpi(response, 'MEP')
                    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'MEAN_AMP', sprintf('MEAN_%s_%s_%s_%s_%s_%s_%s.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, 'alpha')), 'Resolution', 600);
                end
            end
        end
    end
end
