%% T2TVAR SPLIT TRIALS BASED ON PRE-STIMULUS FEATURES
% By Timo van Hattem (timo.van-hattem@uni-tuebingen.de)
% Brain Networks & Plasticity Lab, University of Tuebingen
%
% Description: Assign trials to brain-state bins based on pre-stimulation
% EEG features, including spectral power, interhemispheric connectivity,
% and oscillatory phase. Supports different binning approaches, including
% quantile-based splits and phase bins. Trial assignments can be saved
% and visualized to assess their distribution across recording blocks.

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
eeglab nogui

addpath('D:\2025-08 T2TVAR\Analysis\toolboxes\CircStat2012a');

%% ----------- USER SETTINGS -----------

region              = 'M1'; %M1 %DLPFC
mode                = 'hjorth'; %hjorth %allchan
bands_of_interest   = {'theta','alpha','beta', 'gamma'}; 
prestim_feature     = 'Connectivity'; %Power %Connectivity %Phase

splits              = 'tails25'; %tails25 %tertiles %median %quintiles %sextiles %octiles %deciles '4phasebins' '8phasebins'
nPhaseBins          = 4;   %4 %8
if strcmpi(prestim_feature, 'Phase') && nPhaseBins == 4
    splits          = '4phasebins';
elseif strcmpi(prestim_feature, 'Phase') && nPhaseBins == 8
    splits          = '8phasebins';
end

blockwise           = false;

save_output         = true;
generate_plots      = false;

%% Select input data

% % REFTEP
baseDir = '';
DATAIN_EP = '';
DATAOUT_VIS = '';
all_subjects = load('');
load(fullfile(DATAIN_EP, 'alldata_MEP_master.mat'), 'alldata_MEP_master');

% % REVISION BRAIN STIMULATION NO NON-RESPONSE REJECTION
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_without_nonresponse_rejection.mat'), 'alldata_MEP_master');

% % REVISION BRAIN STIMULATION WITH OUTLIER REJECTION
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_with_outlier_rejection.mat'), 'alldata_MEP_master');

% % REVISION BRAIN STIMULATION APB FDI separately
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_APB.mat'), 'alldata_MEP_master');
% load(fullfile(DATAIN_EP, 'alldata_MEP_master_FDI.mat'), 'alldata_MEP_master');

% % REFTEP REPEAT
% baseDir = '';
% DATAIN_EP = '';
% DATAOUT_VIS = '';
% all_subjects = load('');
% alldata_MEP_master = load(fullfile(DATAIN_EP, 'alldata_MEP_master_repeat.mat'), 'alldata_MEP_master_repeat').alldata_MEP_master_repeat;

%% HJORTH FEATURES

if (strcmpi(region, 'M1') && strcmpi(mode, 'hjorth'))

    for i = 1:size(all_subjects,2)
        subID = all_subjects{i};
        subids = strrep(subID, 'REFTEP_', '');
        subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
        if contains(subID, 'screen5_')
            subID = ['SCREEN5_', subids];
            subids = num2str(str2num(subids) + 500);
        end
        sub_num = horzcat('sub-', subids);

        fprintf('\nSplitting features: %s\n', subID)

        if strcmpi(prestim_feature, 'Power')
            feature_name = horzcat('Hjorth_', prestim_feature); %Power %IRASA
            T_roi = hjorth_power_extract(sub_num, baseDir, feature_name); % Hjorth C3

            % Remove rejected trials
            rejected_trials = find(alldata_MEP_master(i).total_trials_rejected);
            T_roi(ismember(T_roi.Trial_index, rejected_trials), :) = [];
            % Renumber trials
            [~, ~, new_idx] = unique(T_roi.Trial_index, 'stable');
            T_roi.Trial_index = new_idx;
            assign_tbl = feature_split(T_roi, 'Method', splits);

            if save_output
                output_filename = [sub_num, '_sensor_', feature_name, '_df_splits_', splits, '.mat'];
                output_path = fullfile(baseDir, output_filename);
                save(output_path, 'assign_tbl');
            end

        elseif strcmpi(prestim_feature, 'Connectivity')
            feature_name = horzcat('Hjorth_', prestim_feature); %Power %IRASA
            T_roi = hjorth_connectivity_extract(sub_num, baseDir, feature_name);

            % Remove rejected trials
            rejected_trials = find(alldata_MEP_master(i).total_trials_rejected);
            T_roi(ismember(T_roi.Trial_index, rejected_trials), :) = [];
            % Renumber trials
            [~, ~, new_idx] = unique(T_roi.Trial_index, 'stable');
            T_roi.Trial_index = new_idx;
            assign_tbl = feature_split(T_roi, 'Method', splits);

            % If NaN in connectivity computation, replace NaN with mean of
            % that band
            if sum(isnan(T_roi.Hjorth_C3C4_Connectivity))
                warning('%s: Found %.f NaN values. Taking mean for this trial.', subID, sum(isnan(T_roi.Hjorth_C3C4_Connectivity)));
                if sum(isnan(T_roi.Hjorth_C3C4_Connectivity)) > 10
                    return
                end         
                
                NaNs_idx = find(isnan(T_roi.Hjorth_C3C4_Connectivity));
                for i = 1:sum(isnan(T_roi.Hjorth_C3C4_Connectivity))
                    thisNaN_idx = NaNs_idx(i);
                    thisNaN_label = T_roi.Freq_band(thisNaN_idx);
                    thisLabel_mean = nanmean(T_roi.Hjorth_C3C4_Connectivity(T_roi.Freq_band == thisNaN_label));
                    T_roi.Hjorth_C3C4_Connectivity(thisNaN_idx) = thisLabel_mean;
                end
            end

            if save_output
                output_filename = [sub_num, '_sensor_', feature_name, '_df_splits_', splits, '.mat'];
                output_path = fullfile(baseDir, output_filename);
                save(output_path, 'assign_tbl');
            end
            
        elseif strcmpi(prestim_feature, 'Phase')
            feature_name = horzcat('Hjorth_', prestim_feature);
            T_roi = hjorth_phase_extract(sub_num, baseDir, feature_name);

            % Remove rejected trials
            rejected_trials = find(alldata_MEP_master(i).total_trials_rejected);
            T_roi(ismember(T_roi.Trial_index, rejected_trials), :) = [];
            % Renumber trials
            [~, ~, new_idx] = unique(T_roi.Trial_index, 'stable');
            T_roi.Trial_index = new_idx;

            % Extract phase values
            phase_values = double(T_roi.Phase_C3H);
            NaN_idx = find(isnan(phase_values));
            if sum(isnan(phase_values))
                warning('%s: Found %.f NaN phase values that got a random phase value', subID, sum(isnan(phase_values)));
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
            assign_tbl = table(T_roi.Trial_index, phase_values, phase_bin,  phase_labels,...
                'VariableNames', {'Trial_index','Phase_C3H','Phase_Bin', 'Phase_Label'});
            
            if save_output
                output_filename = [sub_num, '_sensor_', feature_name, '_df_splits_', splits, '.mat'];
                output_path = fullfile(baseDir, output_filename);
                save(output_path, 'assign_tbl');
            end

        end
    end
end

return

%% Plotting splits distribution

if generate_plots
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')

        for i = 1:size(all_subjects,2)
            subID = all_subjects{i};
            subids = strrep(subID, 'REFTEP_', '');
            subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
            if contains(subID, 'screen5_')
                subID = ['SCREEN5_', subids];
                subids = num2str(str2num(subids) + 500);
            end
            sub_num = horzcat('sub-', subids);

            input_path = fullfile(baseDir, [sub_num, '_sensor_', horzcat('Hjorth_', prestim_feature), '_df_splits_', splits, '.mat']);
            load(input_path)

            % Get block borders
            EEG = pop_loadset('filename', [sub_num,'_EEG_pre_processed.set'], ...
                'filepath', fullfile('Y:\VANHATTEM\2025-08 T2TVAR\Data\REFTEP_M1\preprocessed', subID));
            block_identifier = [EEG.epoch.block_identifier];  % convert from cell to numeric array
            % Remove rejected trials
            load(fullfile(DATAIN_EP, 'alldata_MEP_master.mat'), 'alldata_MEP_master');
            rejected_trials = find(alldata_MEP_master(strcmpi({alldata_MEP_master.subID}, subID)).total_trials_rejected);
            block_identifier(rejected_trials) = [];
            block_changes = find(diff(block_identifier) ~= 0);
            block_borders = block_changes + 0.5;  % between trials

            % Plot splits
            [fig, stats] = plot_assign_splits_hjorth(assign_tbl, 'Band', bands_of_interest, 'Method', splits,  'BlockBorder', block_borders);
            
            if save_output
                filename = fullfile(DATAOUT_VIS, 'Splits_distribution', [subID, '_', prestim_feature, '_split_distribution_', splits, '.png']);
                exportgraphics(fig, filename);
            end

            pause(0.5)
            close all
        end

    elseif strcmpi(prestim_feature, 'Phase')

        nonuniform_phase_distrubutions = {};
        NaN_in_table = {};
        for i = 1:size(all_subjects,2)
            subID = all_subjects{i};
            subids = strrep(subID, 'REFTEP_', '');
            subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
            if contains(subID, 'screen5_')
                subID = ['SCREEN5_', subids];
                subids = num2str(str2num(subids) + 500);
            end
            sub_num = horzcat('sub-', subids);

            files = dir(fullfile(baseDir, [sub_num, '_sensor_', 'Hjorth_Phase', '_df_splits_', splits, '.mat']));
            inPath = fullfile(files.folder, files.name);
            phase_splits = load(inPath).assign_tbl;

            if sum(isnan(phase_splits.Phase_C3H))
                fprintf('%d NaN in table: %s\n', sum(isnan(phase_splits.Phase_C3H)), sub_num);
                NaN_in_table = [NaN_in_table subID];
                NaN_trial = find(isnan(phase_splits.Phase_C3H));
                phase_splits(NaN_trial,:) = [];
            end

            p_rtest = circ_rtest(phase_splits.Phase_C3H);
            if p_rtest < .05
                fprintf('Uneven phase distribution: %s\n', sub_num);
                nonuniform_phase_distrubutions = [nonuniform_phase_distrubutions subID];
            end

            phase_distrubtions_plot = figure('units','normalized','outerposition',[0.5 0 0.5 1]);
            polarhistogram(phase_splits.Phase_C3H, 30)
            title(sprintf('%s | Rayleigh p = %.3g', sub_num, p_rtest))

            if save_output
                filename = fullfile(DATAOUT_VIS, 'Splits_distribution', [subID, '_Phase_split_distribution_', splits, '.png']);
                exportgraphics(phase_distrubtions_plot, filename);
            end

            close all
        end
    end
end
