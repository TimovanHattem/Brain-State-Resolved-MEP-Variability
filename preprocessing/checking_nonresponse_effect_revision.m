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

baseDir = '';

if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')
    percentMatrix = nan(29, 8);  % preallocate
else
    percentMatrix = nan(29, 4);  % preallocate
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


    feature_name = horzcat('Hjorth_', prestim_feature); %Power %IRASA
    if strcmpi(prestim_feature, 'Power')
        T_roi = hjorth_power_extract(sub_num, baseDir, feature_name); % Hjorth C3
    elseif strcmpi(prestim_feature, 'Connectivity')
        T_roi = hjorth_connectivity_extract(sub_num, baseDir, feature_name);
    elseif strcmpi(prestim_feature, 'Phase')
        T_roi = hjorth_phase_extract(sub_num, baseDir, feature_name);
    end

    % Remove rejected trials ONLY NON RESPONSES
    rejected_trials_all = alldata_MEP_master(subidx).total_trials_rejected;
    rejected_trials_nonresponse = alldata_MEP_master(subidx).idx_meps_rejected;
    rejected_trials_coil = rejected_trials_all-rejected_trials_nonresponse;

    T_roi(ismember(T_roi.Trial_index, rejected_trials_coil), :) = [];
    % Renumber trials
    [~, ~, new_idx] = unique(T_roi.Trial_index, 'stable');
    T_roi.Trial_index = new_idx;

    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')
        assign_tbl = feature_split(T_roi, 'Method', splits);

        thisBand = 'alpha';
        % Load splits table of band
        Tb = assign_tbl(strcmpi(assign_tbl.Freq_band, thisBand), :);

        rejected_bins = Tb(rejected_trials_nonresponse,:);

    elseif strcmpi(prestim_feature, 'Phase')
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


            rejected_bins = assign_tbl(rejected_trials_nonresponse,:);
            rejected_bins.Index = rejected_bins.Phase_Bin;
    end

    

    
    
    if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')
        allBins = 1:8;
    else
        allBins = 1:4;
    end
    counts = histcounts(rejected_bins.Index, [allBins, allBins(end)+1]);  % edges must be nBins+1 long
    percent = 100 * counts / sum(counts);

    if sum(counts) == 0
        percentMatrix(subidx, :) = counts;
    else
        percentMatrix(subidx, :) = percent;
    end

end


%%

meanPercent = mean(percentMatrix, 1, 'omitnan');
semPercent = std(percentMatrix, 0, 1, 'omitnan') / sqrt(size(percentMatrix, 1));

figure;
hold on;

% Bar plot of means
b = bar(allBins, meanPercent, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Error bars (SEM)
errorbar(allBins, meanPercent, semPercent, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 25);

% Overlay individual subject data points (jittered for visibility)
nSubjects = size(percentMatrix, 1);
jitter = 0.05;
for bIdx = 1:length(allBins)
    xJit = allBins(bIdx) + (rand(nSubjects,1) - 0.5) * jitter * 2;
    scatter(xJit, percentMatrix(:, bIdx), 25, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
end

hold off;

%xlabel('Bin', 'FontWeight', 'bold', 'FontSize', 14);
%ylabel('Percent of non-response trials (%)', 'FontWeight', 'bold', 'FontSize', 14);
%title(sprintf('%s', prestim_feature), 'FontSize', 14);
xticks(allBins);
xlim([0.5, 8.5]);
set(gca, 'FontSize', 12, 'Box', 'off');

if strcmpi(prestim_feature, 'Phase')
    xlim([0.5, 4.5]);
    xticklabels({'Trough', 'Rising', 'Peak', 'Falling'});
else
    xticklabels({'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8'});
end


set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.516666666666667,0.3471875,0.408333333333333]);
ax = gca;
set(ax, 'FontSize', 14);
ax.Box = 'off';
ax.TickDir = 'out';
ax.LineWidth = 2;
fontsize(gcf, scale=1.5)
fontname('Arial')

if save_output
    exportgraphics(gcf, fullfile(DATAOUT_VIS, 'VAR', sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_nonresponses.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
end
