clear;          % Clear workspace variables
close all;      % Close all figure windows
clc;            % Clear command window

% Get full path of current script and move to its directory
filePath = matlab.desktop.editor.getActiveFilename;
mainpath = fileparts(filePath);   % Adjust to project root folder structure
cd(mainpath);

%% Hjorth Features

% ----------- USER SETTINGS -----------
region = 'M1'; %DLPFC %M1
response = 'MEP'; %MEP %TEP
prestim_feature = 'Power'; %Phase %Connectivity %Power
TEP_feature = 'p2p'; %p2p %AUC %LMFP %GMFP
TEP_window = 'allearly'; %% allearly %N15P30 %P30N45 %N45P60
band_of_interest = {'theta','alpha','low_beta','high_beta','low_gamma','high_gamma'};
nPerState        = 45;
fs               = 1000;
rng(51);
save             = false;
generateplots    = true;

if strcmpi(region, 'M1')
    % % REFTEP
    DATAIN_SPLIT = '';
    DATAOUT_VIS = '';
    %SCREEN
    % DATAIN_SPLIT = '';
    % DATAOUT_VIS = '';
    if strcmpi(response, 'MEP')
        DATAIN_MEP = '';
        %DATAIN_MEP = '';
        
        % Load your MEP data if not already in workspace
        if ~exist('alldata_MEP_master','var')
            load(fullfile(DATAIN_MEP, 'alldata_MEP_master.mat'), 'alldata_MEP_master');
        end
        all_subjects = {alldata_MEP_master.subID};
    elseif strcmpi(response, 'TEP')
        DATAIN_TEP   = '';

        roi_channels = {'C3'};

        % Load your TEP data if not already in workspace
        if ~exist('alldata_sensor_TEP_master','var')
            load(fullfile(DATAIN_TEP, 'alldata_sensor_TEP_master.mat'), 'alldata_sensor_TEP_master');
        end
        load(fullfile(DATAIN_TEP, 'times_eeg_post.mat'))
        load(fullfile(DATAIN_TEP, 'chanlocs_eeg.mat'))
        timeVec = times_eeg_post;
        chanlocs = chanlocs_eeg;
        all_subjects = {alldata_sensor_TEP_master.subject};
    end
end

count = 0;
resultsTable = struct([]);
for subidx = 1:numel(all_subjects)
    subID = all_subjects{subidx};
    subids = strrep(subID, 'REFTEP_', '');
    %subids = strrep(strrep(subID, 'SCREEN3_', ''), 'SCREEN4_', '');
    sub_num = horzcat('sub-', subids);

    % --- Load split table (Power or Phase)
    splitsTbl = load_splits_table(DATAIN_SPLIT, sub_num, 1, prestim_feature);
    if isempty(splitsTbl)
        warning('No splits found for %s', subID);
        continue
    end

    % --------- Power/Connectivity Analysis (high/low) ---------
    if strcmpi(prestim_feature, 'Power') | strcmpi(prestim_feature, 'Connectivity')
        for fb = 1:numel(band_of_interest)
            thisBand = band_of_interest{fb};

            % --- Load splits table of band
            Tb = splitsTbl(strcmpi(splitsTbl.Freq_band, thisBand), :);
            if isempty(Tb)
                continue
            end

            hi_idx  = double(Tb.Trial_index(Tb.Index == 1));
            lo_idx  = double(Tb.Trial_index(Tb.Index == 3));

            if numel(hi_idx) < nPerState || numel(lo_idx) < nPerState
                disp(['Not enough trials for ' subID])
                continue
            end

            rejects = alldata_MEP_master(subidx).meps_rejected;

            reject_strat = Tb.Index(rejects);
            kept_strat = Tb.Index(~rejects);

            % Save to results
            count = count + 1;
            resultsTable(count).subject = subID;
            resultsTable(count).type    = horzcat('Hjorth_', prestim_feature);
            resultsTable(count).band    = thisBand;
            resultsTable(count).rejected_high       = sum(reject_strat == 1);
            resultsTable(count).rejected_low        = sum(reject_strat == 3);
            resultsTable(count).rejected_rand       = sum(reject_strat == 2);
            resultsTable(count).kept_high       = sum(kept_strat == 1);
            resultsTable(count).kept_low        = sum(kept_strat == 3);
            resultsTable(count).kept_rand       = sum(kept_strat == 2);
        end
    end
end

% Convert to tidy table
T_all = struct2table(resultsTable);

%%
% Total counts across all subjects/bands
total_rejected_high = sum(T_all.rejected_high);
total_rejected_low  = sum(T_all.rejected_low);
total_rejected_rand = sum(T_all.rejected_rand);
total_kept_high = sum(T_all.kept_high);
total_kept_low  = sum(T_all.kept_low);
total_kept_rand = sum(T_all.kept_rand);

fprintf('Total Rejected: High=%d, Low=%d, Random=%d\n', ...
    total_rejected_high, total_rejected_low, total_rejected_rand);
fprintf('Total Kept: High=%d, Low=%d, Random=%d\n', ...
    total_kept_high, total_kept_low, total_kept_rand);

%% -------- Group-level summary --------

% ---------- Power/Connectivity summary ----------
if strcmpi(prestim_feature, 'Power') | strcmpi(prestim_feature, 'Connectivity')

    % Variables to average across subjects
    vars2avg = {'rejected_high', 'rejected_low', 'rejected_rand', ...
            'kept_high', 'kept_low', 'kept_rand'};

    % Group by frequency band
    G = grpstats(T_all, "band", {@mean,@std,@numel}, 'DataVars', vars2avg);

    % Compute SEM columns
    for v = vars2avg
        var = v{1};
        meanName = ['mean_' var];
        stdName  = ['std_' var];
        nName    = ['numel_' var];
        semName  = [var '_sem'];

        if ismember(stdName, G.Properties.VariableNames) && ismember(nName, G.Properties.VariableNames)
            G.(semName) = G.(stdName) ./ sqrt(G.(nName));
        else
            warning('Missing columns for %s', var);
        end
    end

    % Keep only relevant columns
    keepCols = ["band"];
    for v = vars2avg
        var = v{1};
        meanCol = ['mean_' var];
        semCol  = [var '_sem'];
        if ismember(meanCol, G.Properties.VariableNames) && ismember(semCol, G.Properties.VariableNames)
            keepCols = [keepCols, meanCol, semCol];
        end
    end
    G = G(:, keepCols);
end

%% Plot
% Calculate SEM
vars = {'rejected_high', 'rejected_low', 'rejected_rand', ...
        'kept_high', 'kept_low', 'kept_rand'};

% Prepare data
rejected_mean = [G.mean_rejected_high, G.mean_rejected_low, G.mean_rejected_rand];
rejected_sem = [G.rejected_high_sem, G.rejected_low_sem, G.rejected_rand_sem];
kept_mean = [G.mean_kept_high, G.mean_kept_low, G.mean_kept_rand];
kept_sem = [G.kept_high_sem, G.kept_low_sem, G.kept_rand_sem];

% Plot with error bars
figure('Position', [100 100 1000 400]);

subplot(1,2,1);
b = bar(rejected_mean);
hold on;
ngroups = size(rejected_mean, 1);
nbars = size(rejected_mean, 2);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
errorbar(x', rejected_mean, rejected_sem, 'k.', 'LineWidth', 1);
hold off;
set(gca, 'XTickLabel', G.band);
xlabel('Frequency Band');
ylabel('Mean Number of Trials');
title('Rejected Trials by Strategy');
legend({'High', 'Low', 'Random'}, 'Location', 'best');
grid on;

subplot(1,2,2);
b = bar(kept_mean);
hold on;
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
errorbar(x', kept_mean, kept_sem, 'k.', 'LineWidth', 1);
hold off;
set(gca, 'XTickLabel', G.band);
xlabel('Frequency Band');
ylabel('Mean Number of Trials');
title('Kept Trials by Strategy');
legend({'High', 'Low', 'Random'}, 'Location', 'best');
grid on;