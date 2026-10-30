%% T2TVAR BRAIN-STATE-RESOLVED VARIABILITY QUANTIFICATION BETWEEN SESSIONS
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

prestim_feature     = 'Power'; %Power %Connectivity %Phase %PowerPhase %PowerConnectivity
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

do_bootstrap        = false;
n_boot_iter         = 10;

plots_show          = 'VARrel'; %VARrel %VAR

generate_plots      = true;

save_output         = true;

rng(51)

%% ----------- SELECT DATA INPUT -----------

% % REFTEP
reftep = true; screen = false;
DATAIN_SPLIT = '';
DATAOUT_VIS = '';
DATAIN_EP = '';
all_subjects = load('');

% % REFTEP REPEAT
DATAIN_SPLIT = '';
DATAOUT_VIS = '';
DATAIN_EP = '';
all_subjects = load('');
load(fullfile(DATAIN_EP, 'alldata_MEP_master.mat'), 'alldata_MEP_master');
load(fullfile(DATAIN_EP, 'alldata_MEP_master_repeat.mat'), 'alldata_MEP_master_repeat');
all_subjects_ses2 = {alldata_MEP_master_repeat.subID};
all_subjects_ses1 = strrep(all_subjects_ses2, '_rep', '');
all_subjects = [all_subjects_ses1 all_subjects_ses2];

%% ----------- EXTRACT HJORTH FEATURES -----------

if strcmpi(response, 'MEP')
    alldata_MEP_master(~contains({alldata_MEP_master.subID}, all_subjects_ses1)) = [];
    alldata_MEP_master = [alldata_MEP_master alldata_MEP_master_repeat];
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

    alldata_MEP_master(~contains({alldata_MEP_master.subID}, all_subjects_ses1)) = [];
    alldata_MEP_master = [alldata_MEP_master alldata_MEP_master_repeat];
end

sub_range = 1:numel(all_subjects);
count = 0;
resultsTable = struct([]);
for subidx = sub_range
    subID = all_subjects{subidx};
    subids = strrep(subID, 'REFTEP_', '');
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
    end
end
T_all = struct2table(resultsTable);

% Save output table
if save_output
    if strcmpi(response, 'TEP')
        outfile = sprintf('T_VAR_allblocks_%s_%s_betweensessions_%s_%s_%s_%s_%s_%s.csv', upper(region), horzcat('Hjorth_', prestim_feature), upper(response), upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits);
        outpath = DATAIN_EP;
    elseif strcmpi(response, 'MEP')
        outfile = sprintf('T_VAR_allblocks_%s_%s_betweensessions_%s_%s_%s_%s_%s.csv', ...
            upper(region), horzcat('Hjorth_', prestim_feature), upper(response), ...
            upper(var_method), site, splits, strjoin(bands_of_interest, '_'));
        outpath = DATAIN_EP;
    end
    writetable(T_all, fullfile(outpath, outfile));
end

T_ses1 = T_all(ismember(T_all.subject, all_subjects_ses1), :);  % first session
T_ses2 = T_all(ismember(T_all.subject, all_subjects_ses2), :);   % repeat session

%% -------- Group-level summary --------

% ---------- Power/Connectivity summary ----------
if strcmpi(prestim_feature, 'Power') || strcmpi(prestim_feature, 'Connectivity')

    % Variables to average across subjects
    vars2avg = {'MEAN_all', 'MEAN_high', 'MEAN_low', 'MEAN_rand', ...
        'VAR_all', 'VAR_high', 'VAR_low', 'VAR_rand', ...
        'VARrel_high', 'VARrel_low', 'VARrel_rand', ...
        'ratio_high', 'ratio_low', 'logratio_high', 'logratio_low'};

    G_all  = makeGroupSummary(T_all,  vars2avg);
    G_ses1 = makeGroupSummary(T_ses1, vars2avg);
    G_ses2 = makeGroupSummary(T_ses2, vars2avg);
    G = G_all;
end

%% Calcualte CCC
metrics = {'low', 'high','rand'};   % what you want to compare
metrics = strcat(plots_show, '_', metrics);

CCC = nan(numel(bands_of_interest), numel(metrics));
CCC_CI  = cell(numel(bands_of_interest), numel(metrics));
CCC_SEM = nan(numel(bands_of_interest), numel(metrics));

for b = 1:numel(bands_of_interest)
    bandName = bands_of_interest(b);

    for m = 1:numel(metrics)
        metricName = metrics{m};

        x = [];  % session 1
        y = [];  % session 2

        for s = 1:(length(all_subjects)/2)
            idx1 = strcmp(T_ses1.subject, all_subjects_ses1{s}) & strcmp(T_ses1.band, bandName);
            idx2 = strcmp(T_ses2.subject, all_subjects_ses2{s}) & strcmp(T_ses2.band, bandName);

            if ~any(idx1) || ~any(idx2)
                continue;
            end

            x(end+1) = T_ses1.(metricName)(idx1);
            y(end+1) = T_ses2.(metricName)(idx2);
        end

        if numel(x) >= 3   % need at least a few points
            CCC(b,m) = computeCCC(x, y);
        else
            CCC(b,m) = NaN;
        end

        % Bootstrap confidence intervals for CCC
        nBoot = 1000;
        bootVals = zeros(nBoot,1);

        for bi = 1:nBoot
            idx = randi(numel(x), numel(x), 1);  % resample with replacement
            bootVals(bi) = computeCCC(x(idx), y(idx));
        end

        CI_low  = quantile(bootVals, 0.025);
        CI_high = quantile(bootVals, 0.975);

        CCC_CI{b,m} = [CI_low, CI_high];
        CCC_SEM(b,m) = std(bootVals);
    end
end

% Put into a nice table
CCC_table = array2table(CCC, ...
    'VariableNames', metrics, ...
    'RowNames', cellstr(bands_of_interest));

%% ----------- PLOTTING -----------

if generate_plots
    if strcmpi(prestim_feature, 'Power') | strcmpi(prestim_feature, 'Connectivity')

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
        [~, order] = ismember(bands, string(G.band));
        G = G(order,:);
        nBands = height(G);

        % Define plot colours
        col_high = [253, 70, 101] / 255;    % Vibrant red
        col_low = [59, 82, 139] / 255;    % Medium blue
        col_rand = [120, 120, 130] / 255;  % Cool gray with slight blue tint
        col_line = [0.6 0.6 0.6];          % light grey connecting line

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure('Name','MEP dispersion across bands (Session 1+2)','Color','w'); hold on;
        tiledlayout(1, 3, 'TileSpacing','loose','Padding','compact');

        % SUBPLOT 1: VARrel across bands (Session 1+2)
        nexttile; hold on;

        cond_sep = 0.35;   % offset between Session 1 and Session 2
        mk_ses1  = 'o';
        mk_ses2  = '^';

        % X positions for Low=1, High=2, Random=3
        x_centers = [1, 2, 3];
        x_s1 = x_centers - cond_sep/2;   % Session 1 slightly left
        x_s2 = x_centers + cond_sep/2;   % Session 2 slightly right

        % ---- Extract data for the single band ----
        data_low_s1  = T_ses1.VARrel_low;
        data_high_s1 = T_ses1.VARrel_high;
        data_low_s2  = T_ses2.VARrel_low;
        data_high_s2 = T_ses2.VARrel_high;

        data_rand_s1 = unique(T_ses1.VARrel_rand);
        data_rand_s2 = unique(T_ses2.VARrel_rand);

        % ---- Means & SEMs ----
        mu  = @(x) mean(x, 'omitnan');
        sem = @(x) std(x, 'omitnan') / sqrt(sum(~isnan(x)));

        vals_s1 = [mu(data_low_s1),  mu(data_high_s1),  mu(data_rand_s1)];
        sems_s1 = [sem(data_low_s1), sem(data_high_s1), sem(data_rand_s1)];
        vals_s2 = [mu(data_low_s2),  mu(data_high_s2),  mu(data_rand_s2)];
        sems_s2 = [sem(data_low_s2), sem(data_high_s2), sem(data_rand_s2)];

        cols = [col_low; col_high; col_rand];

        % ---- Plot Session 1 (filled markers) ----
        for c = 1:3
            errorbar(x_s1(c), vals_s1(c), sems_s1(c), mk_ses1, ...
                'Color', cols(c,:), 'MarkerFaceColor', cols(c,:), ...
                'MarkerEdgeColor', cols(c,:), ...
                'MarkerSize', 12, 'LineWidth', 4, 'CapSize', 25);
        end

        % ---- Plot Session 2 (open markers) ----
        for c = 1:3
            errorbar(x_s2(c), vals_s2(c), sems_s2(c), mk_ses2, ...
                'Color', cols(c,:), 'MarkerFaceColor', cols(c,:), ...
                'MarkerEdgeColor', cols(c,:), ...
                'MarkerSize', 12, 'LineWidth', 4, 'CapSize', 25);
        end

        % ---- Separating dashed lines between conditions ----
        yl = [85 115]; ylim(yl);
        % for sep = [1.5, 2.5]
        %     plot([sep sep], yl, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'HandleVisibility','off');
        % end

        % ---- Axes formatting ----
        set(gca, 'XTick', x_centers, 'XTickLabel', {'Low', 'High', 'Random'}, ...
            'FontSize', 14, 'LineWidth', 1);
        xlim([0.5, 3.5]);
        axis square; box off;
        ax = gca; ax.TickDir = 'out'; ax.LineWidth = 2;
        yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');
        %xlabel('Power bin', 'FontSize', 11, 'FontWeight', 'bold');
        %ylabel(sprintf('Normalized within-state %s', var_method), 'FontSize', 11, 'FontWeight', 'bold');

        % ---- Legend ----
        h1 = plot(nan,nan, mk_ses1, 'Marker', 'o', 'Color','k',  'MarkerFaceColor','k',  'MarkerSize',8,'LineWidth',1.5);
        h2 = plot(nan,nan, mk_ses1, 'Marker', '^','Color','k', 'MarkerFaceColor','k', 'MarkerSize',8,'LineWidth',1.5);
        legend([h1 h2], {'Session 1','Session 2'}, ...
            'Location','northwest', 'Box','off', 'FontSize', 14);

        % SUBPLOTS 2: Individual test-retest scatter plots
        nexttile; hold on;

        xlim([70 130])
        ylim([70 130])
        % identity line

        plot(ylim, ylim, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);        

        for s = 1:(length(all_subjects)/2)


            % ---- indexing ----
            idx1 = strcmp(T_ses1.subject, all_subjects_ses1{s}) & strcmp(T_ses1.band, bands_of_interest{b});
            idx2 = strcmp(T_ses2.subject, all_subjects_ses2{s}) & strcmp(T_ses2.band, bands_of_interest{b});

            if ~any(idx1) || ~any(idx2)
                continue;
            end

            % ---- values ----
            v1_high = T_ses1.VARrel_high(idx1);
            v2_high = T_ses2.VARrel_high(idx2);

            v1_low  = T_ses1.VARrel_low(idx1);
            v2_low  = T_ses2.VARrel_low(idx2);

            % v1_rand  = T_ses1.VARrel_rand(idx1);
            % v2_rand  = T_ses2.VARrel_rand(idx2);

            % ---- thin line connecting LOW → HIGH (for THIS SUBJECT) ----
            plot([v1_low, v1_high], [v2_low, v2_high], '-', ...
                'Color', col_line, 'LineWidth', 1, 'HandleVisibility','off');

        end

        for s = 1:(length(all_subjects)/2)


            % ---- indexing ----
            idx1 = strcmp(T_ses1.subject, all_subjects_ses1{s}) & strcmp(T_ses1.band, bands_of_interest{b});
            idx2 = strcmp(T_ses2.subject, all_subjects_ses2{s}) & strcmp(T_ses2.band, bands_of_interest{b});

            if ~any(idx1) || ~any(idx2)
                continue;
            end

            % ---- values ----
            v1_high = T_ses1.VARrel_high(idx1);
            v2_high = T_ses2.VARrel_high(idx2);

            v1_low  = T_ses1.VARrel_low(idx1);
            v2_low  = T_ses2.VARrel_low(idx2);

            % v1_rand  = T_ses1.VARrel_rand(idx1);
            % v2_rand  = T_ses2.VARrel_rand(idx2);   

            % ---- scatter points ----
            h1 = scatter(v1_low,  v2_low,  220, col_low,  's', 'filled');
            h2 = scatter(v1_high, v2_high, 220, col_high, 'diamond', 'filled');
            %h3 = scatter(v1_rand, v2_rand, 100, col_rand, '^', 'filled');
        end
        
        axis square;
        %xlabel('Session 1', 'FontSize', 16, 'FontWeight','bold');
        %ylabel('Session 2', 'FontSize', 16, 'FontWeight','bold');
        h_legend = zeros(2, 1);
        h_legend(1) = plot(nan, nan, 'diamond', ...
            'Color', col_high, ...
            'MarkerFaceColor', col_high, ...
            'MarkerSize', 10);
        h_legend(2) = plot(nan, nan, 's', ...
            'Color', col_low, ...
            'MarkerFaceColor', col_low, ...
            'MarkerSize', 10);
        legend(h_legend, {'High', 'Low'}, ...
            'Location', 'northwest', 'Box', 'off', 'FontSize', 14);


        set(gca, 'TickDir','out','Box','off','FontSize',14);

        box on;
        set(gcf, 'Units', 'normalized', 'Position', [0.515234375,0.582638888888889,0.389765625,0.342361111111111]);
        ax = gca;
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.LineWidth = 2;
    end

    %yline(100, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    % fontsize(gcf,scale=1.4)
    % fontname('Arial')

    % Subplot 3 : CCC
    nexttile; hold on;

    % X positions
    x_low  = 1;
    x_high = 2;
    x_rand = 3;

    % ---- Extract values ----
    % Assumes CCC and CCC_CI are indexed as (band, metric)
    % metric 1 = Low, 2 = High, 3 = Random
    mu_low  = CCC(1, 1);
    mu_high = CCC(1, 2);
    mu_rand = CCC(1, 3);

    se_low  = CCC_CI{1,1} % CCC_SEM(1, 1); % CCC_CI{1,1}
    se_high = CCC_CI{1,2} % CCC_SEM(1, 2);
    se_rand = CCC_CI{1,3} % CCC_SEM(1, 3);

    se_low  = [mu_low  - CCC_CI{1,1}(1),  CCC_CI{1,1}(2)  - mu_low];
    se_high = [mu_high - CCC_CI{1,2}(1),  CCC_CI{1,2}(2)  - mu_high];
    se_rand = [mu_rand - CCC_CI{1,3}(1),  CCC_CI{1,3}(2)  - mu_rand];

    % ---- Plot errorbars ----
    errorbar(x_low,  mu_low,  se_low(1),  se_low(2),  'o', ...
        'Color', col_low,  'MarkerFaceColor', col_low,  'MarkerEdgeColor', col_low, ...
        'MarkerSize', 12, 'LineWidth', 4, 'CapSize', 25);

    errorbar(x_high, mu_high, se_high(1),  se_high(2), 'o', ...
        'Color', col_high, 'MarkerFaceColor', col_high, 'MarkerEdgeColor', col_high, ...
        'MarkerSize', 12, 'LineWidth', 4, 'CapSize', 25);

    errorbar(x_rand, mu_rand, se_rand(1),  se_rand(2), 'o', ...
        'Color', col_rand, 'MarkerFaceColor', col_rand, 'MarkerEdgeColor', col_rand, ...
        'MarkerSize', 12, 'LineWidth', 4, 'CapSize', 25);

    % ---- Axes formatting ----
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Low', 'High', 'Random'}, ...
        'FontSize', 14, 'LineWidth', 2);
    xlim([0.5 3.5]);
    ylim([-1 1])
    yticks(linspace(-1,1,9))
    fontsize(gcf, scale=1.5);
    fontname('Arial');
    yline(0, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);

    axis square; box off;
    ax = gca; ax.TickDir = 'out'; ax.LineWidth = 2;

    %xlabel('Power bin', 'FontSize', 14, 'FontWeight', 'bold');
    %ylabel('CCC',       'FontSize', 16, 'FontWeight', 'bold');

    % ---- Final styling ----
    set(gcf, 'Units','normalized', ...
        'Position',[0.48984375, 0.215972222222222, 0.41515625, 0.709027777777778]);


    % xlab = ax.XLabel; xlab.Position(2) = xlab.Position(2) - 0.05;
    % ylab = ax.YLabel; ylab.Position(1) = ylab.Position(1) - 0.1;

    set(gcf, 'Units', 'normalized', 'Position', [0.108203125,0.202083333333333,0.85,0.474305555555556]);

    % Save as high-quality PNG (300 dpi)
    if save_output
        if strcmpi(response, 'TEP')
            print(gcf, fullfile(DATAOUT_VIS, sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_%s_%s_betweensessions_CCC.png', region, horzcat('Hjorth_', prestim_feature), response, upper(TEP_feature), upper(TEP_window), upper(var_method), site, splits, strjoin(bands_of_interest, ''))), '-dpng', '-r300');
        elseif strcmpi(response, 'MEP')
            exportgraphics(gcf, fullfile(DATAOUT_VIS, sprintf('VAR_%s_%s_%s_%s_%s_%s_%s_betweensessions_CCC.png', region, horzcat('Hjorth_', prestim_feature), response, upper(var_method), site, splits, strjoin(bands_of_interest, ''))), 'Resolution', 600);
        end
    end
end

%% Helpers

function G_out = makeGroupSummary(T_in, vars2avg)

G = grpstats(T_in, "band", {@mean,@std,@numel}, 'DataVars', vars2avg);

for v = vars2avg
    var = v{1};
    meanName = ['mean_' var];      % e.g. 'mean_mean_high'
    stdName  = ['std_'  var];      % e.g. 'std_mean_high'
    nName    = ['numel_' var];     % e.g. 'numel_mean_high'
    semName  = [var '_sem'];       % e.g. 'mean_high_sem'

    if ismember(stdName, G.Properties.VariableNames) && ...
            ismember(nName,  G.Properties.VariableNames)
        G.(semName) = G.(stdName) ./ sqrt(G.(nName));
    else
        warning('Missing columns for %s', var);
    end
end

keepCols = "band";
for v = vars2avg
    var = v{1};
    meanCol = ['mean_' var];   % e.g. 'mean_mean_high'
    semCol  = [var '_sem'];    % e.g. 'mean_high_sem'
    if ismember(meanCol, G.Properties.VariableNames) && ...
            ismember(semCol,  G.Properties.VariableNames)
        keepCols = [keepCols, meanCol, semCol];
    end
end

G_out = G(:, keepCols);

end