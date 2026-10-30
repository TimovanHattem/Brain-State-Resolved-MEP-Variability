clear;          % Clear workspace variables
close all;      % Close all figure windows
clc;            % Clear command window

% Get full path of current script and move to its directory
filePath = matlab.desktop.editor.getActiveFilename;
mainpath = fileparts(filePath);   % Adjust to project root folder structure
cd(mainpath);

% Set paths to relevant toolboxes
TOOLBOXPATH = 'D:\2025-08 T2TVAR\Analysis\toolboxes';
addpath(fullfile(TOOLBOXPATH, 'eeglab2024.2'));
eeglab nogui

%% Define input data

bands_of_interest = {'alpha'};
splits = 'octiles';
nBins = 8;  % <-- number of power bins
prestim_feature = 'Power';

% % REFTEP
DATAIN = '';
DATAIN_SPLIT = '';
DATAOUT_VIS = '';
DATAIN_EP = '';
all_subjects = load('');
load(fullfile(DATAIN_EP, 'alldata_MEP_master.mat'), 'alldata_MEP_master'); % for synchronizing rejected trials

%% PSDs

% ---- Welch settings (fixed across all subjects)
% (uses EEG.srate from each subject; if it differs across subjects you should resample first)
winSec    = 0.5;
plotRange = [2 80];

% ---- containers for group-level PSDs: cell array of [nFreq x nSubjects] per bin per band
P_allBins_allBands = cell(numel(bands_of_interest), nBins);
f_group            = cell(numel(bands_of_interest), 1);

for fb = 1:numel(bands_of_interest)
    for bn = 1:nBins
        P_allBins_allBands{fb, bn} = [];
    end
end

% Optional: flat trial-level PSD storage per band x bin
flat_allBins_allBands = cell(numel(bands_of_interest), nBins);
f_flat                = cell(numel(bands_of_interest), 1);

for fb = 1:numel(bands_of_interest)
    for bn = 1:nBins
        flat_allBins_allBands{fb, bn} = [];
    end
end

% =========================
% SUBJECT LOOP
% =========================
for i = 1:29 %length(all_subjects)
    subID = all_subjects{i};
    subids = strrep(subID, 'REFTEP_', '');
    subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
    if contains(subID, 'screen5_')
        subID = ['SCREEN5_', subids];
        subids = num2str(str2double(subids) + 500);
    end
    sub_num = ['sub-', subids];


    % Load pre-stimulus EEG data
    EEG = pop_loadset('filename', [sub_num,'_EEG_pre_processed.set'], ...
        'filepath', fullfile(DATAIN, subID));

    % Remove rejected trials
    rejected_trials = find(alldata_MEP_master(strcmpi({alldata_MEP_master.subID}, subID)).total_trials_rejected);
    EEG = pop_select(EEG, 'notrial', rejected_trials);

    fs = EEG.srate;

    win      = max(round(winSec*fs), 32);
    noverlap = round(0.5*win);
    nfft     = 2^nextpow2(win);

    labels = {EEG.chanlocs.labels};
    labels = strtrim(labels);

    channels_C3 = {'C3', 'FC1', 'FC5', 'CP1', 'CP5'};
    weights     = [1, -1/4, -1/4, -1/4, -1/4];

    % Select splits table once per subject
    splitsTbl = load_splits_table(DATAIN_SPLIT, sub_num, 1, prestim_feature, splits);

    % =========================
    % BAND LOOP
    % =========================
    for fb = 1:numel(bands_of_interest)
        thisBand = bands_of_interest{fb};
        Tb = splitsTbl(strcmpi(splitsTbl.Freq_band, thisBand), :);

        % =========================
        % BAND LOOP
        % =========================
        for fb = 1:numel(bands_of_interest)

            thisBand = bands_of_interest{fb};
            Tb = splitsTbl(strcmpi(splitsTbl.Freq_band, thisBand), :);

            % =========================
            % BIN LOOP
            % =========================
            for bn = 1:nBins

                % Trial indices for this bin
                idx_bin = double(Tb.Trial_index(Tb.Index == bn));

                % Skip empty bins
                if isempty(idx_bin)
                    warning('No trials found for subject %s, band %s, bin %d', sub_num, thisBand, bn);
                    continue;
                end

                % Select trials
                EEG_bin = pop_select(EEG, 'trial', idx_bin);

                % Data: chan x time x trials
                data_bin = EEG_bin.data;

                % Hjorth spatial filter
                data_bin = spatialFilter_new(data_bin, labels, channels_C3, weights);

                % Squeeze to time x trials
                data_bin = squeeze(data_bin);

                % Make sure shape is [time x trials]
                if isvector(data_bin)
                    data_bin = data_bin(:);
                end

                % Welch PSD per trial
                [Pbin, f] = welch_trials_impl(data_bin, win, noverlap, nfft, fs);

                % Subject-level mean PSD
                Pbin_m = mean(Pbin, 2);   % freq x 1

                % Restrict frequency range
                mask = (f >= plotRange(1) & f <= plotRange(2));
                f2   = f(mask);
                Pbin_m = Pbin_m(mask);

                % Frequency grid consistency check per band
                if isempty(f_group{fb})
                    f_group{fb} = f2;
                else
                    if numel(f_group{fb}) ~= numel(f2) || any(abs(f_group{fb} - f2) > 1e-6)
                        error('Frequency grid mismatch across subjects in band %s. Consistent fs/resampling needed.', thisBand);
                    end
                end

                % Append subject mean PSD as new column
                P_allBins_allBands{fb, bn}(:, end+1) = Pbin_m;

                % Optional: store all single-trial PSDs (flattened across subjects)
                Pbin_flat = Pbin(mask, :);   % freq x trials
                flat_allBins_allBands{fb, bn} = [flat_allBins_allBands{fb, bn}, Pbin_flat];

                if isempty(f_flat{fb})
                    f_flat{fb} = f2;
                end
            end
        end
    end
end

%% Plotting PSD for all 8 bins

figure('Color','w','Name','PSD');
tiledlayout(numel(bands_of_interest), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
baseColor = [150, 125, 160] / 255;

baseColor = [150, 125, 160] / 255;

brightColor = baseColor + (1-baseColor)*0.6;   % lighter but still purple

cmap = zeros(nBins,3);
for bn = 1:nBins
    t = (bn-1)/(nBins-1);
    cmap(bn,:) = (1-t)*brightColor + t*baseColor;
end
for fb = 1:numel(bands_of_interest)
    nexttile;
    hold on;

    f = f_group{fb};

    for bn = 1:nBins
        Pbin = P_allBins_allBands{fb, bn};   % freq x subjects

        if isempty(Pbin)
            continue;
        end

        % group mean in linear power, then convert to dB
        bin_g = mean(Pbin, 2, 'omitnan');
        plot(f,10*log10(bin_g), ...
            'LineWidth',2,...
            'Color',cmap(bn,:));

    end

    % Formatting
    grid on;
    xlim(plotRange);
    ylim([-30 0]);
    %xlabel('Frequency (Hz)', 'FontWeight', 'bold');
    %ylabel('PSD (dB)', 'FontWeight', 'bold');
   
end
legend(arrayfun(@(x) sprintf('Q%d', x), 1:nBins, 'UniformOutput', false), ...
    'Location', 'best');
ylim([-30 5])
set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
ax = gca;
set(ax, 'FontSize', 12);
ax.Box = 'off';
ax.TickDir = 'out';
ax.LineWidth = 1.6;
fontsize(gcf, scale=1.4)
fontname('Arial')

%% Plotting 1/f corrected PSD

figure('Color','w','Name','PSD 1/f corrected');

band_ranges = {
    'alpha', 8.0, 13.0;
    'beta',  13.0, 30.0
    };

baseColor = [150, 125, 160] / 255;

brightColor = baseColor + (1-baseColor)*0.6;   % lighter but still purple

cmap = zeros(nBins,3);
for bn = 1:nBins
    t = (bn-1)/(nBins-1);
    cmap(bn,:) = (1-t)*brightColor + t*baseColor;
end

for fb = 1:numel(bands_of_interest)
    thisBand = bands_of_interest{fb};

    nexttile;
    hold on;

    F = f_group{fb};

    % Keep only 2–80 Hz
    mask_full = (F >= 2) & (F <= 80);
    Ffull = F(mask_full);

    % Fit on background bands only
    mask_1f = ((Ffull >= 2 & Ffull <= 8) | (Ffull >= 30 & Ffull <= 80));
    mask_1f = mask_1f & ~(Ffull >= 49 & Ffull <= 51);

    x = log10(Ffull(mask_1f));

    validBins = [];
    h = gobjects(nBins,1);

    for bn = 1:nBins
        Pbin = P_allBins_allBands{fb, bn};   % freq x subjects

        if isempty(Pbin)
            continue
        end

        validBins(end+1) = bn;

        % group mean in linear power
        spectrum_bin = mean(Pbin, 2, 'omitnan');
        Sfull_bin = spectrum_bin(mask_full);

        % log-log 1/f fit
        y_bin = log10(Sfull_bin(mask_1f));
        p_bin = polyfit(x, y_bin, 1);
        yfit_log_bin = polyval(p_bin, log10(Ffull));
        yfit_lin_bin = 10.^yfit_log_bin;

        % subtract aperiodic component
        new_spec_db_bin = 10*log10(Sfull_bin) - 10*log10(yfit_lin_bin);

        % plot bin
        h(bn) = plot(Ffull, new_spec_db_bin, ...
            'LineWidth', 4, ...
            'Color', cmap(bn,:), ...
            'DisplayName', sprintf('Q%d', bn));
    end

    % title(thisBand);
    % xlabel('Frequency (Hz)', 'FontWeight', 'bold');
    % ylabel('Flattened PSD (dB, aperiodic removed)', 'FontWeight', 'bold');

end

yline(0, '--', 'HandleVisibility','off', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'Layer', 'bottom');

freq_low  = band_ranges{1, 2};
freq_high = band_ranges{1, 3};
fill([freq_low freq_high freq_high freq_low], ...
    [-4 -4 12 12], ...
    [0.7 0.7 0.7], ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');


valid = h(h~=0);          % remove empty handles
legend(flip(valid), ...
       flip(arrayfun(@(x) sprintf('Q%d',x),1:nBins,'UniformOutput',false)), ...
       'Location','northeast', 'Box', 'off', 'FontSize', 14);
grid off;
xlim([2 45]);
ylim([-4 12]);

set(gcf, 'Units', 'normalized', 'Position', [0.5578125,0.152083333333333,0.3471875,0.772916666666667]);
ax = gca;
set(ax, 'FontSize', 14);
ax.Box = 'off';
ax.TickDir = 'out';
ax.LineWidth = 2;
fontsize(gcf, scale=1.5)
fontname('Arial')

% Save as high-quality PNG (300 dpi)
exportgraphics(gcf, fullfile(DATAOUT_VIS, sprintf('PSDs_1f_correced.png')), 'Resolution', 600);

%% -------- local function --------
function [P, f] = welch_trials_impl(X, win, noverlap, nfft, fs)
% X expected: time x trials
X = double(X);

% ensure 2D [time x trials]
X = squeeze(X);
if ndims(X) ~= 2
    error("Expected X to be [time x trials], got %s", mat2str(size(X)));
end

% remove DC per trial
X = X - mean(X,1);

nTrials = size(X,2);
for k = 1:nTrials
    [p, f] = pwelch(X(:,k), win, noverlap, nfft, fs);
    if k == 1
        P = zeros(numel(p), nTrials);
    end
    P(:,k) = p;
end
end