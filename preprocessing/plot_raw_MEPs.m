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
eeglab nogui;

%% Select input data
% % REFTEP
DATAIN_BASE  = '';
DATAOUT_BASE = '';
DATAOUT_VIS = '';
all_subjects = load('');

% % REFTEP REPEAT
% DATAIN_BASE  = '';
% DATAOUT_BASE = '';
% DATAOUT_VIS = '';
% all_subjects = load('');

%% Plot RAW individual MEPs

subID = 'REFTEP_018'; % INPUT

subids = strrep(subID, 'REFTEP_', '');
subids = strrep(strrep(subids, 'SCREEN3_', ''), 'screen5_', '');
if contains(subID, 'screen5_')
    subID = ['SCREEN5_', subids];
    subids = num2str(str2num(subids) + 500);
end
sub_num = horzcat('sub-', subids);

set_file = fullfile(DATAIN_BASE, subID, horzcat(sub_num, '_EMG_processed.set'));

% Load EMG dataset
EMG   = pop_loadset(set_file);
data  = EMG.data;
times = EMG.times;

% MEP window
mep_win = [10 50];
[~, idx_start] = min(abs(times - mep_win(1)));
[~, idx_end]   = min(abs(times - mep_win(2)));

x = times(idx_start:idx_end);
mep_data = data(:, idx_start:idx_end, :); % chans x time x trials
nTrials  = size(mep_data, 3);

% Figure
fig = figure('Name','MEP browser');
tiledlayout(1,2,'Padding','compact');

ax1 = nexttile; grid on;
ax2 = nexttile; grid on;

trial_idx = 1;

% Initial plot
plot(ax1, x, squeeze(mep_data(1,:,trial_idx)), 'k', 'LineWidth', 1.5);
ylabel(ax1,'EMG respones [\muV]'); xlabel(ax1,'Time [ms]');
title(ax1, sprintf('%s | %s | Trial %d / %d', ...
    sub_num, EMG.chanlocs(1).labels, trial_idx, nTrials));

plot(ax2, x, squeeze(mep_data(2,:,trial_idx)), 'k', 'LineWidth', 1.5);
ylabel(ax2,'EMG response [\muV]'); xlabel(ax2,'Time [ms]');
title(ax2, sprintf('%s | %s | Trial %d / %d', ...
    sub_num, EMG.chanlocs(2).labels, trial_idx, nTrials));
drawnow;

% Advance only with SPACE
while ishandle(fig) && trial_idx < nTrials
    waitforbuttonpress;
    if ~ishandle(fig); break; end

    if strcmpi(get(fig,'CurrentKey'),'space')
        trial_idx = trial_idx + 1;

        cla(ax1); cla(ax2);

        plot(ax1, x, squeeze(mep_data(1,:,trial_idx)), 'k', 'LineWidth', 1.5);
        ylabel(ax1,'EMG [\muV]'); xlabel(ax1,'Time [ms]');
        title(ax1, sprintf('%s | %s | Trial %d / %d', ...
            sub_num, EMG.chanlocs(1).labels, trial_idx, nTrials));

        plot(ax2, x, squeeze(mep_data(2,:,trial_idx)), 'k', 'LineWidth', 1.5);
        ylabel(ax2,'EMG [\muV]'); xlabel(ax2,'Time [ms]');
        title(ax2, sprintf('%s | %s | Trial %d / %d', ...
            sub_num, EMG.chanlocs(2).labels, trial_idx, nTrials));

        drawnow;
    end
end
