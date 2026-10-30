function [Troi_mean] = hjorth_power_extract(subID, baseDir, feature, varargin)
%% Simplified extractor for Hjorth power tables
% Works on tables with columns:
% {'Trial_index','Condition','Freq_band','Power_C3','Power_C4','Power_FCC3H','Power_FCC4H'}
%
% Usage:
%   Troi_mean = hjorth_power_extract(subID, baseDir, feature)
%   Troi_mean = hjorth_power_extract(subID, baseDir, feature, 'Block', 2)

% ----------- USER SETTINGS -----------
onlyCondition = 'TMS';   % set [] to skip condition filtering
targetChannel = 'Power_C3'; % <- fixed channel for analysis %'Power_C3_relative'
% -------------------------------------

% ---- Parse optional inputs ----
p = inputParser;
addParameter(p, 'Block', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
parse(p, varargin{:});
block = p.Results.Block;

% ---- File selection ----
if isempty(block)
    files = dir(fullfile(baseDir, [subID, '_sensor_', feature, '_df.mat']));
else
    files = dir(fullfile(baseDir, [subID, '_Block', num2str(block), '_sensor_', feature, '_df.mat']));
end
fprintf('Found %d files.\n', numel(files));

inPath = fullfile(files.folder, files.name);
try
    S = load(inPath);
    if ~isfield(S, 'hjorth_power_table')
        error('Variable "hjorth_power_table" not found.');
    end
    X = S.hjorth_power_table;

    % If cell array with headers in first row
    if iscell(X)
        headers = X(1,:);
        data = X(2:end,:);
        headers = matlab.lang.makeUniqueStrings(string(headers));
        T = cell2table(data, 'VariableNames', cellstr(headers));
    elseif istable(X)
        T = X;
    else
        error('Unsupported type for hjorth_power_table: %s', class(X));
    end

    % Convert types
    if ismember('Trial_index', T.Properties.VariableNames)
        T.Trial_index = double(string(T.Trial_index));
    end
    for c = ["Condition","Freq_band"]
        if ismember(c, T.Properties.VariableNames)
            T.(c) = string(T.(c));
        end
    end

    % Filter condition
    if ~isempty(onlyCondition)
        T = T(strcmpi(T.Condition, onlyCondition), :);
    end

    % Check channel
    if ~ismember(targetChannel, T.Properties.VariableNames)
        error('Channel %s not found in table.', targetChannel);
    end

    % Keep relevant columns
    needCols = {'Trial_index','Freq_band',targetChannel};
    T = T(:, needCols);

    % Rename for clarity
    T.Properties.VariableNames{end} = 'Power';

    % Compute mean power per trial and frequency band
    grpVars = {'Trial_index','Freq_band'};
    Troi_mean = groupsummary(T, grpVars, 'mean', 'Power');
    Troi_mean.Properties.VariableNames(strcmp(Troi_mean.Properties.VariableNames,'mean_Power')) = {'Hjorth_C3_Power'};

    % Remove the GroupCount column if present
    if ismember('GroupCount', Troi_mean.Properties.VariableNames)
        Troi_mean = removevars(Troi_mean, {'GroupCount'});
    end

    fprintf('OK: %s\n', files.name);
catch ME
    warning('FAILED: %s\n  -> %s', files.name, ME.message);
end

fprintf('Done.\n');
end
