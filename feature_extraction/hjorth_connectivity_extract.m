function [Tconn_mean] = hjorth_connectivity_extract(subID, baseDir, feature, varargin)
%% Simplified extractor for Hjorth connectivity tables
% Works on tables with columns:
% {'Trial_index','Condition','Freq_band','Signal_comb','iPLV','wPLI','Coh','imCoh'}
%
% Usage:
%   Tconn_mean = hjorth_connectivity_extract(subID, baseDir, feature)
%   Tconn_mean = hjorth_connectivity_extract(subID, baseDir, feature, 'Block', 2)

% ----------- USER SETTINGS -----------
onlyCondition = 'TMS';   % set [] to skip condition filtering
targetPair = 'C3-C4';    % connection of interest
targetMetric = 'iPLV';   % connectivity metric to extract
% -------------------------------------

% ---- Parse optional inputs ----
p = inputParser;
addParameter(p, 'Block', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
%addParameter(p, 'TargetPair', targetPair, @(x) ischar(x) || isstring(x));
%addParameter(p, 'TargetMetric', targetMetric, @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
block = p.Results.Block;
%targetPair = string(p.Results.TargetPair);
%targetMetric = string(p.Results.TargetMetric);

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
    if ~isfield(S, 'hjorth_connectivity_df')
        error('Variable "hjorth_connectivity_df" not found.');
    end
    X = S.hjorth_connectivity_df;

    % If cell array with headers in first row
    if iscell(X)
        headers = X(1,:);
        data = X(2:end,:);
        headers = matlab.lang.makeUniqueStrings(string(headers));
        T = cell2table(data, 'VariableNames', cellstr(headers));
    elseif istable(X)
        T = X;
    else
        error('Unsupported type for hjorth_connectivity_df: %s', class(X));
    end

    % Convert types
    if ismember('Trial_index', T.Properties.VariableNames)
        T.Trial_index = double(string(T.Trial_index));
    end
    for c = ["Condition","Freq_band","Signal_comb"]
        if ismember(c, T.Properties.VariableNames)
            T.(c) = string(T.(c));
        end
    end

    % Filter by condition
    if ~isempty(onlyCondition)
        T = T(strcmpi(T.Condition, onlyCondition), :);
    end

    % Filter by connection of interest
    if ~ismember(targetPair, T.Signal_comb)
        warning('Target pair "%s" not found for subject %s. Returning empty table.', targetPair, subID);
        Tconn_mean = table();
        return
    end
    T = T(strcmpi(T.Signal_comb, targetPair), :);

    % Check connectivity metric
    if ~ismember(targetMetric, T.Properties.VariableNames)
        error('Metric "%s" not found in table.', targetMetric);
    end

    % Keep relevant columns
    needCols = {'Trial_index','Freq_band',targetMetric};
    T = T(:, needCols);

    % Rename for clarity
    T.Properties.VariableNames{end} = 'Connectivity';

    % Compute mean connectivity per trial and frequency band
    grpVars = {'Trial_index','Freq_band'};
    Tconn_mean = groupsummary(T, grpVars, 'mean', 'Connectivity');
    Tconn_mean.Properties.VariableNames(strcmp(Tconn_mean.Properties.VariableNames,'mean_Connectivity')) = {'Hjorth_C3C4_Connectivity'};

    % Remove GroupCount
    if ismember('GroupCount', Tconn_mean.Properties.VariableNames)
        Tconn_mean = removevars(Tconn_mean, {'GroupCount'});
    end

    % Sort by frequency band
    if ismember('Freq_band', Tconn_mean.Properties.VariableNames)
        Tconn_mean = sortrows(Tconn_mean, 'Freq_band');
    end

    fprintf('OK: %s\n', files.name);
catch ME
    warning('FAILED: %s\n  -> %s', files.name, ME.message);
end

fprintf('Done.\n');
end
