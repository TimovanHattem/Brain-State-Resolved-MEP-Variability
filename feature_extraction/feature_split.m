function [assign_tbl, cuts_tbl] = feature_split(roi_mean_tbl, varargin)
% SPLIT_TRIALS - Assign trials to low/mid/high strata per frequency band
% using ROI-mean (long) power values.
%
% INPUT
%   roi_mean_tbl : table with columns:
%                  - Trial_index
%                  - Freq_band
%                  - (optional) Quant_status    ('Absolute' or 'Relative')
%                  - One numeric column:
%                       PSD_ROImean   (power)
%                       Power_mean    (power)
%                       Connectivity_mean (connectivity)
%
% Name-Value pairs:
%   'QuantStatus' : 'Absolute' or 'Relative' (optional)
%   'Method'      : 'median' (default) | 'tertiles' | 'tails25' | 'quintiles' | 'deciles'
%   'Bands'       : cellstr/string array of bands to include (default: all)
%
% OUTPUT
%   assign_tbl : table with columns:
%                - Trial_index
%                - Freq_band
%                - Quant_status
%                - PSD_ROImean/Power_mean
%                - Index   (1=high, 2=mid, 3=low; for 'median', mid only if == median)
%                - Label   ('low','mid','high')
%                - Method
%   cuts_tbl   : table per band with thresholds used:
%                - Freq_band
%                - Method
%                - Cut1    (median or lower cut)
%                - Cut2    (upper cut; NaN for 'median')
%
% Notes
% - Splitting is performed independently within each Freq_band.
% - For 'median', trials exactly equal to the median get Index=2 ('mid').
% - If a band has <3 valid trials, indices are set to NaN.
%
% Example
%   [assign_tbl, cuts_tbl] = split_trials(Troi_mean, ...
%       'QuantStatus','Absolute','Method','tails25');

% ---------- Parse inputs ----------
p = inputParser;
addParameter(p, 'QuantStatus', 'Absolute', @(x)ischar(x) || isstring(x));
addParameter(p, 'Method', 'median', @(x) any(strcmpi(x, {'median','tertiles','tails25','deciles','quintiles','octiles', 'sextiles'})));
addParameter(p, 'Bands', [], @(x) isempty(x) || iscellstr(x) || isstring(x));
parse(p, varargin{:});

quantSel = string(p.Results.QuantStatus);
method   = lower(string(p.Results.Method));
bandsSel = p.Results.Bands;

% ---------- Basic checks / standardize ----------
T = roi_mean_tbl;

% Accept either PSD_ROImean or Hjorth_C3_Power
if ismember('PSD_ROImean', T.Properties.VariableNames)
    powerCol = 'PSD_ROImean';
elseif ismember('Hjorth_C3_Power', T.Properties.VariableNames)
    powerCol = 'Hjorth_C3_Power';
elseif ismember('Hjorth_C3C4_Connectivity', T.Properties.VariableNames)
    powerCol = 'Hjorth_C3C4_Connectivity';
else
    error('roi_mean_tbl must contain either "PSD_ROImean" or "Hjorth_C3_Power" column.');
end

needCols = {'Trial_index','Freq_band', powerCol};
if ~all(ismember(needCols(1:2), T.Properties.VariableNames))
    error('roi_mean_tbl must contain columns: Trial_index, Freq_band, and a value column.');
end

% Normalize types
T.Trial_index  = double(T.Trial_index);
T.Freq_band    = string(T.Freq_band);
T.(powerCol)  = double(T.(powerCol));

% If Quant_status exists, normalize it
hasQuant = ismember('Quant_status', T.Properties.VariableNames);
if hasQuant
    T.Quant_status = string(T.Quant_status);
end
% Filter quantization type (Absolute/Relative)
if hasQuant && ~isempty(quantSel)
    T = T(strcmpi(T.Quant_status, quantSel), :);
end

% Optional band filter
if ~isempty(bandsSel)
    bandsSel = string(bandsSel);
    T = T(ismember(T.Freq_band, bandsSel), :);
end

if isempty(T)
    warning('No rows remaining after filtering Quant_status/Bands. Returning empty outputs.');
    assign_tbl = T;
    cuts_tbl   = table([],[],[],[], 'VariableNames', {'Freq_band','Method','Cut1','Cut2'});
    return
end

% ---------- Prepare outputs ----------
bands = unique(T.Freq_band, 'stable');
assign_parts = cell(numel(bands), 1);
cuts_parts   = cell(numel(bands), 1);

% ---------- Loop over bands ----------
for b = 1:numel(bands)
    band = bands(b);
    Tb = T(T.Freq_band == band, :);

    x = Tb.(powerCol);
    x_valid = ~isnan(x);
    n_valid = sum(x_valid);

    idx = nan(height(Tb),1);
    cut1 = NaN; cut2 = NaN;

    if n_valid >= 3
        xv = x(x_valid);

        switch method
            case "median"
                medv = prctile(xv, 50);
                cut1 = medv; % report as Cut1
                % Assign: < median -> 3 (low); == median -> 2 (mid); > median -> 1 (high)
                idx(x_valid & (x <  medv)) = 3;
                idx(x_valid & (x >  medv)) = 1;
                idx(x_valid & (x == medv)) = 2;

            case "tertiles"
                qs = prctile(xv, [33.333, 66.667]);
                cut1 = qs(1); cut2 = qs(2);
                % <= q1 -> 3 (low); q1< x <q2 -> 2 (mid); >= q2 -> 1 (high)
                idx(x_valid & (x <= qs(1))) = 3;
                idx(x_valid & (x >= qs(2))) = 1;
                idx(x_valid & (x >  qs(1) & x < qs(2))) = 2;

            case "tails25"
                qs = prctile(xv, [25, 75]);
                cut1 = qs(1); cut2 = qs(2);
                % <= 25th -> 3 (low); 25th<x<75th -> 2 (mid); >=75th -> 1 (high)
                idx(x_valid & (x <= qs(1))) = 3;
                idx(x_valid & (x >= qs(2))) = 1;
                idx(x_valid & (x >  qs(1) & x < qs(2))) = 2;

            case "quintiles"
                edges = prctile(xv, 20:20:80);   % 1x4
                cut1  = edges(1);                % 20th percentile
                cut2  = edges(end);              % 80th percentile

                % Initial assignment: 1 = lowest, 5 = highest
                idxAsc = discretize(x, [-Inf edges Inf]);
                idx = idxAsc;

            case "octiles"
                edges = prctile(xv, 12.5:12.5:87.5);   % 1x7 octile boundaries
                cut1  = edges(1);                      % 12.5th percentile
                cut2  = edges(end);                    % 87.5th percentile

                % Initial octile assignment: 1 = lowest, 8 = highest
                idxOct = discretize(x, [-Inf edges Inf]);
                idx = idxOct;

            case "deciles"
                edges = prctile(xv, 10:10:90);   % 1x9 decile boundaries
                cut1  = edges(1);                % 10th percentile
                cut2  = edges(end);              % 90th percentile

                % Initial decile assignment: 1 = lowest, 10 = highest
                idxDec = discretize(x, [-Inf edges Inf]);
                idx = idxDec;
        end
    else
        % Not enough data to split reliably
        idx(:) = NaN;
    end

    % Labels

    lbl = strings(size(idx));
    if method == "deciles"
        for k = 1:10
            lbl(idx == k) = sprintf("d%02d", k);
        end
    
    elseif method == "octiles"
    for k = 1:8
        lbl(idx == k) = sprintf("o%02d", k);
    end

    elseif method == "quintiles"
        for k = 1:5
            lbl(idx == k) = sprintf("q%02d", k); % q01=highest, q05=lowest
        end
    else
        lbl(idx == 1) = "high";
        lbl(idx == 2) = "mid";
        lbl(idx == 3) = "low";
    end

    % Assemble band results
    Tb.Index = idx;
    Tb.Label = lbl;
    Tb.Method = repmat(method, height(Tb), 1);
    assign_parts{b} = Tb;

    cuts_parts{b} = table(band, method, cut1, cut2, ...
        'VariableNames', {'Freq_band','Method','Cut1','Cut2'});
end

% Concatenate
assign_tbl = vertcat(assign_parts{:});
cuts_tbl   = vertcat(cuts_parts{:});

% Nice ordering
assign_tbl = movevars(assign_tbl, {'Index','Label','Method'}, 'After', (powerCol));

end
