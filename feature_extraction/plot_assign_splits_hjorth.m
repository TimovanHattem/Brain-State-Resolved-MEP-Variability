function [fig, stats] = plot_assign_splits_hjorth(assign_tbl, varargin)
% PLOT_ASSIGN_SPLITS - Visualize split assignments (low/mid/high) per band.
%
% Multiple bands are plotted in a 2×grid layout:
%   - Odd rows = histograms
%   - Even rows = trial time courses
%   - Bands run left-to-right across columns

% -------------------- Parse inputs --------------------
p = inputParser;
addParameter(p, 'Band', "", @(x) isstring(x) || ischar(x) || iscell(x));
addParameter(p, 'Method', "", @(x) isstring(x) || ischar(x));
addParameter(p, 'Bins', 20, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'Title', "", @(x) isstring(x) || ischar(x));
addParameter(p, 'SavePath', "", @(x) isstring(x) || ischar(x));
addParameter(p, 'ShowSummary', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'BlockBorder', [], @(x) isnumeric(x) || isempty(x));
parse(p, varargin{:});

bandArg   = string(p.Results.Band);
methodArg = lower(string(p.Results.Method));
nBins     = p.Results.Bins;
titleArg  = string(p.Results.Title);
savePath  = string(p.Results.SavePath);
doPrint   = p.Results.ShowSummary;

% -------------------- Checks --------------------
needCols = {'Trial_index','Freq_band','Index','Label','Method'};
if ~all(ismember(needCols, assign_tbl.Properties.VariableNames))
    error('assign_tbl must contain: %s', strjoin(needCols, ', '));
end

T = assign_tbl;
if ismember('Hjorth_C3_Power', T.Properties.VariableNames)
    powerCol = 'Hjorth_C3_Power';
elseif ismember('Hjorth_C3C4_Connectivity', T.Properties.VariableNames)
    powerCol = 'Hjorth_C3C4_Connectivity';
else
    error('roi_mean_tbl must contain either "Hjorth_C3C4_Connectivity" or "Hjorth_C3_Power" column.');
end

T.Trial_index  = double(T.Trial_index);
T.Freq_band    = string(T.Freq_band);
T.Value  = double(T.(powerCol));
T.Value = 10 * log10(T.Value);
T.Label        = string(T.Label);
T.Method       = string(T.Method);

% Determine bands to plot
allBands = unique(T.Freq_band, 'stable');
if isempty(bandArg)
    bandList = allBands(1);
elseif strcmpi(bandArg, "all")
    bandList = allBands;
else
    bandList = string(bandArg);
end

nBands = numel(bandList);
stats  = repmat(struct(), nBands, 1);

% -------- Layout parameters --------
nCols = min(3, nBands);                     % up to 3 bands per row
nBlocks = ceil(nBands / nCols);             % number of row-pairs
nRows = nBlocks * 2;                        % two rows per block

fig = figure('Color','w','Name', sprintf('Splits'), 'units','normalized','Position',[0.13,0.11,0.775,0.815]);
tl = tiledlayout(2, 4, 'TileSpacing','compact','Padding','compact');

% bandInfo = { ...
%     'theta',      'theta',      4.0,   8.0;
%     'alpha',      'alpha',      8.0,  12.0;
%     'low_beta',   'low-beta',  13.0,  18.0;
%     'high_beta',  'high-beta', 19.0,  35.0;
%     'low_gamma',  'low-gamma', 36.0,  58.0;
%     'high_gamma', 'high-gamma',58.0, 100.0};
bandInfo = { ...
    'theta',      'theta',      4.0,   8.0;
    'alpha',      '\mu-alpha',      8.0,  12.0;
    'beta',   'beta',  13.0,  30.0;
    'gamma', 'gamma', 30.0, 80.0};

bandInfoTbl = cell2table(bandInfo, ...
    'VariableNames', {'key','pretty','fmin','fmax'});
% -------------------- Iterate bands --------------------
for bi = 1:nBands
    band = bandList(bi);
    Tb = T(T.Freq_band == band, :);
    if isempty(Tb)
        warning('No rows for band "%s". Skipping.', band);
        continue
    end

    % Determine method
    methodList = unique(lower(Tb.Method));
    if methodArg ~= ""
        method = methodArg;
    elseif numel(methodList) == 1
        method = methodList(1);
    else
        method = "fromIndex";
    end

    % Extract data
    x = Tb.Value;
    idx = Tb.Index;
    valid = ~(isnan(x) | isnan(idx));
    x = x(valid);
    idx = idx(valid);
    trl = Tb.Trial_index(valid);

    % Compute thresholds
    cut1 = NaN; cut2 = NaN;
    if method == "median"
        cut1 = prctile(x, 50);
    elseif method == "tertiles"
        qq = prctile(x, [33.333 66.667]); cut1 = qq(1); cut2 = qq(2);
    elseif method == "tails25"
        qq = prctile(x, [25 75]); cut1 = qq(1); cut2 = qq(2);
    elseif method == "fromIndex"
        xLow  = x(idx==1); xMid = x(idx==2); xHigh = x(idx==3);
        if ~isempty(xLow) && ~isempty(xMid),  cut1 = (max(xLow)+min(xMid))/2; end
        if ~isempty(xHigh)&& ~isempty(xMid),  cut2 = (min(xHigh)+max(xMid))/2; end
    end

    % Counts & means
    nTotal = numel(x);
    nLow   = sum(idx==3);
    nMid   = sum(idx==2);
    nHigh  = sum(idx==1);
    meanLow  = mean(x(idx==3), 'omitnan');
    meanHigh = mean(x(idx==1), 'omitnan');
    diffMean = meanHigh - meanLow;

    % Colors
    colLow  = [0.2 0.45 0.95];
    colHigh = [0.9 0.35 0.2];
    colMid  = [0.6 0.6 0.6];
    edgeCol = [0 0 0];

    % Histogram edges
    if all(isfinite(x))
        edges = linspace(min(x), max(x), nBins+1);
    else
        edges = nBins;
    end

    % Compute tile positions
    col = mod(bi-1, nCols) + 1;
    block = ceil(bi / nCols);
    rowHist = 1;
    rowTime = 2;

    % --------- (1) Histogram ---------
    ax1 = nexttile(tl); hold(ax1,'on');
    k1  = ax1.Layout.Tile;   
    histogram(ax1, x, edges, 'FaceColor',[0.85 0.85 0.85], 'EdgeColor', edgeCol, 'FaceAlpha',0.7);
    histogram(ax1, x(idx==3), edges, 'FaceColor', colLow,  'EdgeColor','none', 'FaceAlpha',0.8);
    histogram(ax1, x(idx==1), edges, 'FaceColor', colHigh, 'EdgeColor','none', 'FaceAlpha',0.8);
    if ~isnan(cut1), xline(ax1, cut1, '--', 'Color', colLow,  'LineWidth', 2); end
    if ~isnan(cut2), xline(ax1, cut2, '--', 'Color', colHigh, 'LineWidth', 2); end
    xlabel(ax1,'Power'); ylabel(ax1,'# trials');
    matchIdx = strcmpi(bandInfoTbl.key, band);
    methodplot = method;
    if strcmpi(method, 'tails25')
        methodplot = 'quartiles';
    end
    if any(matchIdx)
        prettyBand = bandInfoTbl.pretty{matchIdx};
        freqRange  = sprintf('%.0f–%.0f Hz', bandInfoTbl.fmin(matchIdx), bandInfoTbl.fmax(matchIdx));
        ttl = sprintf('%s (%s) | %s', prettyBand, freqRange, upper(char(methodplot)));
    else
        ttl = sprintf('%s | %s', upper(char(band)), upper(char(methodplot)));
    end
    title(ax1, ttl);
    grid(ax1,'on');

    % --------- (2) Trial time course ---------
    ax2 = nexttile(tl, k1+4); hold(ax2,'on');
    cvec = repmat(colMid, numel(x), 1);
    cvec(idx==3,:) = repmat(colLow,  sum(idx==3), 1);
    cvec(idx==1,:) = repmat(colHigh, sum(idx==1), 1);
    [trlSort, ord] = sort(trl);
    scatter(ax2, trlSort, x(ord), 16, cvec(ord,:), 'filled', 'MarkerFaceAlpha', 0.6);
    if ~isnan(cut1), yline(ax2, cut1, '--', 'Color', colLow,  'LineWidth',1.5); end
    if ~isnan(cut2), yline(ax2, cut2, '--', 'Color', colHigh, 'LineWidth',1.5); end
    xlabel(ax2,'Trial'); ylabel(ax2,'Power');
    title(ax2,'Trial Course');
    grid(ax2,'on');

    % --- Add block borders if provided ---
    if ~isempty(p.Results.BlockBorder)
        yl = ylim(ax2);
        for bb = p.Results.BlockBorder
            xline(ax2, bb, 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '--');
        end
        ylim(ax2, yl); % restore limits
    end

    % --------- Stats ---------
    stats(bi).band      = band;
    stats(bi).method    = method;
    stats(bi).nTotal    = nTotal;
    stats(bi).nLow      = nLow;
    stats(bi).nMid      = nMid;
    stats(bi).nHigh     = nHigh;
    stats(bi).cut1      = cut1;
    stats(bi).cut2      = cut2;
    stats(bi).meanLow   = meanLow;
    stats(bi).meanHigh  = meanHigh;
    stats(bi).diffMean  = diffMean;
end
end
