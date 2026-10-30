function splitsTbl = load_splits_table(baseDir, subID, block, prestim_feature, splits)
% Load splits table for Hjorth Power or Phase data

baseName = sprintf('%s_sensor_Hjorth_%s_df', subID, prestim_feature);

path = fullfile(baseDir, [baseName '_splits_', splits, '.mat']);

splitsTbl = [];

if exist(path, 'file')
    S = load(path);
    fn = fieldnames(S);
    splitsTbl = S.(fn{1});
    return
end

if isempty(splitsTbl)
    warning('No Power/Phase split file found for %s', subID);
end
end
