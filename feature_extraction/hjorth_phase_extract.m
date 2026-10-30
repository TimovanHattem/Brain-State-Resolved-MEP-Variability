function [T_roi] = hjorth_phase_extract(sub_num, baseDir, feature_name)

files = dir(fullfile(baseDir, [sub_num, '_sensor_', feature_name, '_df.mat']));
inPath = fullfile(files.folder, files.name);
S = load(inPath).hjorth_phase_df;  % <-- phase table
headers = S(1,:);
data    = S(2:end,:);
headers = matlab.lang.makeUniqueStrings(string(headers));
T = cell2table(data, 'VariableNames', cellstr(headers));

T_roi = T(:, {'Trial_index', 'Phase_C3H'});

end