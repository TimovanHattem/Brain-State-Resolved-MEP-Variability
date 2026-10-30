function[new_spec, freqs] = compute_PSD(EEG)
%% Computes power density spectrum for subject
% Inputs:
%   - EEG, struct = EEGLAB structured data
% Outputs:
%   - new_spec(f), array = power of each frequency
%   - freqs(f), array = contains all the necessary frequencies to plot

addpath('D:\2025-08 T2TVAR\Analysis\inspo\CheckSNR');

% Get relevant channel data
C3_data = [];
C3_data.fs = EEG.properties.samplingRate;
C3_data.x = spatialFilterCh(EEG, {'C3', 'FC1', 'FC5', 'CP1', 'CP5'}, [1 -0.25 -0.25 -0.25 -0.25]);
segment_length = 4*EEG.properties.samplingRate;

% Spectrum
C3_spec_TEMP = proc_spectrum(C3_data, [.5 45], 'win', hann(segment_length), 'step', segment_length-10);
F = C3_spec_TEMP.t;
spectrum = C3_spec_TEMP.x;

% Get pink noise mask
mask_1f = (F>=0.5 & F <= 8) | (F >= 30 & F <= 45);
mask_full = F>=.5 & F<45;

% Fit 1/f to data and remove
warning('off')
ft = fittype(@(m,b,x)b*x.^m);   % power function
fit_1f = fit(F(mask_1f)', spectrum(mask_1f), ft);
warning('on')
yfit = feval(fit_1f, F);
new_spec = spectrum - yfit;       % subject 1/f noise
new_spec = new_spec(mask_full);
freqs = F(mask_full);

fig_spectrum = figure('Renderer', 'painters', 'Position', [100 100 800 500]);
hold on;
plot(fit_1f, F(mask_1f), spectrum(mask_1f))
plot(F, spectrum, 'DisplayName','raw spectrum')
plot(freqs,new_spec, 'DisplayName','without 1/f')
yline(5, ':', 'DisplayName', 'SNR threshold')
ylabel('Power [dB]')
xlabel('F [Hz]')
title('Spectrum of C3-Hjorth')

ax = gca;
set(ax, 'Color', 'none')
set(ax.XAxis, 'TickDirection', 'out')
set(ax.YAxis, 'TickDirection', 'out')
exportgraphics(fig_spectrum, sprintf('%s/%s-C3-Hjorth-spectrum.pdf', OUT_, subjName), ...
    'BackgroundColor', 'none', 'ContentType', 'vector')

end