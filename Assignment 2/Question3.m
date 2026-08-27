clear;
clc;
close all;

rng(5143, "twister");
params = getSystemParameters();
theme = getPlotTheme();

snrDbList = 0:3:30;
numTrials = 100;
pilotSymbols = randomQpsk(params.numPilots);

pilotBins0 = mod(params.pilotSubcarriers, params.N);
usefulBins0 = mod(params.usefulSubcarriers, params.N);

aPilotCp = exp(-1j * 2 * pi / params.N * (pilotBins0(:) * (0:params.cpLen-1)));
aUsefulCp = exp(-1j * 2 * pi / params.N * (usefulBins0(:) * (0:params.cpLen-1)));
aPilotKnown = exp(-1j * 2 * pi / params.N * (pilotBins0(:) * params.pathDelays));
aUsefulKnown = exp(-1j * 2 * pi / params.N * (usefulBins0(:) * params.pathDelays));

msePilotZf = zeros(size(snrDbList));
mseLinear = zeros(size(snrDbList));
mseFft = zeros(size(snrDbList));
mseFftWindowed = zeros(size(snrDbList));
mseMls = zeros(size(snrDbList));
mseDelayAware = zeros(size(snrDbList));
representativeSNR = 12;
representativeCaptured = false;

for snrIdx = 1:numel(snrDbList)
    snrDb = snrDbList(snrIdx);
    noiseVariance = 10^(-snrDb / 10);
    errorAccumulator = zeros(1, 6);

    for trialIdx = 1:numTrials
        channelImpulseResponse = drawChannel(params);
        trueCfrFull = fft(channelImpulseResponse, params.N);
        trueCfrPilots = trueCfrFull(params.pilotCarrierBins);
        trueCfrUseful = trueCfrFull(params.usefulCarrierBins);

        noise = sqrt(noiseVariance / 2) * ...
            (randn(params.numPilots, 1) + 1j * randn(params.numPilots, 1));
        receivedPilots = pilotSymbols .* trueCfrPilots + noise;
        zfPilots = receivedPilots ./ pilotSymbols;

        linearEstimate = interpolatePilots(zfPilots, params);
        fftEstimate = dftDenoise(zfPilots, params, false);
        fftWindowedEstimate = dftDenoise(zfPilots, params, true);
        mlsTapEstimate = aPilotCp \ zfPilots;
        mlsEstimate = aUsefulCp * mlsTapEstimate;
        delayAwareTapEstimate = aPilotKnown \ zfPilots;
        delayAwareEstimate = aUsefulKnown * delayAwareTapEstimate;

        if ~representativeCaptured && snrDb == representativeSNR && trialIdx == 1
            representativeCaptured = true;
            representativeTrueCfr = trueCfrUseful;
            representativeZfPilots = zfPilots;
            representativeLinear = linearEstimate;
            representativeFftWindowed = fftWindowedEstimate;
            representativeMls = mlsEstimate;
            representativeDelayAware = delayAwareEstimate;
        end

        errorAccumulator(1) = errorAccumulator(1) + mean(abs(trueCfrPilots - zfPilots).^2);
        errorAccumulator(2) = errorAccumulator(2) + mean(abs(trueCfrUseful - linearEstimate).^2);
        errorAccumulator(3) = errorAccumulator(3) + mean(abs(trueCfrUseful - fftEstimate).^2);
        errorAccumulator(4) = errorAccumulator(4) + mean(abs(trueCfrUseful - fftWindowedEstimate).^2);
        errorAccumulator(5) = errorAccumulator(5) + mean(abs(trueCfrUseful - mlsEstimate).^2);
        errorAccumulator(6) = errorAccumulator(6) + mean(abs(trueCfrUseful - delayAwareEstimate).^2);
    end

    msePilotZf(snrIdx) = errorAccumulator(1) / numTrials;
    mseLinear(snrIdx) = errorAccumulator(2) / numTrials;
    mseFft(snrIdx) = errorAccumulator(3) / numTrials;
    mseFftWindowed(snrIdx) = errorAccumulator(4) / numTrials;
    mseMls(snrIdx) = errorAccumulator(5) / numTrials;
    mseDelayAware(snrIdx) = errorAccumulator(6) / numTrials;
end

figure("Name", "Question 3", "Color", theme.figureColor);
ax1 = axes();
plot(ax1, snrDbList, 10 * log10(msePilotZf + eps), "o-", "LineWidth", 1.1, ...
    "MarkerSize", 5, "Color", theme.palette(1, :), ...
    "MarkerFaceColor", theme.palette(1, :), "MarkerEdgeColor", theme.palette(1, :)); hold on;
plot(ax1, snrDbList, 10 * log10(mseLinear + eps), "s-", "LineWidth", 1.1, ...
    "MarkerSize", 5, "Color", theme.palette(2, :), ...
    "MarkerFaceColor", theme.palette(2, :), "MarkerEdgeColor", theme.palette(2, :));
plot(ax1, snrDbList, 10 * log10(mseFft + eps), "d-", "LineWidth", 1.1, ...
    "MarkerSize", 5, "Color", theme.palette(3, :), ...
    "MarkerFaceColor", theme.palette(3, :), "MarkerEdgeColor", theme.palette(3, :));
plot(ax1, snrDbList, 10 * log10(mseFftWindowed + eps), "^-", "LineWidth", 1.1, ...
    "MarkerSize", 5, "Color", theme.palette(4, :), ...
    "MarkerFaceColor", theme.palette(4, :), "MarkerEdgeColor", theme.palette(4, :));
plot(ax1, snrDbList, 10 * log10(mseMls + eps), "v-", "LineWidth", 1.1, ...
    "MarkerSize", 5, "Color", theme.palette(5, :), ...
    "MarkerFaceColor", theme.palette(5, :), "MarkerEdgeColor", theme.palette(5, :));
plot(ax1, snrDbList, 10 * log10(mseDelayAware + eps), "*-", "LineWidth", 1.1, ...
    "MarkerSize", 6, "Color", theme.palette(6, :), ...
    "MarkerFaceColor", theme.palette(6, :), "MarkerEdgeColor", theme.palette(6, :));
grid on;
xlabelHandle = xlabel("SNR (dB)");
ylabelHandle = ylabel("10 log_{10}(MSE)");
titleHandle = title("Question 3: OFDM channel-estimation performance");
applyDarkAxes(ax1, theme);
styleAxisLabel(xlabelHandle, theme);
styleAxisLabel(ylabelHandle, theme);
styleTitle(titleHandle, theme);
legendHandle = legend({ ...
    "ZF on pilot tones only", ...
    "Linear interpolation", ...
    "FFT-based interpolation", ...
    "FFT-based interpolation with tapering", ...
    "Modified LS using tap support < N_{CP}", ...
    "Delay-aware LS with known tap locations"}, ...
    "Location", "southwest");
styleLegend(legendHandle, theme);
hold off;

figure("Name", "Question 3 Representative CFR", "Color", theme.figureColor);
ax2 = axes();
plot(ax2, params.usefulSubcarriers, 20 * log10(abs(representativeTrueCfr) + eps), "-", ...
    "LineWidth", 1.4, "Color", theme.palette(1, :)); hold on;
plot(params.pilotSubcarriers, 20 * log10(abs(representativeZfPilots) + eps), "o", ...
    "LineWidth", 1.0, "MarkerSize", 4, "Color", theme.palette(2, :), ...
    "MarkerFaceColor", theme.palette(2, :), "MarkerEdgeColor", theme.palette(2, :));
plot(params.usefulSubcarriers, 20 * log10(abs(representativeLinear) + eps), "--", ...
    "LineWidth", 1.0, "Color", theme.palette(3, :));
plot(params.usefulSubcarriers, 20 * log10(abs(representativeFftWindowed) + eps), "-.", ...
    "LineWidth", 1.1, "Color", theme.palette(4, :));
plot(params.usefulSubcarriers, 20 * log10(abs(representativeMls) + eps), ":", ...
    "LineWidth", 1.3, "Color", theme.palette(5, :));
plot(params.usefulSubcarriers, 20 * log10(abs(representativeDelayAware) + eps), "-", ...
    "LineWidth", 1.2, "Color", theme.palette(6, :));
grid on;
xlabelHandle = xlabel("Subcarrier label n");
ylabelHandle = ylabel("20 log_{10}|G[n]|");
titleHandle = title("Representative CFR estimate comparison at 12 dB");
applyDarkAxes(ax2, theme);
styleAxisLabel(xlabelHandle, theme);
styleAxisLabel(ylabelHandle, theme);
styleTitle(titleHandle, theme);
legendHandle = legend({ ...
    "True CFR", ...
    "ZF at pilot tones", ...
    "Linear interpolation", ...
    "FFT interpolation with tapering", ...
    "Modified LS", ...
    "Delay-aware LS"}, ...
    "Location", "best");
styleLegend(legendHandle, theme);
hold off;

fprintf("Question 3\n");
fprintf("Pilots are placed every 8 subcarriers at n = +/-8, +/-16, ..., +/-240.\n");
fprintf("The plotted MSE values are averaged per evaluated carrier and then averaged across J = 100 channel realizations.\n");
fprintf("Expected ordering at moderate-to-high SNR: delay-aware LS <= mLS <= FFT-based <= linear interpolation.\n");

function params = getSystemParameters()
params.N = 512;
params.bandwidth = 5.12e6;
params.Ts = 1 / params.bandwidth;
params.cpLen = round(6.25e-6 / params.Ts);
params.usefulSubcarriers = [-240:-1, 1:240];
params.usefulCarrierBins = mod(params.usefulSubcarriers, params.N) + 1;
params.pilotSubcarriers = [-240:8:-8, 8:8:240];
params.pilotCarrierBins = mod(params.pilotSubcarriers, params.N) + 1;
params.numPilots = numel(params.pilotSubcarriers);
params.pathDelays = [0 7 13 18 21 26];
pathPowersDb = [-3 0 -1 -4 -9 -17];
pathPowersLinear = 10.^(pathPowersDb / 10);
params.pathPowers = pathPowersLinear / sum(pathPowersLinear);
end

function linearEstimate = interpolatePilots(zfPilots, params)
realPart = interp1(params.pilotSubcarriers, real(zfPilots).', params.usefulSubcarriers, "linear");
imagPart = interp1(params.pilotSubcarriers, imag(zfPilots).', params.usefulSubcarriers, "linear");
linearEstimate = realPart(:) + 1j * imagPart(:);
end

function cfrEstimate = dftDenoise(zfPilots, params, useWindow)
pilotGrid = zeros(params.N, 1);
pilotGrid(params.pilotCarrierBins) = zfPilots;
impulseResponse = ifft(pilotGrid);
scaleFactor = 8;
filteredImpulseResponse = zeros(params.N, 1);
filteredImpulseResponse(1:params.cpLen) = impulseResponse(1:params.cpLen) * scaleFactor;

if useWindow
    taperLength = min(8, params.cpLen);
    window = ones(params.cpLen, 1);
    taperIndex = (0:taperLength-1).';
    window(end-taperLength+1:end) = 0.5 * (1 + cos(pi * (taperIndex + 1) / taperLength));
    filteredImpulseResponse(1:params.cpLen) = filteredImpulseResponse(1:params.cpLen) .* window;
end

fullCfrEstimate = fft(filteredImpulseResponse);
cfrEstimate = fullCfrEstimate(params.usefulCarrierBins);
end

function qpskSymbols = randomQpsk(symbolCount)
qpskSymbols = exp(1j * (pi / 4 + (pi / 2) * randi([0 3], symbolCount, 1)));
end

function channelImpulseResponse = drawChannel(params)
channelImpulseResponse = zeros(params.pathDelays(end) + 1, 1);
taps = sqrt(params.pathPowers(:) / 2) .* ...
    (randn(numel(params.pathDelays), 1) + 1j * randn(numel(params.pathDelays), 1));
channelImpulseResponse(params.pathDelays + 1) = taps;
end

function theme = getPlotTheme()
theme.figureColor = [0.00 0.00 0.00];
theme.axesColor = [0.05 0.05 0.05];
theme.textColor = [1.00 1.00 1.00];
theme.gridColor = [0.45 0.45 0.45];
theme.palette = [ ...
    0.15 0.85 1.00; ...
    1.00 0.72 0.20; ...
    0.35 1.00 0.45; ...
    1.00 0.40 0.75; ...
    1.00 0.32 0.32; ...
    0.60 0.65 1.00; ...
    0.20 1.00 0.85; ...
    1.00 0.92 0.30];
end

function applyDarkAxes(ax, theme)
set(ax, "Color", theme.axesColor, ...
    "XColor", theme.textColor, ...
    "YColor", theme.textColor, ...
    "GridColor", theme.gridColor, ...
    "MinorGridColor", theme.gridColor, ...
    "GridAlpha", 0.35, ...
    "MinorGridAlpha", 0.18, ...
    "LineWidth", 1.0);
box(ax, "on");
end

function styleTitle(titleHandle, theme)
set(titleHandle, "Color", theme.textColor, ...
    "BackgroundColor", [0.00 0.00 0.00], ...
    "EdgeColor", [0.00 0.00 0.00], ...
    "Margin", 6);
end

function styleAxisLabel(labelHandle, theme)
set(labelHandle, "Color", theme.textColor);
end

function styleLegend(legendHandle, theme)
set(legendHandle, "TextColor", theme.textColor, ...
    "Color", theme.axesColor, ...
    "EdgeColor", theme.gridColor);
end