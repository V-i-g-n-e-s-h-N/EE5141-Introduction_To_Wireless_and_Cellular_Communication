clear;
clc;
close all;

rng(5141, "twister");
params = getSystemParameters();
theme = getPlotTheme();
[outputRoot, figureDir] = prepareReportOutput();
[preambleFd, preambleTd, preambleWithCp, activeSubcarriers] = buildRepeatedPreamble(params);

repetitionLag = params.N / 8;
maxUnambiguousCfoHz = 1 / (2 * repetitionLag * params.Ts);

figure("Name", "Question 1(a)", "Color", theme.figureColor);
tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");

ax1 = nexttile;
ax1.ColorOrder = theme.palette;
freqProfile = zeros(size(params.usefulSubcarriers));
activeMask = mod(params.usefulSubcarriers, 8) == 0;
freqProfile(activeMask) = abs(preambleFd(params.usefulCarrierBins(activeMask)));
stem(params.usefulSubcarriers, freqProfile, "filled", "LineWidth", 1.0, ...
    "Color", theme.palette(1, :), "MarkerFaceColor", theme.palette(1, :), ...
    "MarkerEdgeColor", theme.palette(1, :));
grid on;
xlim([-240 240]);
xlabelHandle = xlabel("Subcarrier label n");
ylabelHandle = ylabel("|D_p[n]|");
titleHandle = title("Frequency-domain preamble: every 8th useful carrier carries QPSK");
applyDarkAxes(ax1, theme);
styleAxisLabel(xlabelHandle, theme);
styleAxisLabel(ylabelHandle, theme);
styleTitle(titleHandle, theme);

ax2 = nexttile;
ax2.ColorOrder = theme.palette;
repeatedBlocks = reshape(preambleTd, repetitionLag, []);


plot(0:repetitionLag-1, real(repeatedBlocks), "LineWidth", 1.0);
grid on;
xlabelHandle = xlabel("Sample index within repeated block");
ylabelHandle = ylabel("Real{x_p[m]}");
titleHandle = title("Time-domain property: all eight 64-sample blocks overlap exactly");
applyDarkAxes(ax2, theme);
styleAxisLabel(xlabelHandle, theme);
styleAxisLabel(ylabelHandle, theme);
styleTitle(titleHandle, theme);
saveFigure(gcf, fullfile(figureDir, "Q1a_PreambleDesign.png"));

fprintf("Question 1(a)\n");
fprintf("Active preamble subcarriers = %d\n", numel(activeSubcarriers));
fprintf("Repeated block length = %d samples\n", repetitionLag);
fprintf("Maximum unambiguous CFO = +/- %.2f kHz\n\n", maxUnambiguousCfoHz / 1e3);

snrDbList = 0:2:20;
numTrials = 10;
mseVsSnr = zeros(size(snrDbList));

for snrIdx = 1:numel(snrDbList)
    snrDb = snrDbList(snrIdx);
    squaredErrors = zeros(numTrials, 1);

    for trialIdx = 1:numTrials
        channelImpulseResponse = drawChannel(params);
        cfoTrueHz = -params.maxOffsetHz + 2 * params.maxOffsetHz * rand;

        receivedClean = conv(preambleWithCp, channelImpulseResponse);
        noiseVariance = 10^(-snrDb / 10);
        noise = sqrt(noiseVariance / 2) * ...
            (randn(size(receivedClean)) + 1j * randn(size(receivedClean)));
        sampleIndex = (0:numel(receivedClean)-1).';
        received = exp(1j * 2 * pi * cfoTrueHz * sampleIndex * params.Ts) .* ...
            (receivedClean + noise);

        receivedWithoutCp = received(params.cpLen+1:params.cpLen+params.N);
        cfoEstimateHz = estimateScCfo(receivedWithoutCp, params, repetitionLag);
        squaredErrors(trialIdx) = abs(cfoTrueHz - cfoEstimateHz)^2;
    end

    mseVsSnr(snrIdx) = mean(squaredErrors);
end

figure("Name", "Question 1(b)", "Color", theme.figureColor);
ax3 = axes();
plot(ax3, snrDbList, 10 * log10(mseVsSnr + eps), "o-", ...
    "LineWidth", 1.2, "MarkerSize", 6, "Color", theme.palette(2, :), ...
    "MarkerFaceColor", theme.palette(2, :), "MarkerEdgeColor", theme.palette(2, :));
grid on;
xlabelHandle = xlabel("SNR (dB)");
ylabelHandle = ylabel("10 log_{10}(MSE of CFO estimate in Hz^2)");
titleHandle = title("SC CFO-estimation MSE for J = 10 trials");
applyDarkAxes(ax3, theme);
styleAxisLabel(xlabelHandle, theme);
styleAxisLabel(ylabelHandle, theme);
styleTitle(titleHandle, theme);
saveFigure(gcf, fullfile(figureDir, "Q1b_CfoMse.png"));

fprintf("Question 1(b)\n");
fprintf("The script draws CFO uniformly from [-28.65, 28.65] kHz in every trial.\n");
fprintf("The plotted MSE is in Hz^2 before conversion to dB.\n");
writeQuestion1Outputs(outputRoot, activeSubcarriers, repetitionLag, maxUnambiguousCfoHz, snrDbList, mseVsSnr);

function params = getSystemParameters()
params.N = 512;
params.subcarrierSpacing = 10e3;
params.bandwidth = 5.12e6;
params.Ts = 1 / params.bandwidth;
params.cpLen = round(6.25e-6 / params.Ts);
params.frameSymbols = 5;
params.symbolLen = params.N + params.cpLen;
params.maxOffsetHz = 28.65e3;
params.usefulSubcarriers = [-240:-1, 1:240];
params.usefulCarrierBins = mod(params.usefulSubcarriers, params.N) + 1;
params.pathDelays = [0 7 13 18 21 26];
pathPowersDb = [-3 0 -1 -4 -9 -17];
pathPowersLinear = 10.^(pathPowersDb / 10);
params.pathPowers = pathPowersLinear / sum(pathPowersLinear);
end

function [freqGrid, timeBlock, blockWithCp, activeSubcarriers] = buildRepeatedPreamble(params)
activeSubcarriers = params.usefulSubcarriers(mod(params.usefulSubcarriers, 8) == 0);
freqGrid = zeros(params.N, 1);
freqGrid(mod(activeSubcarriers, params.N) + 1) = randomQpsk(numel(activeSubcarriers));
timeBlock = ifft(freqGrid) * sqrt(params.N);
blockWithCp = [timeBlock(end-params.cpLen+1:end); timeBlock];
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

function cfoEstimateHz = estimateScCfo(receivedBlock, params, lagSamples)
correlationValue = sum(conj(receivedBlock(1:end-lagSamples)) .* receivedBlock(1+lagSamples:end));
cfoEstimateHz = angle(correlationValue) / (2 * pi * lagSamples * params.Ts);
end

function [outputRoot, figureDir] = prepareReportOutput()
outputRoot = fullfile(pwd, "report_outputs");
figureDir = fullfile(outputRoot, "figures");
if ~exist(outputRoot, "dir")
    mkdir(outputRoot);
end
if ~exist(figureDir, "dir")
    mkdir(figureDir);
end
end

function saveFigure(figHandle, filePath)
if exist("exportgraphics", "file") == 2
    exportgraphics(figHandle, filePath, "Resolution", 300);
else
    saveas(figHandle, filePath);
end
end

function writeQuestion1Outputs(outputRoot, activeSubcarriers, repetitionLag, maxUnambiguousCfoHz, snrDbList, mseVsSnr)
mseDb = 10 * log10(mseVsSnr(:) + eps);
csvData = [snrDbList(:), mseDb];
csvPath = fullfile(outputRoot, "Question1_MSE.csv");
if exist("writematrix", "file") == 2
    writematrix(csvData, csvPath);
else
    dlmwrite(csvPath, csvData, ",");
end

summaryPath = fullfile(outputRoot, "Question1_summary.txt");
fileId = fopen(summaryPath, "w");
fprintf(fileId, "Question 1 summary\n");
fprintf(fileId, "Active preamble subcarriers: %d\n", numel(activeSubcarriers));
fprintf(fileId, "Repetition lag in SC estimator: %d samples\n", repetitionLag);
fprintf(fileId, "Maximum unambiguous CFO: +/- %.2f kHz\n", maxUnambiguousCfoHz / 1e3);
fprintf(fileId, "Reported MSE is in Hz^2 before dB conversion.\n");
fprintf(fileId, "SNR_dB,CFO_MSE_Hz2_dB\n");
for idx = 1:numel(snrDbList)
    fprintf(fileId, "%g,%.6f\n", snrDbList(idx), mseDb(idx));
end
fclose(fileId);
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