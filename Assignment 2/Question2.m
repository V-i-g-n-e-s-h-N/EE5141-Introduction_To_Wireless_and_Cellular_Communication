clear;
clc;
close all;

rng(5142, "twister");
params = getSystemParameters();
theme = getPlotTheme();
[outputRoot, figureDir] = prepareReportOutput();
[~, preambleTd] = buildRepeatedPreamble(params);

numFrames = 2;
snrDb = 6;
timingLag = params.N / 2;
trueWindowStarts0 = params.cpLen + (0:numFrames-1) * params.frameSymbols * params.symbolLen;

txStream = buildFrameStream(params, preambleTd, numFrames);

[metricSingle, decisionSingle, detectedSingle0] = simulateTimingCase(txStream, 1, snrDb, params, timingLag, numFrames);
channelImpulseResponse = drawChannel(params);
[metricMultipath, decisionMultipath, detectedMultipath0] = simulateTimingCase(txStream, channelImpulseResponse, snrDb, params, timingLag, numFrames);

figure("Name", "Question 2", "Color", theme.figureColor);
tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");

nexttile;
plotTimingMetric(metricSingle, decisionSingle, detectedSingle0, trueWindowStarts0, ...
    "Question 2(a): SC timing metric for a single-tap channel", theme);

nexttile;
plotTimingMetric(metricMultipath, decisionMultipath, detectedMultipath0, trueWindowStarts0, ...
    "Question 2(b): SC timing metric for the given multipath channel", theme);
saveFigure(gcf, fullfile(figureDir, "Q2_SC_TimingMetric.png"));

fprintf("Question 2\n");
fprintf("Single-tap detected FFT window starts (0-based samples): %s\n", mat2str(detectedSingle0));
fprintf("Multipath detected FFT window starts (0-based samples): %s\n", mat2str(detectedMultipath0));
fprintf("In multipath, the plateau broadens because delayed paths keep the correlation high over a wider region.\n");
fprintf("A safe FFT window is still chosen inside the CP-supported high-metric region around each preamble.\n");
writeQuestion2Outputs(outputRoot, trueWindowStarts0, detectedSingle0, detectedMultipath0);

function params = getSystemParameters()
params.N = 512;
params.bandwidth = 5.12e6;
params.Ts = 1 / params.bandwidth;
params.cpLen = round(6.25e-6 / params.Ts);
params.frameSymbols = 5;
params.symbolLen = params.N + params.cpLen;
params.usefulSubcarriers = [-240:-1, 1:240];
params.usefulCarrierBins = mod(params.usefulSubcarriers, params.N) + 1;
params.pathDelays = [0 7 13 18 21 26];
pathPowersDb = [-3 0 -1 -4 -9 -17];
pathPowersLinear = 10.^(pathPowersDb / 10);
params.pathPowers = pathPowersLinear / sum(pathPowersLinear);
end

function [freqGrid, timeBlock, blockWithCp] = buildRepeatedPreamble(params)
activeSubcarriers = params.usefulSubcarriers(mod(params.usefulSubcarriers, 8) == 0);
freqGrid = zeros(params.N, 1);
freqGrid(mod(activeSubcarriers, params.N) + 1) = randomQpsk(numel(activeSubcarriers));
timeBlock = ifft(freqGrid) * sqrt(params.N);
blockWithCp = [timeBlock(end-params.cpLen+1:end); timeBlock];
end

function txStream = buildFrameStream(params, preambleTd, numFrames)
preambleWithCp = [preambleTd(end-params.cpLen+1:end); preambleTd];
txStream = [];

for frameIdx = 1:numFrames
    txStream = [txStream; preambleWithCp];

    for symbolIdx = 1:params.frameSymbols-1
        dataGrid = zeros(params.N, 1);
        dataGrid(params.usefulCarrierBins) = randomQpsk(numel(params.usefulCarrierBins));
        timeBlock = ifft(dataGrid) * sqrt(params.N);
        txStream = [txStream; timeBlock(end-params.cpLen+1:end); timeBlock];
    end
end
end

function [metric, decisionMetric, detectedWindows0] = simulateTimingCase(txStream, channelImpulseResponse, snrDb, params, timingLag, numFrames)
receivedClean = conv(txStream, channelImpulseResponse(:));
noiseVariance = 10^(-snrDb / 10);
noise = sqrt(noiseVariance / 2) * ...
    (randn(size(receivedClean)) + 1j * randn(size(receivedClean)));
received = receivedClean + noise;

metric = schmidlCoxMetric(received, timingLag);
decisionMetric = conv(metric, ones(params.cpLen, 1) / params.cpLen, "same");
detectedWindows0 = detectFrameWindows(decisionMetric, params, numFrames);
end

function metric = schmidlCoxMetric(received, timingLag)
numPositions = numel(received) - 2 * timingLag + 1;
metric = zeros(numPositions, 1);

for startIdx = 1:numPositions
    firstHalf = received(startIdx:startIdx+timingLag-1);
    secondHalf = received(startIdx+timingLag:startIdx+2*timingLag-1);
    pValue = sum(conj(firstHalf) .* secondHalf);
    rValue = sum(abs(secondHalf).^2);
    metric(startIdx) = abs(pValue)^2 / (rValue^2 + eps);
end
end

function detectedWindows0 = detectFrameWindows(decisionMetric, params, numFrames)
detectedWindows0 = zeros(1, numFrames);
centerToWindowOffset0 = floor((params.cpLen - 1) / 2);

for frameIdx = 1:numFrames
    frameStart0 = (frameIdx - 1) * params.frameSymbols * params.symbolLen;
    searchStart = frameStart0 + 1;
    searchStop = min(numel(decisionMetric), frameStart0 + params.symbolLen + params.cpLen + params.pathDelays(end));
    [~, localPeakIdx] = max(decisionMetric(searchStart:searchStop));
    detectedWindows0(frameIdx) = searchStart + localPeakIdx - 2 + centerToWindowOffset0;
end
end

function plotTimingMetric(metric, decisionMetric, detectedWindows0, trueWindowStarts0, plotTitle, theme)
sampleAxis0 = 0:numel(metric)-1;
metric = metric / max(metric + eps);
decisionMetric = decisionMetric / max(decisionMetric + eps);

ax = gca;
ax.ColorOrder = theme.palette;
plot(sampleAxis0, metric, "LineWidth", 1.0, "Color", theme.palette(1, :), ...
    "DisplayName", "SC timing metric");
hold on;
plot(sampleAxis0, decisionMetric, "LineWidth", 1.2, "Color", theme.palette(2, :), ...
    "DisplayName", "CP-smoothed decision metric");
grid on;
xlabelHandle = xlabel("Sample index");
ylabelHandle = ylabel("Normalized metric");
titleHandle = title(plotTitle);
applyDarkAxes(ax, theme);
styleAxisLabel(xlabelHandle, theme);
styleAxisLabel(ylabelHandle, theme);
styleTitle(titleHandle, theme);

yl = ylim(ax);
for idx = 1:numel(trueWindowStarts0)
    if idx == 1
        line(ax, [trueWindowStarts0(idx) trueWindowStarts0(idx)], yl, ...
            "Color", theme.palette(3, :), "LineStyle", "--", "LineWidth", 1.0, ...
            "DisplayName", "True FFT window");
    else
        line(ax, [trueWindowStarts0(idx) trueWindowStarts0(idx)], yl, ...
            "Color", theme.palette(3, :), "LineStyle", "--", "LineWidth", 1.0, ...
            "HandleVisibility", "off");
    end
end

for idx = 1:numel(detectedWindows0)
    if idx == 1
        line(ax, [detectedWindows0(idx) detectedWindows0(idx)], yl, ...
            "Color", theme.palette(5, :), "LineStyle", ":", "LineWidth", 1.0, ...
            "DisplayName", "Derived FFT window");
    else
        line(ax, [detectedWindows0(idx) detectedWindows0(idx)], yl, ...
            "Color", theme.palette(5, :), "LineStyle", ":", "LineWidth", 1.0, ...
            "HandleVisibility", "off");
    end
end

legendHandle = legend("show", "Location", "best");
styleLegend(legendHandle, theme);
hold off;
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

function writeQuestion2Outputs(outputRoot, trueWindowStarts0, detectedSingle0, detectedMultipath0)
summaryPath = fullfile(outputRoot, "Question2_summary.txt");
fileId = fopen(summaryPath, "w");
fprintf(fileId, "Question 2 summary\n");
fprintf(fileId, "True FFT windows (0-based): %s\n", mat2str(trueWindowStarts0));
fprintf(fileId, "Detected FFT windows in single-tap case (0-based): %s\n", mat2str(detectedSingle0));
fprintf(fileId, "Detected FFT windows in multipath case (0-based): %s\n", mat2str(detectedMultipath0));
fprintf(fileId, "The reported derived timing corresponds to the CP-smoothed plateau center shifted by roughly N_CP/2 toward the FFT boundary.\n");
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

function styleLegend(legendHandle, theme)
set(legendHandle, "TextColor", theme.textColor, ...
    "Color", theme.axesColor, ...
    "EdgeColor", theme.gridColor);
end