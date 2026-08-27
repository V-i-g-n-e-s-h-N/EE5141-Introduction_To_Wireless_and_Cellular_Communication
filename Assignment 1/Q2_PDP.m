%% Problem 2: Frequency-Selective Fading — Power Delay Profile
% 2W = 10 MHz wideband signal, 6-path PDP, 2048-point FFT, 3 realizations

clear; clc; close all;

%% ---- System Parameters ----
W_bw  = 10e6;           % pass-band bandwidth 2W = 10 MHz
Ts    = 1 / W_bw;       % sampling period = 1/(2W) = 0.1 us
NFFT  = 2048;           % FFT size as specified
fs    = W_bw;           % sampling frequency = 10 MHz

%% ---- PDP Definition ----
gain_dB  = [-2, 0, -1, -6, -9, -14];       % path gains in dB
delay_us = [0, 1.8, 3.5, 5.7, 8.1, 12.3]; % path delays in microseconds
L = length(gain_dB);

%% ---- Step 1: Convert dB → linear variances ----
sigma2_lin = 10.^(gain_dB / 10);

%% ---- Step 2: Normalize so sum(sigma^2) = 1 ----
sigma2 = sigma2_lin / sum(sigma2_lin);

fprintf('\nNormalized variances:\n');
fprintf('  sigma^2 = '); fprintf('%.4f  ', sigma2); fprintf('\n');
fprintf('  Sum     = %.4f  (should be 1)\n\n', sum(sigma2));

%% ---- Step 3: Delay → tap indices ----
% Ts = 0.1 us  →  tap_k = round(delay_us / (Ts*1e6))
tap_idx = round(delay_us / (Ts * 1e6));   % convert us to samples

fprintf('Tap indices:\n');
for i = 1:L
    fprintf('  Path %d: delay=%.1f us → tap %d,  sigma^2=%.4f\n', ...
            i, delay_us(i), tap_idx(i), sigma2(i));
end
fprintf('\n');

%% ---- Step 4–6: Generate 3 realizations, compute H[k], plot ----
rng(42);
colors = {'b', 'r', 'g'};

figure('Name', 'Q2: Frequency Response | H(f) |^2 — 3 Realizations');
hold on;

for run = 1:3

    %% ---- Build sparse impulse response h[n] ----
    % h has length = max_tap + 1, but we zero-pad to NFFT anyway
    h = zeros(1, NFFT);      % pre-allocate at FFT size directly

    for i = 1:L
        % Each path gain a_i ~ CN(0, sigma^2_i)
        % Real and imag parts each ~ N(0, sigma^2_i / 2)
        std_each = sqrt(sigma2(i) / 2);
        a_i = std_each * (randn + 1j * randn);
        h(tap_idx(i) + 1) = a_i;   % +1 for 1-based MATLAB indexing
    end

    %% ---- 2048-point FFT → frequency response ----
    H = fft(h, NFFT);

    %% ---- Frequency axis ----
    f_MHz = (0:NFFT-1) * (fs/NFFT) / 1e6;   % 0 to fs in MHz

    %% ---- Plot 10*log10(|H|^2) ----
    plot(f_MHz, 10*log10(abs(H).^2), colors{run}, 'LineWidth', 1.2);
end

%% ---- Coherence bandwidth marker ----
% Approx Bc ~ 1/tau_max = 1/12.3e-6 ~ 81 kHz
tau_mean = sum(sigma2 .* delay_us);
tau_sq_mean = sum(sigma2 .* (delay_us.^2));
tau_rms = sqrt(tau_sq_mean - tau_mean^2);

Bc_kHz = 1 / (5 * tau_rms * 1e-6) / 1e3;
Bc_MHz = Bc_kHz / 1e3;   % convert to MHz to match f_MHz axis

% Draw vertical xlines at every multiple of Bc across the band
fhB_MHz = Bc_MHz : Bc_MHz : fs/1e6;   % positions in MHz (skip 0)
for k = 1:length(fhB_MHz)
    if k == 1
        xline(fhB_MHz(k), 'w--', 'LineWidth', 1.2, ...
              'DisplayName', sprintf('B_c \\approx %.0f kHz', Bc_kHz), ...
              'Interpreter', 'tex');
    else
        xline(fhB_MHz(k), 'w--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

xlabel('Frequency (MHz)', 'Interpreter', 'tex');
ylabel('10 log_{10} |{\itH}({\itf})|^2  (dB)', 'Interpreter', 'tex');
title({'Frequency Response of Fading Channel with Specified PDP using 3 Independent Channel Realizations'});
legend('Realization 1', 'Realization 2', 'Realization 3', sprintf('B_c \\approx %.0f kHz', Bc_kHz), 'Location', 'best', 'Interpreter', 'tex');
grid on;
xlim([0, fs/1e6]);