%% Problem 1a(a): Levinson AR-based Noise Coloring Filter — Fading Simulator
% Jake's PSD model: AR(P) filter driven by white Gaussian noise

clear; clc; close all;

set(groot, 'DefaultFigureColor', 'k');
set(groot, 'DefaultAxesColor', 'k');
set(groot, 'DefaultTextColor', 'w');    % 'w' is White
set(groot, 'DefaultAxesXColor', 'w');
set(groot, 'DefaultAxesYColor', 'w');
clf;

%% ---- Parameters ----
v   = 30;           % velocity (m/s)
fc  = 2e9;          % carrier frequency (Hz)
c   = 3e8;          % speed of light (m/s)
W   = 50e3;         % one-sided bandwidth (2W = 100 kHz)
Ts  = 1 / (2*W);    % sampling period: 10 us
fd  = v * fc / c;   % max Doppler frequency: 200 Hz
N   = 8000;         % number of samples
P   = 100;          % AR model order (tune between 50-200)

%% ---- Step 1: Jake's autocorrelation (per I or Q component) ----
% Full complex ACF:  R_h[k] = J0(2*pi*fd*k*Ts)   (unit total power)
% I (or Q) component ACF: R_comp[k] = 0.5 * J0(...)  (half the power each)
lags   = 0:P;
R_comp = 0.5 * besselj(0, 2*pi * fd * Ts * lags);
R_comp(1) = R_comp(1) + 1e-6; 

%% ---- Step 2: Levinson-Durbin — solve Yule-Walker equations ----
% Returns AR polynomial a (length P+1, a(1)=1) and prediction error sigma2
[a, sigma2] = levinson(R_comp, P);

%% ---- Step 3: White Gaussian noise (independent I and Q) ----
rng(42);   % fix seed for reproducibility
w_I = sqrt(sigma2) * randn(1, N);
w_Q = sqrt(sigma2) * randn(1, N);

%% ---- Step 4: Color the noise — filter through 1/A(z) ----
% filter(B, A, x): B=1 (no numerator zeros), A=a (denominator polynomial)
h_I = filter(1, a, w_I);
h_Q = filter(1, a, w_Q);

h = h_I + 1j * h_Q;    % complex path gain

%% ---- Plot: Time-domain fading realization ----
t_ms = (0:N-1) * Ts * 1e3;   % time in ms

figure('Name', 'Levinson AR Fading — v=30 m/s');

subplot(4,1,1);
plot(t_ms, 20*log10(abs(h)), 'b');
xlabel('Time (ms)'); ylabel('Instantaneous power |{\it h}[n]|^2 (dB)', 'Interpreter', 'tex');
title_str = sprintf('Path Gain Envelope for an AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method', P, v, fd);
title(title_str, 'Interpreter', 'tex');
grid on;

subplot(4,1,2);
plot(t_ms, angle(h) * 180/pi, 'm');
xlabel('Time (ms)'); 
ylabel('Phase (deg)');
title_str_phase1 = sprintf('Instantaneous Phase of the AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method', P, v, fd);
title(title_str_phase1, 'Interpreter', 'tex');
grid on;


subplot(4,1,3);
plot(t_ms, h_I, 'r');
xlabel('Time (ms)'); ylabel('Amplitude');
title_str2 = sprintf('Real part (I component) of the AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method ', P, v, fd'); 
title(title_str2, 'Interpreter', 'tex');
grid on;

subplot(4,1,4);
plot(t_ms, h_Q, 'g');
xlabel('Time (ms)'); ylabel('Amplitude');
title_str3 = sprintf('Imaginary part (Q component) of the AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method', P, v, fd'); 
title(title_str3, 'Interpreter', 'tex');
grid on;

%% Figure 2: ACF and Cross-correlation
Ns   = length(h);
Nlag     = min(40000, Ns-1);                               % lag range to display
[acf, lags] = xcorr(h,   Nlag, 'normalized');            % ACF of complex h
[xc,  ~   ] = xcorr(h_I, h_Q, Nlag, 'normalized');      % cross-corr I vs Q

tau_ms   = lags * Ts * 1e3;                              % lag axis in ms
R_theory = besselj(0, 2*pi * fd * abs(lags) * Ts);      % theoretical: J0(2*pi*fd*tau)

figure('Name', 'ACF and Cross-Correlation: Levinson Method');

% --- Subplot 1: Autocorrelation (ACF) ---
subplot(1,2,1);
plot(tau_ms, real(acf), 'b', 'LineWidth', 1.2); hold on;
plot(tau_ms, R_theory,  'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex'); 
ylabel('{\itR_h}[{\it\tau}]', 'Interpreter', 'tex');
title_str1 = sprintf('Autocorrelation (ACF) of the AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method', P, v, fd);
title(title_str1, 'Interpreter', 'tex');
legend('Simulated (Real Part)', 'Theoretical: {\itJ}_0(2\pi{\itf}_d{\it\tau})', 'Interpreter', 'tex');
grid on;

% --- Subplot 2: Cross-Correlation ---
xc_real  = real(xc);
rms_xc   = sqrt(mean(xc_real.^2));

subplot(1,2,2);
plot(tau_ms, xc_real, 'b', 'LineWidth', 1.2); hold on;
yline(0, 'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex'); 
ylabel('{\itR}_{IQ}[{\it\tau}]', 'Interpreter', 'tex');

title_str2 = sprintf('Cross-Correlation (I vs Q) of the AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method', P, v, fd);
title(title_str2, 'Interpreter', 'tex');
legend('Simulated', 'Ideal (0)', 'Interpreter', 'tex');
grid on;

% Add RMS Text Box
text(0.02, 0.95, sprintf('RMS = %.4f', rms_xc), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'Color', 'w', 'BackgroundColor', 'k', 'EdgeColor', 'k');


%% ---- Second Independent Path (different random draw) ----
% Remark 2 in Zheng & Xiao: using distinct (theta_k, phi_nk, varphi_nk)
% creates uncorrelated faders with identical statistical properties

rng(99);   % different seed -> independent realization
w_I2 = sqrt(sigma2) * randn(1, N);
w_Q2 = sqrt(sigma2) * randn(1, N);
h_I2 = filter(1, a, w_I2);
h_Q2 = filter(1, a, w_Q2);

h2 = h_I2 + 1j * h_Q2;    % complex path gain

%% ---- Cross-correlation: h (path 1) vs h2 (path 2) ----
[xc_paths, lags_p] = xcorr(h, h2, Nlag, 'normalized');
xc_paths_real = real(xc_paths);
rms_paths     = sqrt(mean(xc_paths_real.^2));
tau_ms_p      = lags_p * Ts * 1e3;

%% ---- Figure 3: Inter-Path Correlation Plot ----
figure('Name', 'Inter-Path Correlation: Path 1 vs Path 2');

subplot(1,2,1);
plot(t_ms, 20*log10(abs(h)),  'b', 'LineWidth', 1);  hold on;
plot(t_ms, 20*log10(abs(h2)), 'r', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('|h[n]| (dB)');
title_str_p = sprintf('Two Independent Path Gains of the AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method', P, v, fd);
title(title_str_p, 'Interpreter', 'tex');
legend('Path 1', 'Path 2');
grid on;

subplot(1,2,2);
plot(tau_ms_p, xc_paths_real, 'w', 'LineWidth', 1.2);  hold on;
yline(0, 'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex');
ylabel('R_{h_1 h_2}[\tau]', 'Interpreter', 'tex');
title_str_xc = sprintf('Cross-Correlation of Path 1 vs Path 2, each of the AR(%d) process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Levinson Method', P, v, fd);
title(title_str_xc, 'Interpreter', 'tex');
legend('Simulated', 'Ideal (0)', 'Interpreter', 'tex');
grid on;
text(0.02, 0.95, sprintf('RMS = %.4f', rms_paths), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'Color', 'w', 'BackgroundColor', 'k', 'EdgeColor', 'k');