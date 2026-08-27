%% Problem 1a(c): Modified Sum-of-Sinusoids — Zheng & Xiao (2002)
% Reference: IEEE Comm. Letters, Vol.6, No.6, June 2002, eqs. (3)-(4)

clear; clc; close all;

%% ---- Parameters ----
v  = 30;            % velocity (m/s)
fc = 2e9;           % carrier frequency (Hz)
c  = 3e8;           % speed of light (m/s)
W  = 50e3;          % one-sided bandwidth
Ts = 1/(2*W);       % sampling period: 10 us
fd = v*fc/c;        % max Doppler: 200 Hz
wd = 2*pi*fd;       % max angular Doppler (rad/s)
N  = 8000;          % number of samples
M  = 20;            % sinusoids per component (M>=8 sufficient per paper)

%% ---- Time vector ----
t = (0:N-1) * Ts;       % (1 x N) in seconds

%% ---- Random variables — one draw per channel realization ----
rng(42);
theta    = unifrnd(-pi, pi);           % scalar: randomizes Doppler angles
phi_n    = unifrnd(-pi, pi, 1, M);    % (1 x M): I component phases
varphi_n = unifrnd(-pi, pi, 1, M);    % (1 x M): Q component phases

%% ---- Doppler arrival angles: eq. (4) ----
n       = (1:M);                            % (1 x M)
alpha_n = (2*pi*n - pi + theta) / (4*M);   % (1 x M)  conditionally random

%% ---- Generate quadrature components: eqs. (3b), (3c) ----
% Vectorized: t is (N x 1), alpha/phases are (1 x M) -> outer product (N x M)
t_col = t(:);                                                % (N x 1)

Z_c = sqrt(2/M) * sum(cos(wd * t_col * cos(alpha_n) + phi_n),    2)';  % (1 x N)
Z_s = sqrt(2/M) * sum(cos(wd * t_col * sin(alpha_n) + varphi_n), 2)';  % (1 x N)

h = Z_c + 1j * Z_s;    % complex path gain, eq. (3a)

%% ---- Plot: Time-domain fading realization ----
t_ms = t * 1e3;    % convert to ms

figure('Name', 'Zheng-Xiao SOS Fading — v=30 m/s');

% --- Subplot 1: Magnitude in dB ---
subplot(4,1,1);
plot(t_ms, 20*log10(abs(h)), 'b');
xlabel('Time (ms)'); ylabel('Instantaneous power |{\it h}[n]|^2 (dB)', 'Interpreter', 'tex');
title_str = sprintf('Path Gain Envelope of the channel with M=%d sinusoids for {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with the modified SOS Method', M, v, fd);
title(title_str, 'Interpreter', 'tex');
grid on;

% --- Phase Plot: Zheng-Xiao SOS Method ---
subplot(4,1,2);
plot(t_ms, angle(h) * 180/pi, 'm');
xlabel('Time (ms)'); 
ylabel('Phase (deg)');
title_str_phase3 = sprintf('Instantaneous Phase of the channel with M=%d sinusoids for {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with the modified SOS Method', M, v, fd);
title(title_str_phase3, 'Interpreter', 'tex');
grid on;

% --- Subplot 2: In-phase Component ---
subplot(4,1,3);
plot(t_ms, Z_c, 'r');
xlabel('Time (ms)'); ylabel('Amplitude');
title_str2 = sprintf('Real part (I component) of the channel with M=%d sinusoids for {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with the modified SOS Method', M, v, fd); 
title(title_str2, 'Interpreter', 'tex');
grid on;

% --- Subplot 3: Quadrature Component ---
subplot(4,1,4);
plot(t_ms, Z_s, 'g');
xlabel('Time (ms)'); ylabel('Amplitude');
title_str3 = sprintf('Imaginary part (Q component) of the channel with M=%d sinusoids for {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with the modified SOS Method', M, v, fd); 
title(title_str3, 'Interpreter', 'tex');
grid on;

%% Figure 2: ACF and Cross-correlation
h_I = Z_c;
h_Q = Z_s;
Ns   = length(h);
Nlag     = min(40000, Ns-1);                               % lag range to display
[acf, lags] = xcorr(h,   Nlag, 'normalized');            % ACF of complex h
[xc,  ~   ] = xcorr(h_I, h_Q, Nlag, 'normalized');      % cross-corr I vs Q

tau_ms   = lags * Ts * 1e3;                              % lag axis in ms
R_theory = besselj(0, 2*pi * fd * abs(lags) * Ts);      % theoretical: J0(2*pi*fd*tau)

figure('Name', 'ACF and Cross-Correlation: Zheng-Xiao SOS');

% --- Subplot 1: Autocorrelation (ACF) ---
subplot(1,2,1);
plot(tau_ms, real(acf), 'b', 'LineWidth', 1.2); hold on;
plot(tau_ms, R_theory,  'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex'); 
ylabel('{\itR_h}[{\it\tau}]', 'Interpreter', 'tex');
title_str1 = sprintf('Autocorrelation (ACF) of the M=%d process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with the modified SOS Method', M, v, fd);
title(title_str1, 'Interpreter', 'tex');
legend('Simulated (Real Part)', 'Theoretical: {\itJ}_0(2\pi{\itf}_d{\it\tau})', 'Interpreter', 'tex');
grid on;

% --- Subplot 2: Cross-Correlation (Zheng-Xiao SOS) ---
xc_real  = real(xc);
rms_xc   = sqrt(mean(xc_real.^2));

subplot(1,2,2);
plot(tau_ms, xc_real, 'b', 'LineWidth', 1.2); hold on;
yline(0, 'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex'); 
ylabel('{\itR}_{IQ}[{\it\tau}]', 'Interpreter', 'tex');

title_str2 = sprintf('Cross-Correlation (I vs Q) of the M=%d process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with the modified SOS Method', M, v, fd);
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
theta2    = unifrnd(-pi, pi);
phi_n2    = unifrnd(-pi, pi, 1, M);
varphi_n2 = unifrnd(-pi, pi, 1, M);

alpha_n2 = (2*pi*n - pi + theta2) / (4*M);

Z_c2 = sqrt(2/M) * sum(cos(wd * t_col * cos(alpha_n2) + phi_n2),    2)';
Z_s2 = sqrt(2/M) * sum(cos(wd * t_col * sin(alpha_n2) + varphi_n2), 2)';

h2 = Z_c2 + 1j * Z_s2;

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
title_str_p = sprintf('Two Independent Path Gains with M=%d, {\\itv}=%d m/s, {\\itf}_d=%.0f Hz with the modified SOS Method', M, v, fd);
title(title_str_p, 'Interpreter', 'tex');
legend('Path 1', 'Path 2');
grid on;

subplot(1,2,2);
plot(tau_ms_p, xc_paths_real, 'w', 'LineWidth', 1.2);  hold on;
yline(0, 'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex');
ylabel('R_{h_1 h_2}[\tau]', 'Interpreter', 'tex');
title_str_xc = sprintf('Cross-Correlation of Path 1 vs Path 2, each with M=%d, {\\itv}=%d m/s with the modified SOS Method', M, v);
title(title_str_xc, 'Interpreter', 'tex');
legend('Simulated', 'Ideal (0)', 'Interpreter', 'tex');
grid on;
text(0.02, 0.95, sprintf('RMS = %.4f', rms_paths), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'Color', 'w', 'BackgroundColor', 'k', 'EdgeColor', 'k');