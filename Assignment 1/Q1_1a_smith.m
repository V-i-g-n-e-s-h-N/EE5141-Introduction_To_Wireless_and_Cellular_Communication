%% Problem 1a(b): Smith's Model (FFT-based) Fading Simulator
% Reference: Rappaport's Wireless Communications, Ch.3

clear; clc; close all;

%% ---- Parameters ----
v  = 30;            % velocity (m/s)
fc = 2e9;           % carrier frequency (Hz)
c  = 3e8;           % speed of light (m/s)
W  = 50e3;          % one-sided bandwidth
Ts = 1/(2*W);       % sampling period: 10 us
fd = v*fc/c;        % max Doppler: 200 Hz
N  = 8192;          % FFT size as specified
fs = 1/Ts;          % sampling frequency
df = fs/N;          % frequency resolution

%% ---- Step 1: Frequency axis (centered at zero) ----
f = (-N/2 : N/2-1) * df;   % runs from -fs/2 to fs/2 - df

%% ---- Step 2: Jake's PSD ----
% S(f) = 1 / (pi*fd*sqrt(1-(f/fd)^2))  for |f| < fd, else 0
S = zeros(1, N);
idx = abs(f) < fd;
S(idx) = 1 ./ (pi * fd * sqrt(1 - (f(idx)/fd).^2));


% Singularity at ±fd: replace any Inf/NaN with the max finite value
S(~isfinite(S)) = max(S(isfinite(S)));

% Normalize so that sum(S)*df = 1  (unit power)
S = S / (sum(S) * df);

%% ---- Step 3: Complex Gaussian in frequency domain ----
rng(42);
X = sqrt(N) * (randn(1,N) + 1j*randn(1,N)) / sqrt(2);   % X[k] ~ CN(0,1)

%% ---- Step 4: Shape spectrum by sqrt(Jake's PSD) ----
Y = X .* sqrt(S * df);    % E[|Y[k]|^2] = S[k]*df

%% ---- Step 5: IFFT → complex fading envelope ----
% ifftshift: reorder from [-N/2..N/2-1] to [0..N-1] before IFFT
h = N * ifft(ifftshift(Y));    % N factor restores Parseval scaling

% Normalize to unit average power
h = h / sqrt(mean(abs(h).^2));
h_I = real(h);
h_Q = imag(h);

%% ---- Plot: Time-domain fading realization ----
t_ms = (0:N-1) * Ts * 1e3;

figure('Name', "Smith's FFT Fading — v=30 m/s");

% --- Subplot 1: Magnitude in dB ---
subplot(4,1,1);
plot(t_ms, 20*log10(abs(h)), 'b');
xlabel('Time (ms)'); ylabel('Instantaneous power |{\it h}[n]|^2 (dB)', 'Interpreter', 'tex');
title_str = sprintf('Path Gain Envelope of the chnanel with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', v, fd);
title(title_str, 'Interpreter', 'tex');
grid on;

subplot(4,1,2);
plot(t_ms, angle(h) * 180/pi, 'm');
xlabel('Time (ms)'); 
ylabel('Phase (deg)');
title_str_phase2 = sprintf('Instantaneous Phase of the channel with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', v, fd);
title(title_str_phase2, 'Interpreter', 'tex');
grid on;

% --- Subplot 2: In-phase Component ---
subplot(4,1,3);
plot(t_ms, h_I, 'r');
xlabel('Time (ms)'); ylabel('Amplitude');
title_str2 = sprintf('Real part (I component) of the channel with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', v, fd); 
title(title_str2, 'Interpreter', 'tex');
grid on;

% --- Subplot 3: Quadrature Component ---
subplot(4,1,4);
plot(t_ms, h_Q, 'g');
xlabel('Time (ms)'); ylabel('Amplitude');
title_str3 = sprintf('Imaginary part (Q component) of the channel with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', v, fd); 
title(title_str3, 'Interpreter', 'tex');
grid on;

%% Figure 2: ACF and Cross-correlation
Ns   = length(h);
Nlag     = min(20000, Ns-1);                               % lag range to display
[acf, lags] = xcorr(h,   Nlag, 'normalized');            % ACF of complex h
[xc,  ~   ] = xcorr(h_I, h_Q, Nlag, 'normalized');      % cross-corr I vs Q

tau_ms   = lags * Ts * 1e3;                              % lag axis in ms
R_theory = besselj(0, 2*pi * fd * abs(lags) * Ts);      % theoretical: J0(2*pi*fd*tau)

figure('Name', 'ACF and Cross-Correlation: Smith''s Method');

% --- Subplot 1: Autocorrelation (ACF) ---
subplot(1,2,1);
plot(tau_ms, real(acf), 'b', 'LineWidth', 1.2); hold on;
plot(tau_ms, R_theory,  'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex'); 
ylabel('{\itR_h}[{\it\tau}]', 'Interpreter', 'tex');
title_str1 = sprintf('Autocorrelation (ACF) of the N=%d process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', N, v, fd);
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
title_str2 = sprintf('Cross-Correlation (I vs Q) of the N=%d process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', N, v, fd);
title(title_str2, 'Interpreter', 'tex');
legend('Simulated', 'Ideal (0)', 'Interpreter', 'tex');
grid on;
text(0.02, 0.95, sprintf('RMS = %.4f', rms_xc), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'BackgroundColor', 'k', 'EdgeColor', 'k');


%% ---- Second Independent Path (different random draw) ----
% Remark 2 in Zheng & Xiao: using distinct (theta_k, phi_nk, varphi_nk)
% creates uncorrelated faders with identical statistical properties

rng(99);   % different seed -> independent realization
X2 = (randn(1,N) + 1j*randn(1,N)) / sqrt(2);   % X[k] ~ CN(0,1)
Y2 = X2 .* sqrt(S * df);    % E[|Y[k]|^2] = S[k]*df
h2 = N * ifft(ifftshift(Y2));
h2 = h2 / sqrt(mean(abs(h2).^2));


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
title_str_p = sprintf('Two Independent Path Gains of 2 N=%d processes with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', N, v, fd);
title(title_str_p, 'Interpreter', 'tex');
legend('Path 1', 'Path 2');
grid on;

subplot(1,2,2);
plot(tau_ms_p, xc_paths_real, 'w', 'LineWidth', 1.2);  hold on;
yline(0, 'r--', 'LineWidth', 1.2);
xlabel('Time Lag, {\it\tau} (ms)', 'Interpreter', 'tex');
ylabel('R_{h_1 h_2}[\tau]', 'Interpreter', 'tex');
title_str_xc = sprintf('Cross-Correlation of Path 1 vs Path 2, each of a N=%d process with {\\itv} = %d m/s or {\\itf}_d = %.0f Hz with Smith''s Method', N, v, fd);
title(title_str_xc, 'Interpreter', 'tex');
legend('Simulated', 'Ideal (0)', 'Interpreter', 'tex');
grid on;
text(0.02, 0.95, sprintf('RMS = %.4f', rms_paths), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'Color', 'w', 'BackgroundColor', 'k', 'EdgeColor', 'k');