%% Problem 3: Time-Varying Frequency-Selective Fading — 3D Plot
% Same 6-path PDP as Q2, fd = 90 Hz per path
% 7 snapshots at t = 0, 5, 10, 15, 20, 25, 30 ms

clear; clc; close all;

%% ---- PDP Parameters (identical to Q2) ----
W_bw     = 10e6;              % 2W = 10 MHz
Ts_tap   = 1/W_bw;            % tap spacing = 0.1 us
NFFT     = 2048;

gain_dB  = [-2, 0, -1, -6, -9, -14];
delay_us = [0, 1.8, 3.5, 5.7, 8.1, 12.3];
L        = length(gain_dB);

% Linear variances, normalized to sum = 1
sigma2   = 10.^(gain_dB/10);
sigma2   = sigma2 / sum(sigma2);

% Tap indices (delay / tap_spacing)
tap_idx  = round(delay_us / (Ts_tap * 1e6));

fprintf('Normalized sigma^2 per path:\n');
for i = 1:L
    fprintf('  Path %d: tap=%3d, sigma^2=%.4f\n', i, tap_idx(i), sigma2(i));
end

%% ---- Doppler and Snapshot Parameters ----
fd        = 90;                      % max Doppler per path (Hz)
wd        = 2*pi*fd;
M         = 20;                      % SOS sinusoids per path (Zheng-Xiao)
t_snap_ms = 0:5:30;                  % [0 5 10 15 20 25 30] ms
t_snap    = t_snap_ms * 1e-3;        % in seconds
Nsnap     = length(t_snap);          % 7 snapshots

fprintf('Coherence time Tc ~ 1/fd = %.1f ms\n', 1/fd*1e3);
fprintf('Snapshot spacing = 5 ms  (= %.1f * Tc)\n\n', 5e-3 * fd);

%% ---- Fix SOS random parameters for each path (drawn once) ----
% Each path has its own independent fading process
rng(42);
theta_all  = unifrnd(-pi, pi, 1, L);      % (1 x L)
phi_all    = unifrnd(-pi, pi, M, L);      % (M x L)
varphi_all = unifrnd(-pi, pi, M, L);      % (M x L)
n_vec      = (1:M)';                      % (M x 1)

%% ---- Frequency axis for plotting ----
f_MHz = (0:NFFT-1) * (W_bw/NFFT) / 1e6;  % 0 to 10 MHz

%% ---- Compute H[k, t] for all snapshots ----
% H_dB: (NFFT x Nsnap) matrix of |H(f, t)|^2 in dB
H_dB = zeros(NFFT, Nsnap);

for s = 1:Nsnap
    t_now  = t_snap(s);
    h_snap = zeros(1, NFFT);     % sparse impulse response at this snapshot

    for i = 1:L
        %% Zheng-Xiao SOS: evaluate path i at time t_now
        % Doppler arrival angles (conditionally random via theta)
        alpha_n = (2*pi*n_vec - pi + theta_all(i)) / (4*M);   % (M x 1)

        % Quadrature components (scalar values at t_now)
        Z_c = sqrt(2/M) * sum(cos(wd * t_now * cos(alpha_n) + phi_all(:,i)));
        Z_s = sqrt(2/M) * sum(cos(wd * t_now * sin(alpha_n) + varphi_all(:,i)));

        % Scale to CN(0, sigma2_i):
        % E[Z_c^2]=1, E[Z_s^2]=1  →  sqrt(sigma2/2) gives correct variance
        a_i = sqrt(sigma2(i)/2) * (Z_c + 1j*Z_s);

        % Accumulate at correct delay tap
        h_snap(tap_idx(i) + 1) = h_snap(tap_idx(i) + 1) + a_i;
    end

    % 2048-pt FFT (zero-padding already implicit since NFFT >> max tap)
    H_snap       = fft(h_snap, NFFT);
    H_dB(:, s)   = 10*log10(abs(H_snap).^2);
end

%% ---- 3D Surface Plot ----
% meshgrid: T_grid and F_grid both (NFFT x Nsnap)
[T_grid, F_grid] = meshgrid(t_snap_ms, f_MHz);

figure('Name', 'Q3: Time-Varying Frequency Response — 3D Surface');
surf(T_grid, F_grid, H_dB, 'EdgeColor', 'none');
colormap jet;
colorbar;
xlabel('Time  (ms)');
ylabel('Frequency  (MHz)');
zlabel('10 log_{10} |H(f,t)|^2  (dB)');
title({'Time-Varying Frequency Response of fading channel described by the given PDP and f_d = 90 Hz per path'});
view(45, 35);
grid on;

%% ---- 2D Overlay Plot (waterfall style — easier to read) ----
figure('Name', 'Q3: Frequency Snapshots Overlaid');
hold on;
cmap = lines(Nsnap);
for s = 1:Nsnap
    plot(f_MHz, H_dB(:,s), 'Color', cmap(s,:), 'LineWidth', 1.2);
end
xlabel('Frequency  (MHz)');
ylabel('10 log_{10} |H(f,t)|^2  (dB)');
title({'Time-Varying Frequency Response of fading channel at each snapshot time with f_d = 90 Hz per path'});
legend_str = arrayfun(@(x) sprintf('t = %d ms', x), t_snap_ms, 'UniformOutput', false);
legend(legend_str, 'Location', 'best');
grid on;

%% ===================================================================
%% BELLO SYSTEM FUNCTIONS — All 4 Domains
%% ===================================================================
% Recompute storing complex matrices (same SOS params, no new rng draw)
h_all     = zeros(NFFT, Nsnap);   % h(tau, t): delay-time plane
H_complex = zeros(NFFT, Nsnap);   % H(f,   t): frequency-time plane

for s = 1:Nsnap
    t_now = t_snap(s);
    h_tmp = zeros(1, NFFT);
    for i = 1:L
        alpha_n = (2*pi*n_vec - pi + theta_all(i)) / (4*M);
        Z_c = sqrt(2/M) * sum(cos(wd * t_now * cos(alpha_n) + phi_all(:,i)));
        Z_s = sqrt(2/M) * sum(cos(wd * t_now * sin(alpha_n) + varphi_all(:,i)));
        a_i = sqrt(sigma2(i)/2) * (Z_c + 1j*Z_s);
        h_tmp(tap_idx(i)+1) = h_tmp(tap_idx(i)+1) + a_i;
    end
    h_all(:, s)     = h_tmp.';
    tmp             = fft(h_tmp, NFFT);
    H_complex(:, s) = tmp.';
end

%% ---- Derived Bello functions via FFT over time axis ----
% S(f, nu)  = FT_t{ H(f,t) } — each row of H_complex FFT'd over Nsnap points
S_fn = fftshift(fft(H_complex, Nsnap, 2), 2);    % (NFFT x Nsnap)

% B(tau,nu) = FT_t{ h(tau,t) } — each row of h_all FFT'd over Nsnap points
B_tn = fftshift(fft(h_all,     Nsnap, 2), 2);    % (NFFT x Nsnap)

%% ---- Axes ----
% Delay axis — only show meaningful taps (0 to last path + margin)
tau_max_disp = max(tap_idx) + 5;
tau_us_disp  = (0:tau_max_disp) * Ts_tap * 1e6;  % microseconds
tau_rows     = 1 : tau_max_disp+1;               % 1-indexed

% Doppler axis from FFT of 7 snapshots spaced dt=5ms
dt_snap = 5e-3;
dnu     = 1 / (Nsnap * dt_snap);                          % ~28.6 Hz resolution
nu_Hz   = (-floor(Nsnap/2) : ceil(Nsnap/2)-1) * dnu;     % centered

fprintf('\n-- Bello Function Axes --\n');
fprintf('  Doppler resolution : dnu = %.2f Hz\n', dnu);
fprintf('  Doppler range      : [%.1f, %.1f] Hz\n', nu_Hz(1), nu_Hz(end));
fprintf('  Delay range shown  : 0 to %.1f us (%d taps)\n\n', ...
        tau_us_disp(end), tau_max_disp+1);

%% ===== FIGURE 3: Bello Domain 1 — h(tau, t) =====
[T3, TAU3] = meshgrid(t_snap_ms, tau_us_disp);
h_dB3      = 10*log10(abs(h_all(tau_rows, :)).^2 + 1e-12);

figure('Name', 'Bello 1: h(tau,t) — Delay-Time Plane');
surf(T3, TAU3, h_dB3, 'EdgeColor', 'none');
colormap jet; colorbar;
xlabel('Time  t  (ms)');
ylabel('Delay  \tau  (\mus)');
zlabel('|h(\tau,t)|^2  (dB)');
title({'Bello Function 1:  Input Delay-Spread Function  h(\tau, t)', ...
       'Direct time-varying impulse response | f_d = 90 Hz | 6-path PDP'});
view(45, 35); grid on;

%% ===== FIGURE 4: Bello Domain 2 — H(f, t) =====
[T4, F4] = meshgrid(t_snap_ms, f_MHz);
H_dB4    = 10*log10(abs(H_complex).^2 + 1e-12);

figure('Name', 'Bello 2: H(f,t) — Frequency-Time Plane');
surf(T4, F4, H_dB4, 'EdgeColor', 'none');
colormap jet; colorbar;
xlabel('Time  t  (ms)');
ylabel('Frequency  f  (MHz)');
zlabel('|H(f,t)|^2  (dB)');
title({'Bello Function 2:  Time-Varying Transfer Function  H(f, t)', ...
       'FT_\tau\{h(\tau,t)\} | f_d = 90 Hz | 6-path PDP'});
view(45, 35); grid on;

%% ===== FIGURE 5: Bello Domain 3 — S(f, nu) =====
[NU5, F5] = meshgrid(nu_Hz, f_MHz);
S_dB5     = 10*log10(abs(S_fn).^2 + 1e-12);

figure('Name', 'Bello 3: S(f,nu) — Doppler-Spread Function');
surf(NU5, F5, S_dB5, 'EdgeColor', 'none');
colormap jet; colorbar;
xlabel('Doppler  \nu  (Hz)');
ylabel('Frequency  f  (MHz)');
zlabel('|S(f,\nu)|^2  (dB)');
title({'Bello Function 3:  Doppler-Spread Function  S(f, \nu)', ...
       'FT_t\{H(f,t)\} | dnu = 28.6 Hz | f_d = 90 Hz per path'});
view(45, 35); grid on;

%% ===== FIGURE 6: Bello Domain 4 — B(tau, nu) =====
[NU6, TAU6] = meshgrid(nu_Hz, tau_us_disp);
B_dB6      = 10*log10(abs(B_tn(tau_rows, :)).^2 + 1e-12);

figure('Name', 'Bello 4: B(tau,nu) — Delay-Doppler Spread Function');
surf(NU6, TAU6, B_dB6, 'EdgeColor', 'none');
colormap jet; colorbar;
xlabel('Doppler  \nu  (Hz)');
ylabel('Delay  \tau  (\mus)');
zlabel('|B(\tau,\nu)|^2  (dB)');
title({'Bello Function 4:  Delay-Doppler Spread Function  B(\tau, \nu)', ...
       'FT_t\{h(\tau,t)\} | dnu = 28.6 Hz | f_d = 90 Hz per path'});
view(45, 35); grid on;