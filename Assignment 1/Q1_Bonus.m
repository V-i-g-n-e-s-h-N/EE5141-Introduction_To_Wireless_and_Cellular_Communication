%% Flexible Doppler Behavior Model — von Mises PAS + Modified SoS
% =========================================================================
%  Reference: IEEE C802.16m-07/140
%  "A Flexible Doppler Behavior Model"
%  Shuangquan Wang, Xiaodong Wang, Mohammad Madihian (2007)
%
%  Implements the von Mises Doppler model using the Modified Sum of 
%  Sinusoids (SoS) method (Zheng & Xiao, 2002).
%
%  The classical Jakes' model is the special case where kappa = 0.
% =========================================================================

clear; clc; close all;

%% ---- System Parameters ----
v   = 30;           % mobile speed (m/s)
fc  = 2e9;          % carrier frequency (Hz)
c   = 3e8;          % speed of light (m/s)
W   = 50e3;         % one-sided bandwidth (Hz)
Ts  = 1/(2*W);      % sampling period (s)
fd  = v*fc/c;       % max Doppler frequency (Hz)
wd  = 2*pi*fd;      % angular Doppler (rad/s)
N   = 8192;         % number of time samples
fs  = 1/Ts;         % sampling frequency (Hz)

%% ---- von Mises PAS Parameters (Eq. 8, single cluster N=1) ----
%  p(theta) = exp[kappa * cos(theta - phi)] / (2*pi*I0(kappa))
%  kappa = 0  ->  uniform PAS  ->  Jakes' model (special case)
%  kappa > 0  ->  non-isotropic scattering
kappa = 5;           % concentration parameter (angle spread control)
phi   = pi/4;        % mean AoA relative to mobile direction (rad)

%% ---- Common vectors ----
t    = (0:N-1) * Ts;
t_ms = t * 1e3;
Nlag = min(4000, N-1); % Reduced lag length for cleaner plotting

% =========================================================================
%% FIGURE 1: Power Azimuth Spectrum (PAS) Comparison
% =========================================================================
theta_ax = linspace(-pi, pi, 1000);
pas_vm   = exp(kappa * cos(theta_ax - phi)) / (2*pi*besseli(0, kappa));
pas_unif = ones(size(theta_ax)) / (2*pi);

figure('Name', 'Fig 1: PAS Comparison');
plot(theta_ax*180/pi, pas_vm,   'b',  'LineWidth', 1.5); hold on;
plot(theta_ax*180/pi, pas_unif, 'r--','LineWidth', 1.2);
xlabel('Azimuth Angle \theta (deg)');
ylabel('p(\theta)');
title(sprintf('Power Azimuth Spectrum of von Mises (\\kappa=%.1f, \\phi=%.0f°) vs Uniform (Jakes)', kappa, phi*180/pi), 'Interpreter','tex');
legend('von Mises', 'Uniform (Jakes)', 'Location', 'best');
grid on;

% =========================================================================
%% FIGURE 2: Doppler Spectrum Comparison — Eq. (10) vs Eq. (7)
% =========================================================================
df    = fs / N;
f_ax  = (-N/2 : N/2-1) * df;
idx   = abs(f_ax) < fd;
f_in  = f_ax(idx);
ratio = f_in / fd;
sq1mr = sqrt(1 - ratio.^2);

% von Mises Doppler spectrum — Eq. (10)
S_vm = zeros(1, N);
S_vm(idx) = exp(kappa * ratio * cos(phi)) .* ...
            cosh(kappa * sq1mr * sin(phi)) ./ ...
            (pi * fd * sq1mr * besseli(0, kappa));
S_vm(~isfinite(S_vm)) = max(S_vm(isfinite(S_vm)));
S_vm_norm = S_vm / (sum(S_vm)*df);          % normalize to unit power

% Jakes' Doppler spectrum — Eq. (7)
S_jk = zeros(1, N);
S_jk(idx) = 1 ./ (pi * fd * sqrt(1 - (f_in/fd).^2));
S_jk(~isfinite(S_jk)) = max(S_jk(isfinite(S_jk)));
S_jk_norm = S_jk / (sum(S_jk)*df);

figure('Name', 'Fig 2: Doppler Spectrum');
plot(f_ax, S_vm_norm, 'b',  'LineWidth', 1.5); hold on;
plot(f_ax, S_jk_norm, 'r--','LineWidth', 1.2);
xlim([-1.3*fd, 1.3*fd]);
xlabel('Frequency f (Hz)');  ylabel('S_g(f)');
title(sprintf('Doppler Spectrum of von Mises (\\kappa=%.1f, \\phi=%.0f°) vs Jakes', kappa, phi*180/pi), 'Interpreter','tex');
legend('von Mises (Eq. 10)', 'Jakes bathtub (Eq. 7)', 'Location', 'best');
grid on;

% =========================================================================
%% THEORETICAL ACFs for comparison
% =========================================================================
lags_vec   = -Nlag:Nlag;
signed_tau = lags_vec * Ts;          % signed for correct Im{r_g} symmetry
tau_abs    = abs(lags_vec) * Ts;
tau_ms     = lags_vec * Ts * 1e3;

% von Mises ACF — Eq. (9):
% r_g(tau) = I0(sqrt(kappa^2 - 4*pi^2*fd^2*tau^2 + j*4*pi*kappa*fd*tau*cos(phi))) / I0(kappa)
arg_vm  = sqrt(kappa^2 - 4*pi^2*fd^2*signed_tau.^2 + ...
               1j*4*pi*kappa*fd*signed_tau*cos(phi));
R_vm_th = besseli(0, arg_vm) / besseli(0, kappa);

% Jakes ACF — Eq. (6):  r_g(tau) = J0(2*pi*fd*|tau|)
R_jk_th = besselj(0, 2*pi*fd*tau_abs);


% =========================================================================
%% MODIFIED SUM-OF-SINUSOIDS (Zheng & Xiao) + von Mises Generation
% =========================================================================
M_sos = 20; % Number of sinusoids

% 1. Draw arrival angles from von Mises(phi, kappa) instead of uniform
rng(42);
alpha_n = vmrnd(phi, kappa, 1, M_sos);

% 2. Generate TWO independent sets of random phases ~ Uniform[-pi, pi)
phi_n    = unifrnd(-pi, pi, 1, M_sos);  % For In-phase
varphi_n = unifrnd(-pi, pi, 1, M_sos);  % For Quadrature

% 3. Calculate Z_c (In-phase) and Z_s (Quadrature) using Zheng & Xiao formulation
t_col = t(:); % column vector (N x 1) for matrix expansion
Z_c = sqrt(2/M_sos) * sum(cos(wd * t_col * cos(alpha_n) + phi_n), 2).'; % (1 x N)
Z_s = sqrt(2/M_sos) * sum(cos(wd * t_col * sin(alpha_n) + varphi_n), 2).'; % (1 x N)

% 4. Combine into complex baseband channel
h_sos = Z_c + 1j * Z_s;

% Normalize to ensure unit average power
h_sos = h_sos / sqrt(mean(abs(h_sos).^2));

% Extract normalized components for plotting
Z_c = real(h_sos);
Z_s = imag(h_sos);

% ---- Figure 3: Time-domain — SoS ----
figure('Name', 'Fig 3: SoS — von Mises Fading');
subplot(4,1,1);
plot(t_ms, 20*log10(abs(h_sos)), 'b');
xlabel('Time (ms)'); ylabel('|h[n]| (dB)');
title(sprintf('Path Gain of Modified SoS (M=%d), von Mises (\\kappa=%.1f, \\phi=%.0f°), f_d=%.0f Hz', M_sos, kappa, phi*180/pi, fd), 'Interpreter','tex');
grid on;
subplot(4,1,2);
plot(t_ms, angle(h_sos)*180/pi, 'm');
xlabel('Time (ms)'); ylabel('Phase (deg)');
title(sprintf('Instantaneous Phase of Modified SoS (M=%d), von Mises (\\kappa=%.1f, \\phi=%.0f°), f_d=%.0f Hz', M_sos, kappa, phi*180/pi, fd), 'Interpreter','tex');
grid on;
subplot(4,1,3);
plot(t_ms, Z_c, 'r');
xlabel('Time (ms)'); ylabel('Amplitude');
title(sprintf('Real Part of Modified SoS (M=%d), von Mises (\\kappa=%.1f, \\phi=%.0f°), f_d=%.0f Hz', M_sos, kappa, phi*180/pi, fd), 'Interpreter','tex');
grid on;
subplot(4,1,4);
plot(t_ms, Z_s, 'g');
xlabel('Time (ms)'); ylabel('Amplitude');
title(sprintf('Imaginary Part of Modified SoS (M=%d), von Mises (\\kappa=%.1f, \\phi=%.0f°), f_d=%.0f Hz', M_sos, kappa, phi*180/pi, fd), 'Interpreter','tex');
grid on;

% =========================================================================
%% LOCAL FUNCTION: von Mises random variate generator
%  Best & Fisher (1979) rejection algorithm — Ref [10] of the paper
% =========================================================================
function theta = vmrnd(mu, kap, varargin)
%VMRND  Generate von Mises distributed random variates.
%   theta = vmrnd(mu, kappa, m, n)
%   mu    : mean direction (rad)
%   kappa : concentration (>= 0);  kappa=0 -> uniform on [-pi,pi)
    if kap < 1e-6
        theta = -pi + 2*pi*rand(varargin{:});
        return;
    end
    sz   = [varargin{:}];
    nTot = prod(sz);

    tau_bf = 1 + sqrt(1 + 4*kap^2);
    rho    = (tau_bf - sqrt(2*tau_bf)) / (2*kap);
    r_bf   = (1 + rho^2) / (2*rho);

    theta = zeros(1, nTot);
    cnt   = 0;
    while cnt < nTot
        u1 = rand;
        z  = cos(pi * u1);
        f  = (1 + r_bf*z) / (r_bf + z);
        cc = kap * (r_bf - f);
        u2 = rand;
        if cc*(2 - cc) > u2 || log(cc/u2) + 1 >= cc
            u3  = rand;
            cnt = cnt + 1;
            theta(cnt) = sign(u3 - 0.5) * acos(f);
        end
    end
    theta = mu + reshape(theta, sz);
    theta = mod(theta + pi, 2*pi) - pi;      % wrap to [-pi, pi]
end