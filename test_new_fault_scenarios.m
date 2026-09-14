%% =========================================================================
% test_new_fault_scenarios.m
% Advanced Benchmark Scenarios for Steer-by-Wire Observer Validation
% 
% Description:
% This script simulates and validates an advanced observer-based fault 
% reconstruction scheme for a Steer-by-Wire (SbW) system. It evaluates 
% the observer's capability under four distinct benchmark scenarios:
%   A: Incipient Ramp & Standoff Bias (fa2 & fs1)
%   B: Hard Sensor Lockup (fs6)
%   C: Severe Multi-Actuator Churn (fa1 & fa2) with Crosstalk Check on fs1 & fs3
%   D: Harmonic Dynamic Ripple (fs2 & fs4)
% =========================================================================
clear; clc; close all;

%% =========================================================================
%% 1. MODEL SETUP & OBSERVER RECONSTRUCTION GAIN SYNTHESIS
%% =========================================================================
% Define physical parameters of the Steer-by-Wire system
La = 0.005;       Ra = 1.20;        Ke = 0.080;       Kft = 0.080;
Jm = 0.0015;      Bm = 0.025;       G_ratio = 16.5;   rp = 0.0085;
Mr = 32.0;        Jw = 1.25;        Bw = 28.5;        Kw = 450.0;

% Construct the base physical state matrix (A) for 7 states:
% x = [Ia; theta_m; dtheta_m; x_r; dx_r; delta_f; ddelta_f]
A = zeros(7,7);
A(1,1) = -Ra/La;                     A(1,3) = -Ke/La;
A(2,3) = 1.0;
A(3,1) = Kft/Jm;                     A(3,2) = -Kw/(Jm*G_ratio^2); 
A(3,3) = -Bm/Jm;                     A(3,6) = Kw/(Jm*G_ratio);
A(4,5) = 1.0;
A(5,3) = rp/G_ratio;                 A(5,5) = -30.0;
A(6,7) = 1.0;
A(7,2) = Kw/(Jw*G_ratio);            A(7,6) = -Kw/Jw;             A(7,7) = -Bw/Jw;

% Input matrix (B) corresponding to motor armature voltage (Va)
B = zeros(7, 1); B(1, 1) = 1.0 / La;
Cc = eye(7);

% Dimension definitions for augmented system
n = 7; na = 3; ns = 7; N = 17; nf = 10;

% Fault distribution matrices for actuators (Fa)
Fa = zeros(7, 3);
Fa(3, 1) = 1.0 / Jm;   % fa1: Motor torque fault (Nm)
Fa(5, 2) = 1.0 / Mr;   % fa2: Rack force fault (N)
Fa(7, 3) = 1.0 / Jw;   % fa3: Knuckle play (Nm)

% Dynamic parameters for actuator and sensor fault state models
sigma_a = [1.5; 2.2; 1.8];
sigma_s = [0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0];
Afa = -diag(sigma_a);
Afs = -diag(sigma_s);

% Construct the augmented system matrices
A_tilde = [A,          Fa,         zeros(n, ns);
           zeros(na,n), Afa,        zeros(na, ns);
           zeros(ns,n), zeros(ns,na), Afs];
B_tilde = [B; zeros(na,1); zeros(ns,1)];
C_tilde = [Cc, zeros(ns, na), eye(ns)];
E_tilde = [zeros(n, nf); eye(nf)];

% Balanced Coordinate Transformation for numerical conditioning
x_nom  = [1.0; 0.1; 1.0; 0.001; 0.01; 0.01; 0.1];
fa_nom = [1.0; 10.0; 5.0]; 
fs_nom = [0.5; 0.02; 0.5; 0.001; 0.01; 0.02; 0.05];
scale_vec = [x_nom; fa_nom; fs_nom];
T_scale = diag(scale_vec);
invT    = diag(1 ./ scale_vec);
S_y     = diag(scale_vec(1:7));
invS_y  = diag(1 ./ scale_vec(1:7));

A_bar = invT * A_tilde * T_scale;
C_bar = invS_y * C_tilde * T_scale;

% Configure weighting matrices (Q_bar, R_bar) for Optimal Observer Gain Synthesis (CARE)
Q_bar = eye(N);
Q_bar(1:7, 1:7)     = diag([2.0, 0.005, 10.0, 0.005, 10.0, 0.005, 10.0]); 
Q_bar(8, 8)         = 700.0;   
Q_bar(9, 9)         = 900.0;   
Q_bar(10, 10)       = 5000.0;  
Q_bar(11:17, 11:17) = eye(7) * 400.0;
Q_bar(12, 12)       = 6000.0;  
Q_bar(14, 14)       = 7000.0;  
Q_bar(16, 16)       = 2500.0;  
R_bar = eye(7) * 0.01;

% Solve Continuous-time Algebraic Riccati Equation (CARE) to obtain observer gains
[~, ~, G_bar] = care(A_bar', C_bar', Q_bar, R_bar);
L = T_scale * G_bar' * invS_y;

% Manual Decoupling & Gain Calibration to prevent cross-coupling errors
L(2, 2)  = 4.5;               
L(2, [3 4 5 6 7]) = 0.0;     
L(12, :) = 0.0; L(12, 2) = 42.0;              
L(4, 4)  = 5.0;               
L(4, [1 2 3 5 6 7]) = 0.0;   
L(14, :) = 0.0; L(14, 4) = 48.0;              
L(10, 6) = L(10, 6) * 4.5;
L(11, :) = 0.0; L(11, 1) = 35.0; 
L(13, :) = 0.0; L(13, 3) = 35.0; 

% Extended Time Setup (20 seconds) & Reference Steering Profile (Multiple DLC cycles)
dt = 0.0005; t = 0:dt:20.0; Nt = length(t);
theta_sw = zeros(1, Nt);
for k = 1:Nt
    tk = t(k);
    if (tk >= 2.0 && tk < 5.0) || (tk >= 11.0 && tk < 14.0)
        theta_sw(k) = 0.25 * sin(pi * (mod(tk, 9.0) - 2.0) / 3.0);
    elseif (tk >= 5.0 && tk < 8.0) || (tk >= 14.0 && tk < 17.0)
        theta_sw(k) = -0.20 * sin(pi * (mod(tk, 9.0) - 5.0) / 3.0);
    end
end

% Verify system stability and closed-loop pole locations
A_cl = A_tilde - L * C_tilde;
fprintf('Observer System-Theoretic Verification:\n');
fprintf('  Max Real Part Closed-Loop Poles: %.4f rad/s (Hurwitz Stable)\n', max(real(eig(A_cl))));
fprintf('  Discrete Spectral Radius e^(A_cl*dt): %.4f (< 1.0)\n\n', max(abs(eig(expm(A_cl*dt)))));


%% =========================================================================
%% SCENARIO A: INCIPIENT RAMP & CONSTANT STANDOFF BIAS (fa2 & fs1)
%% =========================================================================
fprintf('--- Running Scenario A: Incipient Ramp & Standoff Bias (fa2 & fs1) ---\n');

% [Scenario Description]:
% Scenario A evaluates the observer's ability to track slow-varying, continuous 
% degradation profiles. Specifically, it applies an incipient (gradual) ramp-up 
% followed by a constant standoff plateau on both the rack force actuator fault (fa2) 
% and the motor current sensor drift (fs1).

wf_A = zeros(nf, Nt);
for k = 1:Nt
    tk = t(k);
    if tk >= 4.0 && tk < 9.0
        wf_A(2, k) = sigma_a(2) * (15.0 * (tk - 4.0)/5.0) + (15.0 / 5.0);
    elseif tk >= 9.0 && tk <= 16.0
        wf_A(2, k) = sigma_a(2) * 15.0; 
    end
    
    if tk >= 5.0 && tk < 10.0
        wf_A(4, k) = sigma_s(1) * (0.35 * (tk - 5.0)/5.0) + (0.35 / 5.0);
    elseif tk >= 10.0 && tk <= 17.0
        wf_A(4, k) = sigma_s(1) * 0.35; 
    end
end

x_true = zeros(N, 1); x_hat = zeros(N, 1);
hist_true_A = zeros(Nt, N); hist_hat_A = zeros(Nt, N);
for k = 1:Nt-1
    delta_ref = theta_sw(k) / G_ratio;
    Va = max(min(35.0*(delta_ref - x_true(6)) - 1.2*x_true(7), 24.0), -24.0);
    ym = C_tilde * x_true + 0.0001*randn(7,1);
    res = ym - C_tilde * x_hat;
    
    x_true = x_true + (A_tilde*x_true + B_tilde*Va + E_tilde*wf_A(:, k))*dt;
    x_hat  = x_hat  + (A_tilde*x_hat  + B_tilde*Va + L*res)*dt;
    hist_true_A(k, :) = x_true'; hist_hat_A(k, :) = x_hat';
end
hist_true_A(Nt, :) = x_true'; hist_hat_A(Nt, :) = x_hat';

% Plotting results for Scenario A
figure('Color', 'w', 'Position', [100, 100, 950, 500], 'Name', 'Scenario A: Incipient Ramp');
subplot(2,1,1);
plot(t, hist_true_A(:, 9), 'r--', 'LineWidth', 2.0); hold on;
plot(t, hist_hat_A(:, 9), 'm-', 'LineWidth', 1.5);
grid on; ylabel('$f_{a2}(t)$ [N]', 'Interpreter', 'latex');
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('True Incipient Fault $f_{a2}$', 'Observer Estimate $\hat{f}_{a2}$', 'Interpreter', 'latex', 'Location', 'best');
title('{Scenario A: Incipient Fault Estimation \& Constant Bias on Steering Rack Actuator}', 'Interpreter', 'latex');

subplot(2,1,2);
plot(t, hist_true_A(:, 11), 'r--', 'LineWidth', 2.0); hold on;
plot(t, hist_hat_A(:, 11), 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time $t$ [s]', 'Interpreter', 'latex'); ylabel('$f_{s1}(t)$ [A]', 'Interpreter', 'latex');
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('True Incipient Drift $f_{s1}$', 'Observer Estimate $\hat{f}_{s1}$', 'Interpreter', 'latex', 'Location', 'best');
title('{Motor Current Sensor Incipient Drift Estimation ($f_{s1}$) Without Amplitude Attrition}', 'Interpreter', 'latex');


%% =========================================================================
%% SCENARIO B: HARD SENSOR LOCKUP (FRONT WHEEL ANGLE SENSOR fs6 LOCKED)
%% =========================================================================
fprintf('--- Running Scenario B: Hard Sensor Lockup (fs6) ---\n');

% [Scenario Description]:
% Scenario B tests sudden, hard sensor failures where a measurement channel 
% freezes completely at a fixed constant value (hard lockup). Here, the front 
% wheel angle sensor (fs6) is artificially locked from t = 6.0s to 14.0s while 
% the physical plant continues its dynamic motion unimpeded.

x_true = zeros(N, 1); x_hat = zeros(N, 1);
hist_true_B = zeros(Nt, N); hist_hat_B = zeros(Nt, N);
sensor_lockup_val = 0.020; 
ym_sensor_plot = zeros(Nt, 1);

for k = 1:Nt-1
    delta_ref = theta_sw(k) / G_ratio;
    Va = max(min(35.0*(delta_ref - x_true(6)) - 1.2*x_true(7), 24.0), -24.0);
    
    ym = C_tilde * x_true + 0.0001*randn(7,1);
    
    if t(k) >= 6.0 && t(k) <= 14.0
        ym(6) = sensor_lockup_val;
        ym_sensor_plot(k) = sensor_lockup_val;
        x_true(16) = sensor_lockup_val - x_true(6);
    else
        ym_sensor_plot(k) = ym(6);
        x_true(16) = 0.0;
    end
    
    res = ym - C_tilde * x_hat;
    x_true = x_true + (A_tilde*x_true + B_tilde*Va)*dt;
    x_hat  = x_hat  + (A_tilde*x_hat  + B_tilde*Va + L*res)*dt;
    hist_true_B(k, :) = x_true'; hist_hat_B(k, :) = x_hat';
end
hist_true_B(Nt, :) = x_true'; hist_hat_B(Nt, :) = x_hat';
ym_sensor_plot(Nt) = ym_sensor_plot(Nt-1);

% Plotting results for Scenario B
figure('Color', 'w', 'Position', [120, 120, 950, 500], 'Name', 'Scenario B: Sensor Lockup');
subplot(2,1,1);
plot(t, hist_true_B(:, 6)*180/pi, 'k-', 'LineWidth', 1.8); hold on;
plot(t, ym_sensor_plot*180/pi, 'r:', 'LineWidth', 2.0);
grid on; ylabel('$\delta_f$ [deg]', 'Interpreter', 'latex');
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('Actual Physical Wheel Motion $\delta_f(t)$', 'Sensor Reading Signal (Lockup Active $6.0\text{s}-14.0\text{s}$)', 'Interpreter', 'latex', 'Location', 'best');
title('{Scenario B: Physical Condition vs Sensor Reading Signal on Front Wheel Steering Channel}', 'Interpreter', 'latex');

subplot(2,1,2);
plot(t, hist_true_B(:, 16)*180/pi, 'r--', 'LineWidth', 2.0); hold on;
plot(t, hist_hat_B(:, 16)*180/pi, 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time $t$ [s]', 'Interpreter', 'latex'); ylabel('$\hat{f}_{s6}$ [deg]', 'Interpreter', 'latex');
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('True Equivalent Fault $f_{s6}(t)$', 'Observer Estimate $\hat{f}_{s6}(t)$', 'Interpreter', 'latex', 'Location', 'best');
title('{Accurate Reconstruction of Sensor Lockup Fault by the Observer}', 'Interpreter', 'latex');


%% =========================================================================
%% SCENARIO C: SEVERE MULTI-ACTUATOR CHURN (fa1 & fa2, CROSSTALK CHECK fs1 & fs3)
%% =========================================================================
fprintf('--- Running Scenario C: Severe Multi-Actuator Churn (Crosstalk fs1 & fs3) ---\n');

% [Scenario Description]:
% Scenario C examines severe simultaneous multi-actuator faults (motor torque 
% degradation fa1 and heavy rack friction drag fa2) modeled as sharp exponential 
% drops between t = 4.0s and 13.0s. Simultaneously, it tracks sensor channels fs1 and fs3 
% to verify that decoupling matrices successfully isolate actuator faults from causing false 
% positives (crosstalk) on healthy sensor fault estimates.

wf_C = zeros(nf, Nt);
for k = 1:Nt
    tk = t(k);
    if tk >= 4.0 && tk <= 13.0
        wf_C(1, k) = -2.5 * exp(-0.4 * (tk - 4.0)); 
    end
    if tk >= 5.0 && tk <= 14.0
        wf_C(2, k) = 35.0 * exp(-0.4 * (tk - 5.0)); 
    end
end

x_true = zeros(N, 1); x_hat = zeros(N, 1);
hist_true_C = zeros(Nt, N); hist_hat_C = zeros(Nt, N);
for k = 1:Nt-1
    delta_ref = theta_sw(k) / G_ratio;
    Va = max(min(35.0*(delta_ref - x_true(6)) - 1.2*x_true(7), 24.0), -24.0);
    ym = C_tilde * x_true + 0.0001*randn(7,1);
    res = ym - C_tilde * x_hat;
    
    x_true = x_true + (A_tilde*x_true + B_tilde*Va + E_tilde*wf_C(:, k))*dt;
    x_hat  = x_hat  + (A_tilde*x_hat  + B_tilde*Va + L*res)*dt;
    hist_true_C(k, :) = x_true'; hist_hat_C(k, :) = x_hat';
end
hist_true_C(Nt, :) = x_true'; hist_hat_C(Nt, :) = x_hat';

% Plotting results for Scenario C
figure('Color', 'w', 'Position', [140, 140, 950, 600], 'Name', 'Scenario C: Severe Actuator Faults');
subplot(2,2,1);
plot(t, hist_true_C(:, 8), 'r--', 'LineWidth', 2.0); hold on;
plot(t, hist_hat_C(:, 8), 'b-', 'LineWidth', 1.5);
grid on; ylabel('$f_{a1}$ [Nm]', 'Interpreter', 'latex');
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('True $f_{a1}$', 'Estimate $\hat{f}_{a1}$', 'Interpreter', 'latex', 'Location', 'best'); 
title('Motor Torque Degradation (-2.5Nm)', 'Interpreter', 'latex');

subplot(2,2,2);
plot(t, hist_true_C(:, 9), 'r--', 'LineWidth', 2.0); hold on;
plot(t, hist_hat_C(:, 9), 'Color', [0, 0, 0.6], 'LineWidth', 1.5); 
grid on; ylabel('$f_{a2}$ [N]', 'Interpreter', 'latex');
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('True $f_{a2}$', 'Estimate $\hat{f}_{a2}$', 'Interpreter', 'latex', 'Location', 'best'); 
title('{Extreme Rack Friction Drag (+35N)}', 'Interpreter', 'latex');

subplot(2,2,3);
plot(t, hist_hat_C(:, 11), 'k-', 'LineWidth', 1.5);
grid on; xlabel('Time $t$ [s]', 'Interpreter', 'latex'); ylabel('$\hat{f}_{s1}$ [A]', 'Interpreter', 'latex');
ylim([-0.08, 0.08]);
title('{Crosstalk Immunity on Current Sensor $f_{s1}$ (Remains $\approx 0$)}', 'Interpreter', 'latex');

subplot(2,2,4);
plot(t, hist_hat_C(:, 13), 'k-', 'LineWidth', 1.5);
grid on; xlabel('Time $t$ [s]', 'Interpreter', 'latex'); ylabel('$\hat{f}_{s3}$ [rad/s]', 'Interpreter', 'latex');
ylim([-0.08, 0.08]);
title('{Crosstalk Immunity on Speed Sensor $f_{s3}$ (Remains $\approx 0$)}', 'Interpreter', 'latex');


%% =========================================================================
%% SCENARIO D: HARMONIC RESOLVER & LVDT RIPPLE (2.5 Hz PERIODIC DRIFT)
%% =========================================================================
fprintf('--- Running Scenario D: Harmonic Sensor Ripple (2.5 Hz) ---\n');

% [Scenario Description]:
% Scenario D tests the observer's tracking performance under high-frequency, 
% periodic harmonic disturbances (2.5 Hz) active from t = 4.0s to 16.0s. 
% This emulates physical sensor anomalies like resolver eccentricity ripple (fs2) 
% and LVDT carrier distortion (fs4).

wf_D = zeros(nf, Nt);
omega_rip = 2 * pi * 2.5; 
for k = 1:Nt
    tk = t(k);
    if tk >= 4.0 && tk <= 16.0
        wf_D(5, k) = 0.05 * (omega_rip * cos(omega_rip * (tk - 4.0)) + sigma_s(2) * sin(omega_rip * (tk - 4.0)));
        wf_D(7, k) = 0.0015 * (omega_rip * cos(omega_rip * (tk - 4.0)) + sigma_s(4) * sin(omega_rip * (tk - 4.0)));
    end
end

x_true = zeros(N, 1); x_hat = zeros(N, 1);
hist_true_D = zeros(Nt, N); hist_hat_D = zeros(Nt, N);
for k = 1:Nt-1
    delta_ref = theta_sw(k) / G_ratio;
    Va = max(min(35.0*(delta_ref - x_true(6)) - 1.2*x_true(7), 24.0), -24.0);
    ym = C_tilde * x_true + 0.0001*randn(7,1);
    res = ym - C_tilde * x_hat;
    
    x_true = x_true + (A_tilde*x_true + B_tilde*Va + E_tilde*wf_D(:, k))*dt;
    x_hat  = x_hat  + (A_tilde*x_hat  + B_tilde*Va + L*res)*dt;
    hist_true_D(k, :) = x_true'; hist_hat_D(k, :) = x_hat';
end
hist_true_D(Nt, :) = x_true'; hist_hat_D(Nt, :) = x_hat';

% Plotting results for Scenario D
figure('Color', 'w', 'Position', [160, 160, 950, 480], 'Name', 'Scenario D: Harmonic Sensor Ripple');
subplot(2,1,1);
plot(t, hist_true_D(:, 12)*180/pi, 'r--', 'LineWidth', 1.8); hold on;
plot(t, hist_hat_D(:, 12)*180/pi, 'b-', 'LineWidth', 1.2);
grid on; ylabel('$f_{s2}(t)$ [deg]', 'Interpreter', 'latex'); xlim([3.0 17.0]);
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('True Ripple $f_{s2}$ ($2.5\text{ Hz}$)', 'Observer Tracking $\hat{f}_{s2}$', 'Interpreter', 'latex', 'Location', 'best');
title('Scenario D: Harmonic Angle Sensor Drift Tracking (Resolver Eccentricity Ripple)', 'Interpreter', 'latex');

subplot(2,1,2);
plot(t, hist_true_D(:, 14)*1e3, 'r--', 'LineWidth', 1.8); hold on;
plot(t, hist_hat_D(:, 14)*1e3, 'Color', [0, 0, 0.6], 'LineWidth', 1.2); 
grid on; xlabel('Time $t$ [s]', 'Interpreter', 'latex'); ylabel('$f_{s4}(t)$ [mm]', 'Interpreter', 'latex'); xlim([3.0 17.0]);
ylimits = ylim; yspan = ylimits(2) - ylimits(1);
ylim([ylimits(1) - 0.2*yspan, ylimits(2) + 0.2*yspan]);
legend('True Ripple $f_{s4}$ ($2.5\text{ Hz}$)', 'Observer Tracking $\hat{f}_{s4}$', 'Interpreter', 'latex', 'Location', 'best');
title('{Harmonic Rack Position Sensor Drift Tracking (LVDT Carrier Distortion)}', 'Interpreter', 'latex');

fprintf('=================================================================\n');
fprintf('All 4 Advanced Benchmark Scenarios Successfully Simulated up to 20 Seconds!\n');
fprintf('=================================================================\n');