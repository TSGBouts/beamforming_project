clear;
clc;
close all;

ideal_data = load('photonic_ttd_30deg_ideal.mat');
nonideal_data = load('photonic_ttd_30deg_non-ideal.mat');
thermal_data = load('photonic_ttd_30deg_thermal.mat');

f = ideal_data.f_export;

H_ideal = ideal_data.Hrel_export;
H_nonideal = nonideal_data.Hrel_export;
H_thermal = thermal_data.Hrel_export;

Habs_ideal = ideal_data.H_export;
Habs_nonideal = nonideal_data.H_export;
Habs_thermal = thermal_data.H_export;

N = 8;
c = 3e8;
f0 = 10e9;
theta0 = 30;

lambda0 = c/f0;
d = lambda0/2;
k0 = 2*pi*f0/c;

n = 0:N-1;

theta_deg = -90:0.01:90;
theta = deg2rad(theta_deg);

Nf = length(f);

dtau = d*sind(theta0)/c;
tau = n*dtau;

error_PS = zeros(Nf,1);
error_TTD = zeros(Nf,1);
error_ideal = zeros(Nf,1);
error_nonideal = zeros(Nf,1);
error_thermal = zeros(Nf,1);

theta_peak_TTD = zeros(Nf,1);
theta_peak_ideal = zeros(Nf,1);
theta_peak_nonideal = zeros(Nf,1);
theta_peak_thermal = zeros(Nf,1);

w_PS = exp(-1j*n*k0*d*sind(theta0));

for kf = 1:Nf

    fk = f(kf);
    k = 2*pi*fk/c;

    steering = exp(1j*(n.')*(k*d*sin(theta)));

    AF_PS = abs(w_PS*steering);
    AF_PS = AF_PS/max(AF_PS);

    w_TTD = exp(-1j*2*pi*fk*tau);

    AF_TTD = abs(w_TTD*steering);
    AF_TTD = AF_TTD/max(AF_TTD);

    AF_ideal = abs(H_ideal(kf,:)*steering);
    AF_ideal = AF_ideal/max(AF_ideal);

    AF_nonideal = abs(H_nonideal(kf,:)*steering);
    AF_nonideal = AF_nonideal/max(AF_nonideal);

    AF_thermal = abs(H_thermal(kf,:)*steering);
    AF_thermal = AF_thermal/max(AF_thermal);

    [~,idx_PS] = max(AF_PS);
    [~,idx_TTD] = max(AF_TTD);
    [~,idx_ideal] = max(AF_ideal);
    [~,idx_nonideal] = max(AF_nonideal);
    [~,idx_thermal] = max(AF_thermal);

    error_PS(kf) = theta_deg(idx_PS) - theta0;
    error_TTD(kf) = theta_deg(idx_TTD) - theta0;
    error_ideal(kf) = theta_deg(idx_ideal) - theta0;
    error_nonideal(kf) = theta_deg(idx_nonideal) - theta0;
    error_thermal(kf) = theta_deg(idx_thermal) - theta0;

    theta_peak_TTD(kf) = theta_deg(idx_TTD);
    theta_peak_ideal(kf) = theta_deg(idx_ideal);
    theta_peak_nonideal(kf) = theta_deg(idx_nonideal);
    theta_peak_thermal(kf) = theta_deg(idx_thermal);

end

figure;
hold on;
grid on;

plot(f/1e9,error_PS,'LineWidth',1.5);
plot(f/1e9,error_ideal,':','LineWidth',1.7);
plot(f/1e9,error_nonideal,'-.','LineWidth',1.5);
plot(f/1e9,error_thermal,'LineWidth',1.7);

yline(0,'k--');

xlabel('Frequency (GHz)');
ylabel('Pointing Error (degrees)');
title('Pointing Error vs Frequency');

legend( ...
    'Phase Shifter', ...
    'Ideal Photonic TTD', ...
    'Nonideal Photonic TTD', ...
    'Thermal + Noise Photonic TTD', ...
    'Location','best');

f_plot = [7.5e9 10e9 12.5e9];

for q = 1:length(f_plot)

    [~,kf] = min(abs(f-f_plot(q)));

    fk = f(kf);
    k = 2*pi*fk/c;

    steering = exp(1j*(n.')*(k*d*sin(theta)));

    AF_PS = abs(w_PS*steering);
    AF_PS = AF_PS/max(AF_PS);

    w_TTD = exp(-1j*2*pi*fk*tau);

    AF_TTD = abs(w_TTD*steering);
    AF_TTD = AF_TTD/max(AF_TTD);

    AF_ideal = abs(H_ideal(kf,:)*steering);
    AF_ideal = AF_ideal/max(AF_ideal);

    AF_nonideal = abs(H_nonideal(kf,:)*steering);
    AF_nonideal = AF_nonideal/max(AF_nonideal);

    AF_thermal = abs(H_thermal(kf,:)*steering);
    AF_thermal = AF_thermal/max(AF_thermal);

    figure;
    hold on;
    grid on;

    plot(theta_deg,20*log10(max(AF_PS,1e-12)), ...
        'LineWidth',1.5);

    plot(theta_deg,20*log10(max(AF_ideal,1e-12)), ...
        ':','LineWidth',1.7);

    plot(theta_deg,20*log10(max(AF_nonideal,1e-12)), ...
        '-.','LineWidth',1.5);

    plot(theta_deg,20*log10(max(AF_thermal,1e-12)), ...
        'LineWidth',1.7);

    xlabel('Angle (degrees)');
    ylabel('Normalized Array Factor (dB)');

    title(['Beam Pattern at ', ...
        num2str(fk/1e9,'%.1f'),' GHz']);

    legend( ...
        'Phase Shifter', ...
        'Ideal Photonic TTD', ...
        'Nonideal Photonic TTD', ...
        'Thermal + Noise Photonic TTD', ...
        'Location','best');

    xlim([-90 90]);
    ylim([-40 0]);

end

nonideal_deviation = ...
    abs(theta_peak_nonideal - theta_peak_ideal);

thermal_deviation = ...
    abs(theta_peak_thermal - theta_peak_ideal);

figure;
hold on;
grid on;

plot(f/1e9,nonideal_deviation, ...
    'LineWidth',1.5);

plot(f/1e9,thermal_deviation, ...
    'LineWidth',1.5);

xlabel('Frequency (GHz)');
ylabel('Deviation from Ideal Photonic TTD (degrees)');
title('Photonic TTD Degradation vs Frequency');

legend( ...
    'Nonideal Photonic TTD', ...
    'Thermal + Noise Photonic TTD', ...
    'Location','best');

[~,kf0] = min(abs(f-f0));

relative_nonideal_dB = ...
    20*log10(abs(H_nonideal(kf0,:)));

relative_thermal_dB = ...
    20*log10(abs(H_thermal(kf0,:)));

figure;
hold on;
grid on;

plot(1:N,relative_nonideal_dB, ...
    'o-','LineWidth',1.5);

plot(1:N,relative_thermal_dB, ...
    's-','LineWidth',1.5);

yline(0,'k--');

xlabel('Channel');
ylabel('Relative Magnitude (dB)');
title('Channel Imbalance at 10 GHz');

legend( ...
    'Nonideal Photonic TTD', ...
    'Thermal + Noise Photonic TTD', ...
    'Location','best');

phase_nonideal = ...
    rad2deg(angle(H_nonideal(kf0,:) ./ H_ideal(kf0,:)));

phase_thermal = ...
    rad2deg(angle(H_thermal(kf0,:) ./ H_ideal(kf0,:)));

figure;
hold on;
grid on;

plot(1:N,phase_nonideal, ...
    'o-','LineWidth',1.5);

plot(1:N,phase_thermal, ...
    's-','LineWidth',1.5);

yline(0,'k--');

xlabel('Channel');
ylabel('Phase Error (degrees)');
title('Channel Phase Error at 10 GHz');

legend( ...
    'Nonideal Photonic TTD', ...
    'Thermal + Noise Photonic TTD', ...
    'Location','best');

if isfield(nonideal_data,'tau_err') && ...
        isfield(thermal_data,'tau_err')

    figure;
    hold on;
    grid on;

    plot(1:N,nonideal_data.tau_err*1e12, ...
        'o-','LineWidth',1.5);

    plot(1:N,thermal_data.tau_err*1e12, ...
        's-','LineWidth',1.5);

    yline(0,'k--');

    xlabel('Channel');
    ylabel('Delay Error (ps)');
    title('Photonic TTD Delay Error');

    legend( ...
        'Nonideal Photonic TTD', ...
        'Thermal + Noise Photonic TTD', ...
        'Location','best');

end

absolute_nonideal_dB = ...
    -20*log10(abs( ...
    Habs_nonideal(kf0,:) ./ Habs_ideal(kf0,:)));

absolute_thermal_dB = ...
    -20*log10(abs( ...
    Habs_thermal(kf0,:) ./ Habs_ideal(kf0,:)));

figure;
hold on;
grid on;

plot(1:N,absolute_nonideal_dB, ...
    'o-','LineWidth',1.5);

plot(1:N,absolute_thermal_dB, ...
    's-','LineWidth',1.5);

xlabel('Channel');
ylabel('Degradation Relative to Ideal (dB)');
title('Absolute Photonic Degradation at 10 GHz');

legend( ...
    'Nonideal Photonic TTD', ...
    'Thermal + Noise Photonic TTD', ...
    'Location','best');

max_error_PS = max(abs(error_PS));
max_error_TTD = max(abs(error_TTD));
max_error_ideal = max(abs(error_ideal));
max_error_nonideal = max(abs(error_nonideal));
max_error_thermal = max(abs(error_thermal));

max_nonideal_added = max(nonideal_deviation);
max_thermal_added = max(thermal_deviation);

[~,idx_worst_thermal] = max(thermal_deviation);

imbalance_nonideal = ...
    max(relative_nonideal_dB) - min(relative_nonideal_dB);

imbalance_thermal = ...
    max(relative_thermal_dB) - min(relative_thermal_dB);

fprintf('\n============================================\n');
fprintf('      WIDEBAND BEAMFORMING COMPARISON\n');
fprintf('============================================\n');

fprintf('\nTarget steering angle: %.2f deg\n',theta0);

fprintf('\nMaximum pointing error:\n');
fprintf('Phase Shifter:          %.3f deg\n',max_error_PS);
fprintf('Ideal Photonic TTD:     %.3f deg\n',max_error_ideal);
fprintf('Nonideal Photonic TTD:  %.3f deg\n',max_error_nonideal);
fprintf('Thermal + Noise TTD:    %.3f deg\n',max_error_thermal);

fprintf('\nDeviation from ideal photonic TTD:\n');
fprintf('Nonideal contribution:  %.3f deg\n',max_nonideal_added);
fprintf('Thermal/noise contribution: %.3f deg\n',max_thermal_added);

fprintf('\nWorst thermal/noise frequency: %.3f GHz\n', ...
    f(idx_worst_thermal)/1e9);

fprintf('\nChannel imbalance at 10 GHz:\n');
fprintf('Nonideal:       %.3f dB\n',imbalance_nonideal);
fprintf('Thermal + Noise: %.3f dB\n',imbalance_thermal);

if isfield(nonideal_data,'tau_MAE')

    fprintf('\nNonideal delay error:\n');
    fprintf('MAE:        %.4f ps\n', ...
        nonideal_data.tau_MAE*1e12);
    fprintf('RMS:        %.4f ps\n', ...
        nonideal_data.tau_RMS*1e12);
    fprintf('Worst case: %.4f ps\n', ...
        nonideal_data.tau_MAX*1e12);

end

if isfield(thermal_data,'tau_MAE')

    fprintf('\nThermal + noise delay error:\n');
    fprintf('MAE:        %.4f ps\n', ...
        thermal_data.tau_MAE*1e12);
    fprintf('RMS:        %.4f ps\n', ...
        thermal_data.tau_RMS*1e12);
    fprintf('Worst case: %.4f ps\n', ...
        thermal_data.tau_MAX*1e12);

end

fprintf('\n============================================\n');