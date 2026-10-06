clear;
clc;
close all;

% Parameters
N = 8;
f0 = 1e10;
c = 3*10e8;

lambda = c / f0;
d = lambda / 2;
k = 2*pi/lambda;

% Angle sweep
theta_deg = -180:0.01:180;
theta = deg2rad(theta_deg);

% Beam Squint figure
fig3 = figure;
hold on;
grid on;

theta0_sweep_deg = 0:0.1:60;
theta0_sweep = deg2rad(theta0_sweep_deg);

theta0_deg_values = [0, 20, 40];
theta0_values = deg2rad(theta0_deg_values);

f_values = [7e9, 8e9, 9e9, 1e10, 11e9, 12e9, 13e9];

for f = f_values

    a = f0/f;

    argument = a*sin(theta0_sweep);

    Dtheta = nan(size(theta0_sweep));

    valid = abs(argument) <= 1;

    Dtheta(valid) = asin(argument(valid)) - theta0_sweep(valid);

    figure(fig3);

    plot(rad2deg(theta0_sweep), rad2deg(Dtheta), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('f / f0 = %.1f', 1/a));

end

% Beam Squint settings
figure(fig3);
grid on;
xlabel('Pointed Beam Angle at Center Frequency');
ylabel('Beam Angle Deviation');
title('Beam Squint');
xlim([0 60]);
ylim([-40 40]);
legend('show');

for theta0 = theta0_values

    Dphi = k*d*sin(theta0);

    psi = pi*d*sin(theta)/lambda - Dphi/2;

    AF = sin(N.*psi) ./ (N*sin(psi));

    singular = abs(sin(psi)) < 1e-12;
    AF(singular) = 1;

    AF_dB = 20*log10(abs(AF));

end

% Fixed Phase Shifter Beam Squint figure
fig10 = figure;
hold on;
grid on;

theta0_wideband_deg_values = [10, 30, 50];
theta0_wideband_values = deg2rad(theta0_wideband_deg_values);

f_wideband_values = 9e9:0.25e9:11e9;
DeltaTheta_PS_values = zeros(size(f_wideband_values));

for theta0 = theta0_wideband_values

    % Phase shift designed at center frequency
    Dphi = 2*pi*f0*d*sin(theta0)/c;

    theta_peak_values = zeros(size(f_wideband_values));
    theta_analytical_values = zeros(size(f_wideband_values));

    i = 1;

    for f = f_wideband_values

        kf = 2*pi*f/c;

        psi = (kf*d*sin(theta) - Dphi)/2;

        AF = sin(N.*psi) ./ (N*sin(psi));

        singular = abs(sin(psi)) < 1e-12;
        AF(singular) = 1;

        AF_dB = 20*log10(abs(AF));

        visible = theta_deg >= -90 & theta_deg <= 90;

        AF_visible = AF_dB(visible);
        theta_visible = theta_deg(visible);

        [~, idx_peak] = max(AF_visible);

        theta_peak_values(i) = theta_visible(idx_peak);

        theta_analytical_values(i) = ...
            rad2deg(asin((f0/f)*sin(theta0)));

        if theta0 == theta0_wideband_values(2)
            DeltaTheta_PS_values(i) = ...
            theta_peak_values(i) - rad2deg(theta0);
        end

        i = i + 1;

    end

    figure(fig10);

    plot(f_wideband_values/1e9, theta_peak_values, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('Simulation theta0 = %.0f°', ...
        rad2deg(theta0)));

    plot(f_wideband_values/1e9, theta_analytical_values, '--', ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('Analytical theta0 = %.0f°', ...
        rad2deg(theta0)));

end

% Fixed Phase Shifter Wideband settings
figure(fig10);
grid on;
xlabel('Frequency (GHz)');
ylabel('Beam Direction (degrees)');
title('Fixed Phase Shifter Beam Squint');
xlim([9 11]);
legend('show');

