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

theta0_deg_values = [0, 20, 40];
theta0_values = deg2rad(theta0_deg_values);

theta0_wideband_deg_values = [10, 30, 50];
theta0_wideband_values = deg2rad(theta0_wideband_deg_values);

f_wideband_values = 9e9:0.25e9:11e9;
DeltaTheta_PS_values = zeros(size(f_wideband_values));

visible = theta_deg >= -90 & theta_deg <= 90;

% True Time Delay
fig11 = figure;
hold on;
grid on;

fig12 = figure;
pax6 = polaraxes;
hold(pax6, 'on');

fprintf('--- True Time Delay Implementation ---\n')

for theta0 = theta0_values

    t_delay = d * sin(theta0) / c;

    AF_TTD = sin((N*pi*d*(sin(theta) - sin(theta0))/lambda)) ./ ...
        (N.*sin((pi*d*(sin(theta) - sin(theta0))/lambda)));

    singular = abs(sin(theta) - sin(theta0)) < 1e-12;
    AF_TTD(singular) = 1;

    AF_TTD = AF_TTD./max(AF_TTD);
    AF_TTD_dB = 20*log10(abs(AF_TTD));

    figure(fig11);

    plot(rad2deg(theta), AF_TTD_dB, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('theta0 = %d°', rad2deg(theta0)));

    AF_TTD_dB_polar = AF_TTD_dB;
    AF_TTD_dB_polar(AF_TTD_dB_polar < -40) = -40;

    rho6 = AF_TTD_dB_polar + 40;

    polarplot(pax6, theta, rho6, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('theta0 = %d°', rad2deg(theta0)));

    fprintf('Required time delay between elemens for theta0 = %d deg is: %.2f picoseconds, while maximum is: %.2f picoseconds\n', ...
        rad2deg(theta0), t_delay * 1e12, t_delay * 1e12 * (N-1));

end

fig13 = figure;
hold on;
grid on;

fig14 = figure;
pax7 = polaraxes;
hold(pax7, 'on');

theta0 = deg2rad(30);
DeltaTheta_TTD_values = zeros(size(f_wideband_values));

i = 1;

for f = f_wideband_values

    lambda_f = c / f;

    AF_TTD_f = sin((N*pi*d*(sin(theta) - sin(theta0))/lambda_f)) ./ ...
        (N.*sin((pi*d*(sin(theta) - sin(theta0))/lambda_f)));

    singular = abs(sin(theta) - sin(theta0)) < 1e-12;
    AF_TTD_f(singular) = 1;

    AF_TTD_f = AF_TTD_f./max(AF_TTD_f);
    AF_TTD_f_dB = 20*log10(abs(AF_TTD_f));

    theta_visible = theta(visible);
    AF_TTD_f_dB_visible = AF_TTD_f_dB(visible);

    [~, idx_peak_TTD] = max(AF_TTD_f_dB_visible);

    DeltaTheta_TTD_values(i) = ...
        rad2deg(theta_visible(idx_peak_TTD) - theta0);

    % Phase Shifter Angle Error 
    psi_PS = pi*d*(sin(theta)/lambda_f - sin(theta0)/lambda);

    AF_PS_f = sin(N*psi_PS) ./ (N*sin(psi_PS));
    singular = abs(sin(psi_PS)) < 1e-12;
    AF_PS_f(singular) = 1;

    AF_PS_f = abs(AF_PS_f);
    AF_PS_f = AF_PS_f ./ max(AF_PS_f);
    AF_PS_f_dB = 20*log10(AF_PS_f);

    theta_visible = theta(visible);
    AF_PS_f_dB_visible = AF_PS_f_dB(visible);

    [~, idx_peak_PS] = max(AF_PS_f_dB_visible);

    theta_peak_PS = theta_visible(idx_peak_PS);

    DeltaTheta_PS_values(i) = ...
        rad2deg(theta_peak_PS - theta0);

    % Plotting
    if (f == 9e9 || f == 1e10 || f == 11e9)

    figure(fig13);

    plot(rad2deg(theta), AF_TTD_f_dB, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('f / f0 = %.1f', f / f0));

    AF_TTD_f_dB_polar = AF_TTD_f_dB;
    AF_TTD_f_dB_polar(AF_TTD_f_dB_polar < -40) = -40;

    rho7 = AF_TTD_f_dB_polar + 40;

    polarplot(pax7, theta, rho7, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('f / f0 = %.1f', f / f0));

    end

    i = i + 1;

end

fprintf('\n--- Phase Shifter vs True Time Delay (Steering Angle = 30°) ---\n');

fprintf('f/f0 | PS Error | TTD Error | Difference (deg)\n');

fprintf('----------------------------------------------\n');

for i = 1:length(f_wideband_values)

    fprintf('%.2f | %5.2f    | %5.2f     | %5.2f\n', ...
        f_wideband_values(i)/f0, ...
        DeltaTheta_PS_values(i), ...
        DeltaTheta_TTD_values(i), ...
        abs(DeltaTheta_PS_values(i) - DeltaTheta_TTD_values(i)));

end

% TTD settings
figure(fig11);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title('True Time Delay Implementation at center frequency');
xlim([-90 90]);
ylim([-40 0]);
legend('show');

% TTD Polar settings
figure(fig12);
pax6.RLim = [0 40];
pax6.RTick = [0 10 20 30 40];
pax6.RTickLabel = {'-40','-30','-20','-10','0'};
title(pax6, 'True Time Delay Implementation at center frequency');
legend(pax6, 'show');

figure(fig13);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title('True Time Delay Implementation at different frequencies for a target angle of 30°');
xlim([-90 90]);
ylim([-40 0]);
legend('show');

% TTD Polar settings
figure(fig14);
pax7.RLim = [0 40];
pax7.RTick = [0 10 20 30 40];
pax7.RTickLabel = {'-40','-30','-20','-10','0'};
title(pax7, 'True Time Delay Implementation at different frequencies for a target angle of 30°');
legend(pax7, 'show');