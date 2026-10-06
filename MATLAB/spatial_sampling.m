clear;
clc;
close all;

% Parameters
N = 8;
f0 = 1e10;
c = 3*10e8;

lambda = c / f0;

% Angle sweep
theta_deg = -180:0.01:180;
theta = deg2rad(theta_deg);

% Spatial Sampling figure
fig4 = figure;
hold on;
grid on;

% Grating Lobe Spatial Sampling Polar figure
fig5 = figure;
pax2 = polaraxes;
hold(pax2, 'on');

d_values_grating_lobe = [0.5*lambda, 0.75*lambda, lambda];
theta0_spatial = deg2rad(40);

fprintf('--- Grating Lobe Positions --\n')

for d = d_values_grating_lobe

    Dphi = (2*pi*d*sin(theta0_spatial)) / lambda;

    psi = pi*d*sin(theta)/lambda - Dphi/2;

    AF = sin(N.*psi) ./ (N*sin(psi));

    singular = abs(sin(psi)) < 1e-12;
    AF(singular) = 1;

    AF_dB = 20*log10(abs(AF));

    % Cartesian
    figure(fig4);

    plot(rad2deg(theta), AF_dB, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('d / lambda = %.2f', d/lambda));

    % Polar
    Total_dB_polar = AF_dB;
    Total_dB_polar(Total_dB_polar < -40) = -40;

    rho2 = Total_dB_polar + 40;

    polarplot(pax2, theta, rho2, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('d / lambda = %.2f', d/lambda));

    m_values = -3:3;

    for m = m_values

        if m == 0
            continue
        end

        argument = sin(theta0_spatial) + m*lambda/d;

        if abs(argument) <= 1

            theta_m = asin(argument);

            fprintf('Grating lobe for d/lambda = %.2f for m = %d at %.1f°\n', ...
                d/lambda, m, rad2deg(theta_m));

        end

    end

end

% Spatial Sampling settings
figure(fig4);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title('Spatial Sampling with center angle of 40° for d / lambda > 0.5 (AF only)');
xlim([-90 90]);
ylim([-40 0]);
legend('show');

% Spatial Sampling Polar settings
figure(fig5);
pax2.RLim = [0 40];
pax2.RTick = [0 10 20 30 40];
pax2.RTickLabel = {'-40','-30','-20','-10','0'};
title(pax2, 'Spatial Sampling with center angle of 40° for d / lambda > 0.5 (AF only)');
legend(pax2, 'show');

% Aperture Spatial Sampling figure
fig6 = figure;
hold on;
grid on;

% Aperture Sampling Polar figure
fig7 = figure;
pax3 = polaraxes;
hold(pax3, 'on');

d_values_aperture = [0.1*lambda, 0.25*lambda, 0.5*lambda];

for d = d_values_aperture

    Dphi = (2*pi*d*sin(theta0_spatial)) / lambda;

    psi = pi*d*sin(theta)/lambda - Dphi/2;

    AF = sin(N.*psi) ./ (N*sin(psi));

    singular = abs(sin(psi)) < 1e-12;
    AF(singular) = 1;

    AF_dB = 20*log10(abs(AF));

    % Cartesian
    figure(fig6);

    plot(rad2deg(theta), AF_dB, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('d / lambda = %.2f', d/lambda));

    % Polar
    Total_dB_polar = AF_dB;
    Total_dB_polar(Total_dB_polar < -40) = -40;

    rho3 = Total_dB_polar + 40;

    polarplot(pax3, theta, rho3, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('d / lambda = %.2f', d/lambda));

end

% Aperture Sampling settings
figure(fig6);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title('Spatial Sampling with center angle of 40° for d / lambda < 0.5 (AF only)');
xlim([-90 90]);
ylim([-40 0]);
legend('show');

% Aperture Sampling Polar settings
figure(fig7);
pax3.RLim = [0 40];
pax3.RTick = [0 10 20 30 40];
pax3.RTickLabel = {'-40','-30','-20','-10','0'};
title(pax3, 'Spatial Sampling with center angle of 40° for d / lambda < 0.5 (AF only)');
legend(pax3, 'show');
