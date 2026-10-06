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

% Cartesian figure
fig1 = figure;
hold on;
grid on;

% Polar figure
fig2 = figure;
pax = polaraxes;
hold(pax, 'on');

fprintf('---Phase Shifter Implementation---\n')

for theta0 = theta0_values

    Dphi = k*d*sin(theta0);

    psi = pi*d*sin(theta)/lambda - Dphi/2;

    AF = sin(N.*psi) ./ (N*sin(psi));

    singular = abs(sin(psi)) < 1e-12;
    AF(singular) = 1;

    AF_dB = 20*log10(abs(AF));

    EF = cos(theta);

    Total = abs(AF).*abs(EF);
    Total = Total./max(Total);

    Total_dB = 20*log10(abs(Total));

    fprintf('Phase difference of adjecent antennas for theta0 = %d: %.1f°\n', ...
        rad2deg(theta0), rad2deg(Dphi));

    % Cartesian
    figure(fig1);

    plot(theta_deg, Total_dB, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('\\theta_0 = %d°', rad2deg(theta0)));

    % Polar
    Total_dB_polar = Total_dB;
    Total_dB_polar(Total_dB_polar < -40) = -40;

    rho = Total_dB_polar + 40;

    polarplot(pax, theta, rho, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('\\theta_0 = %d°', rad2deg(theta0)));

    HPBW_AF = 0.886*lambda/(N*d*cos(theta0));
    HPBW_deg_AF = rad2deg(HPBW_AF);

    AFfun = @(th) abs( ...
        sin(N*(pi*d*sin(th)/lambda - Dphi/2)) ./ ...
        (N*sin(pi*d*sin(th)/lambda - Dphi/2)) ...
        ) - 1/sqrt(2);

    theta_right = fzero(AFfun, [theta0 + 0.0001, deg2rad(60)]);
    theta_left = fzero(AFfun, [deg2rad(-60), theta0 - 0.0001]);

    HPBW_exact_AF = rad2deg(theta_right - theta_left);

    fprintf('HPBW_AF approximation: %.4f --- ', HPBW_deg_AF);
    fprintf('HPBW_AF exact: %.4f°\n', HPBW_exact_AF);

end

fprintf('\n')

% Cartesian settings
figure(fig1);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title('Normalized Total Pattern (Array Factor & Element Factor)');
xlim([-90 90]);
ylim([-40 0]);
legend('show');

% Polar settings
figure(fig2);
pax.RLim = [0 40];
pax.RTick = [0 10 20 30 40];
pax.RTickLabel = {'-40','-30','-20','-10','0'};
title(pax, 'Normalized Total Pattern (Array Factor & Element Factor)');
legend(pax, 'show');