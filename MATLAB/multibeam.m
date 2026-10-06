clear;
clc;
close all;

% Parameters
N = 8;
f0 = 10e9;
c = 3e8;

lambda = c/f0;
d = lambda/2;
k = 2*pi/lambda;

theta0_values = [-25 35];

% Angle sweep
theta_deg = -90:0.1:90;
theta = deg2rad(theta_deg);

n = 0:N-1;
w = zeros(1,N);

% Cartesian figures
fig2 = figure;
hold on;
grid on;

fig1 = figure;
hold on;
grid on;

fprintf('---Multibeam Implementation for 2 beams---\n')

for theta0 = theta0_values
    w = w + exp(-1j*n*k*d*sind(theta0));
end

AF = zeros(size(theta));

for m = 1:N
    AF = AF + w(m).*exp(1j*n(m)*k*d*sin(theta));
end

visible = theta_deg >= -90 & theta_deg <= 90;

AF = abs(AF);
AF = AF./max(AF);
AF_dB = 20*log10(AF);

theta_visible = theta_deg(visible);
AF_visible = AF(visible);

[pks, locs] = findpeaks(AF_visible);

peak_angles = theta_visible(locs);

[~, idx1] = min(abs(peak_angles - theta0_values(1)));

[~, idx2] = min(abs(peak_angles - theta0_values(2)));

beam1_angle = peak_angles(idx1);
beam2_angle = peak_angles(idx2);

beam1_level = 20*log10(pks(idx1));
beam2_level = 20*log10(pks(idx2));

fprintf('Beam 1: %.2f deg, %.2f dB\n', ...
    beam1_angle, beam1_level);

fprintf('Beam 2: %.2f deg, %.2f dB\n', ...
    beam2_angle, beam2_level);

fprintf('Difference of strength between two beams: %.d dB\n', abs(beam2_level-beam1_level));

pks([idx1 idx2]) = [];
locs([idx1 idx2]) = [];

[max_sidelobe, idx] = max(pks);

PSLL = 20*log10(max_sidelobe);

fprintf('Largest remaining local peak: %.2f dB\n', PSLL);

figure(fig2);
plot(theta_deg,AF_dB,'LineWidth',1.8);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title(sprintf( ...
    'Multibeam: \\theta_1 = %g°, \\theta_2 = %g°', ...
    theta0_values(1),theta0_values(2)));
xlim([-90 90]);
ylim([-40 0]);

fprintf('\n---Beam Squint Error for varying frequencies---')

f_values = [8e9, f0, 12e9];

for f = f_values

    n = 0:N-1;

    lambda0 = c/f;
    k0 = 2*pi/lambda0;

    AF_f = zeros(size(theta));

    for m = 1:N
        AF_f = AF_f + w(m).*exp(1j*n(m)*k0*d*sin(theta));
    end

    AF_f = abs(AF_f);
    AF_f = AF_f./max(AF_f);
    AF_dB_f = 20*log10(AF_f);

    [pks_f, locs_f] = findpeaks(AF_f);

    peak_angles_f = theta_deg(locs_f);

    % Peak closest to expected squinted position
    expected1 = asind((f0/f)*sind(theta0_values(1)));
    expected2 = asind((f0/f)*sind(theta0_values(2)));

    [~, idx1_f] = min(abs(peak_angles_f - expected1));
    [~, idx2_f] = min(abs(peak_angles_f - expected2));

    beam1_angle_f = peak_angles_f(idx1_f);
    beam2_angle_f = peak_angles_f(idx2_f);

    error1 = beam1_angle_f - theta0_values(1);
    error2 = beam2_angle_f - theta0_values(2);

    fprintf('\nf/f0 = %.1f\n',f/f0);
    fprintf('Beam 1: %.2f deg, error = %.2f deg\n', ...
        beam1_angle_f,error1);
    fprintf('Beam 2: %.2f deg, error = %.2f deg\n', ...
        beam2_angle_f,error2);

    figure(fig1);
    plot(rad2deg(theta), AF_dB_f, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('f/f0 = %.1f', f/f0));
end

% Plot
figure(fig1)
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title(sprintf( ...
    'Beam Sqiunt: \\theta_1 = %g°, \\theta_2 = %g°', ...
    theta0_values(1),theta0_values(2)));
xlim([-90 90]);
ylim([-40 0]);
legend('show');
