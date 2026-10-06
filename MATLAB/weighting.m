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

theta0_weighting = deg2rad(40);

Dphi = k*d*sin(theta0_weighting);

psi = pi*d*sin(theta)/lambda - Dphi/2;

AF = sin(N.*psi) ./ (N*sin(psi));
singular = abs(sin(psi)) < 1e-12;
AF(singular) = 1;
AF_dB = 20*log10(abs(AF));

% Weighting Cartesian figure
fig8 = figure;
hold on;
grid on;

% Weighting Polar figure
fig9 = figure;
pax4 = polaraxes;
hold(pax4, 'on');

% Weighting
Dphi = k*d*sin(theta0_weighting);

n = 0:N-1;

wh = 0.54 - 0.46*cos(2*pi*n/(N-1));

wt = taylorwin(N, 4, -30);

AF_Hamming_weighted = zeros(size(theta));
AF_Taylor_weighted = zeros(size(theta));

for i = 1:N

    AF_Hamming_weighted = AF_Hamming_weighted + ...
        wh(i)*exp(1j*(i-1)*(k*d*sin(theta) - Dphi));

end

for i = 1:N

    AF_Taylor_weighted = AF_Taylor_weighted + ...
        wt(i)*exp(1j*(i-1)*(k*d*sin(theta) - Dphi));

end

AF_Hamming_weighted = abs(AF_Hamming_weighted);
AF_Hamming_weighted = AF_Hamming_weighted/max(AF_Hamming_weighted);
AF_Hamming_weighted_dB = 20*log10(AF_Hamming_weighted);

AF_Taylor_weighted = abs(AF_Taylor_weighted);
AF_Taylor_weighted = AF_Taylor_weighted/max(AF_Taylor_weighted);
AF_Taylor_weighted_dB = 20*log10(AF_Taylor_weighted);

AF_dB_polar = AF_dB;
AF_dB_polar(AF_dB_polar < -40) = -40;
rho_uniform = AF_dB_polar + 40;

[~, idx_peak_Hamming] = max(AF_Hamming_weighted_dB);
[~, idx_peak_Taylor] = max(AF_Taylor_weighted_dB);
[~, idx_peak_Uniform] = max(AF_dB);

idx_left = find(AF_Hamming_weighted_dB(1:idx_peak_Hamming) <= -3, 1, 'last');
idx_right_rel = find(AF_Hamming_weighted_dB(idx_peak_Hamming:end) <= -3, 1, 'first');
idx_right = idx_peak_Hamming + idx_right_rel - 1;
HPBW_Hamming = theta_deg(idx_right) - theta_deg(idx_left);

idx_left = find(AF_Taylor_weighted_dB(1:idx_peak_Taylor) <= -3, 1, 'last');
idx_right_rel = find(AF_Taylor_weighted_dB(idx_peak_Taylor:end) <= -3, 1, 'first');
idx_right = idx_peak_Taylor + idx_right_rel - 1;
HPBW_Taylor = theta_deg(idx_right) - theta_deg(idx_left);

idx_left = find(AF_dB(1:idx_peak_Uniform) <= -3, 1, 'last');
idx_right_rel = find(AF_dB(idx_peak_Uniform:end) <= -3, 1, 'first');
idx_right = idx_peak_Uniform + idx_right_rel - 1;
HPBW_Uniform = theta_deg(idx_right) - theta_deg(idx_left);

visible = theta_deg >= -90 & theta_deg <= 90;

[pks, ~] = findpeaks(AF_Hamming_weighted_dB(visible));
[~, main_idx] = max(pks);
pks(main_idx) = [];
PSLL_Hamming = max(pks);

[pks, ~] = findpeaks(AF_Taylor_weighted_dB(visible));
[~, main_idx] = max(pks);
pks(main_idx) = [];
PSLL_Taylor = max(pks);

[pks, ~] = findpeaks(AF_dB(visible));
[~, main_idx] = max(pks);
pks(main_idx) = [];
PSLL_Uniform = max(pks);

fprintf('--- Weighting Metrics ---\n');

fprintf('Uniform: HPBW = %.2f deg, PSLL = %.2f dB\n', ...
    HPBW_Uniform, PSLL_Uniform);

fprintf('Hamming: HPBW = %.2f deg, PSLL = %.2f dB\n', ...
    HPBW_Hamming, PSLL_Hamming);

fprintf('Taylor:  HPBW = %.2f deg, PSLL = %.2f dB\n', ...
    HPBW_Taylor, PSLL_Taylor);

%Cartesian
figure(fig8);

plot(rad2deg(theta), AF_Hamming_weighted_dB, ...
    'LineWidth', 1.8, ...
    'DisplayName', sprintf('With Hamming weighting'));

plot(rad2deg(theta), AF_Taylor_weighted_dB, ...
    'LineWidth', 1.8, ...
    'DisplayName', sprintf('With Taylor weighting'));

plot(rad2deg(theta), AF_dB, ...
    'LineWidth', 1.8, ...
    'DisplayName', sprintf('Uniform weighting'));

% Polar
Total_Hamming_weighted_dB_polar = AF_Hamming_weighted_dB;
Total_Hamming_weighted_dB_polar(Total_Hamming_weighted_dB_polar < -40) = -40;
rho4 = Total_Hamming_weighted_dB_polar + 40;

Total_Taylor_weighted_dB_polar = AF_Taylor_weighted_dB;
Total_Taylor_weighted_dB_polar(Total_Taylor_weighted_dB_polar < -40) = -40;
rho5 = Total_Taylor_weighted_dB_polar + 40;

polarplot(pax4, theta, rho4, ...
    'LineWidth', 1.8, ...
    'DisplayName', sprintf('With Hamming weighting'));

polarplot(pax4, theta, rho5, ...
    'LineWidth', 1.8, ...
    'DisplayName', sprintf('With Taylor weighting'));

polarplot(pax4, theta, rho_uniform, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Uniform weighting');

% Weighting settings
figure(fig8);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title('Weighting algorithms at center angle of 40° (AF only)');
xlim([-90 90]);
ylim([-40 0]);
legend('show');

% Weighting Polar settings
figure(fig9);
pax4.RLim = [0 40];
pax4.RTick = [0 10 20 30 40];
pax4.RTickLabel = {'-40','-30','-20','-10','0'};
title(pax4, 'Weighting algorithms at center angle of 40° (AF only)');
legend(pax4, 'show');