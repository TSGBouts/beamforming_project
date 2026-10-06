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

n = 0:N-1;

theta0 = 50;

% Cartesian figures
fig1 = figure;
hold on;
grid on;

fig2 = figure;
hold on;
grid on;

% Histograms
fig3 = figure;
hold on;
grid on;

fig4 = figure;
hold on;
grid on;

% Angle sweep
theta_deg = -90:0.01:90;
theta = deg2rad(theta_deg);

visible = theta_deg >= -90 & theta_deg <= 90;
theta_visible = theta(visible);

phase_ideal = -rad2deg(n*k*d*sind(theta0));
phase_ideal = mod(phase_ideal,360);

w_ideal = exp(1j*deg2rad(phase_ideal));
AF_ideal = zeros(size(theta));

for m = 1:N
    AF_ideal = AF_ideal + ...
        w_ideal(m).*exp(1j*n(m)*k*d*sin(theta));
end

AF_ideal = abs(AF_ideal);
AF_ideal = AF_ideal/max(AF_ideal);
AF_ideal_dB = 20*log10(AF_ideal);

% 2, 4 & 6-bit Phase Shifters
bit_values = [2,4,6];

for bits = bit_values

phase_step = 360/(2^bits);

phase_quantized = round(phase_ideal/phase_step)*phase_step;
phase_quantized = mod(phase_quantized,360);

w_quantized = exp(1j*deg2rad(phase_quantized));
AF_quantized = zeros(size(theta));

phase_error = mod(phase_quantized - phase_ideal + 180,360) - 180;

fprintf('---%d-bit Phase Shifter---\n', bits);

fprintf('Element | Ideal | Quantized | Error\n');

for m = 1:N
    fprintf('%d       |%7.2f| %7.2f   |%6.2f\n', ...
        m,phase_ideal(m),phase_quantized(m),phase_error(m));
end

for m = 1:N
    AF_quantized = AF_quantized + ...
        w_quantized(m).*exp(1j*n(m)*k*d*sin(theta));
end

AF_quantized = abs(AF_quantized);
AF_quantized = AF_quantized/max(AF_quantized);
AF_quantized_dB = 20*log10(AF_quantized);

[~, idx_peak] = max(AF_quantized);

theta_peak = theta_deg(idx_peak);

pointing_error = theta_peak - theta0;

fprintf('Peak angle: %.2f deg\n', theta_peak);
fprintf('Pointing error: %.2f deg\n\n', pointing_error);

figure(fig1);

plot(theta_deg,AF_quantized_dB, ...
    'LineWidth',1.8, ...
    'DisplayName', sprintf('%d-bit', bits));

hold on;

end

plot(theta_deg,AF_ideal_dB, ...
    'LineWidth',1.8, ...
    'DisplayName','Ideal');

grid on;
xlabel('Angle (degrees)');
ylabel('Normalized AF (dB)');
title('Effect of Phase Quantization');

xlim([-90 90]);
ylim([-40 0]);

legend('show');

% Random Phase/Amplitude Mismatch
sigma_A = 0.2;
sigma_phi = 20;

num_trials = 5000;

pointing_errors = zeros(1,num_trials);
PSLL_values = zeros(1,num_trials);

for trial = 1:num_trials

    A_e = 1 + sigma_A*randn(1,N);
    phase_e = sigma_phi*randn(1,N);

    w_error = A_e .* ...
        exp(1j*deg2rad(phase_ideal + phase_e));

    AF_e = zeros(size(theta));

    for m = 1:N
        AF_e = AF_e + ...
            w_error(m).*exp(1j*n(m)*k*d*sin(theta));
    end

    AF_e = abs(AF_e);
    AF_e = AF_e/max(AF_e);

    [~,idx_peak] = max(AF_e);
    theta_peak = theta_deg(idx_peak);

    pointing_errors(trial) = theta_peak - theta0;

    AF_e_dB = 20*log10(abs(AF_e));

    [pks, locs] = findpeaks(AF_e_dB);

    [~, main_idx] = max(pks);

    pks(main_idx) = [];
    locs(main_idx) = [];

    [max_secondary_lobe, idx] = max(pks);

    PSLL_values(trial) = max_secondary_lobe;
    
    if trial == 1

        figure(fig2);
        plot(rad2deg(theta),AF_e_dB, ...
            'LineWidth',1.8, ...
            'DisplayName', sprintf('With mismatch'));

        hold on;

        plot(rad2deg(theta),AF_ideal_dB, ...
            'LineWidth',1.8, ...
            'DisplayName', sprintf('Ideal'));
    end

end

fprintf('---Pointing Erros---\n')

fprintf('Mean pointing error: %.3f deg\n', ...
    mean(pointing_errors));

fprintf('RMS pointing error: %.3f deg\n', ...
    rms(pointing_errors));

fprintf('Worst observed pointing error: %.3f deg\n', ...
    max(abs(pointing_errors)));

fprintf('\n---PSLL---\n')

fprintf('Mean PSLL: %.3f dB\n', ...
    mean(PSLL_values));

fprintf('PSLL standard deviation: %.3f dB\n', ...
    std(PSLL_values));

fprintf('Worst observed PSLL: %.3f dB\n', ...
    -min(abs(PSLL_values)));

% Plot settings
figure(fig2);
grid on;
xlabel('Angle (degrees)');
ylabel('Amplitude (dB)');
title('Effect of random Phase (σA = 0.2*Α) -- Amplitude missmatch (σphi = 20°)');
sigma_A = 0.2;
sigma_phi = 20;
xlim([-90 90]);
ylim([-40 0]);
legend('show');

figure(fig3);
histogram(pointing_errors,30);

grid on;
xlabel('Pointing Error (degrees)');
ylabel('Number of Trials');
title('Monte Carlo Pointing Error');

figure(fig4);
histogram(PSLL_values,30);

grid on;
xlabel('Max Sidelobe (dB)');
ylabel('Number of Trials');
title('Monte Carlo PSLL');