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

visible = theta_deg >= -90 & theta_deg <= 90;

% Worst Case Side Lobes for varying BWs & Steering Angles
theta0_values = 0:1:60;
BW_values = [0.05 0.1 0.2 0.3 0.5];

theta0_TTD_sweep = deg2rad(theta0_values);

i = 1;

for BW = BW_values

    j = 1;

    for theta0 = theta0_TTD_sweep

        fmin = f0 * (1 - BW/2);
        fmax = f0 * (1 + BW/2);

        worst_secondary_lobe = -Inf;
        worst_secondary_angle = NaN;

        for f = fmin:0.01*f0:fmax

            lambda_f = c / f;

            AF_TTD_f = sin((N*pi*d*(sin(theta) - sin(theta0))/lambda_f)) ./ ...
                (N.*sin((pi*d*(sin(theta) - sin(theta0))/lambda_f)));

            singular = abs(sin(theta) - sin(theta0)) < 1e-12;
            AF_TTD_f(singular) = 1;

            AF_TTD_f = AF_TTD_f./max(AF_TTD_f);
            AF_TTD_f_dB = 20*log10(abs(AF_TTD_f));

            theta_visible = theta(visible);
            AF_TTD_f_dB_visible = AF_TTD_f_dB(visible);

            [pks, locs] = findpeaks(AF_TTD_f_dB_visible);

            [~, main_idx] = max(pks);

            pks(main_idx) = [];
            locs(main_idx) = [];

            [max_secondary_lobe, idx] = max(pks);

            if max_secondary_lobe > worst_secondary_lobe
                worst_secondary_lobe = max_secondary_lobe;
            end

        end
        worst_case_lobe(i,j) = worst_secondary_lobe;

        j = j + 1;
    end
    i = i + 1;
end

figure;

imagesc(theta0_values, 1:length(BW_values), worst_case_lobe);

set(gca, 'YDir', 'normal');

yticks(1:length(BW_values));
yticklabels(BW_values*100);

colorbar;

xlabel('Steering Angle (degrees)');
ylabel('Fractional Bandwidth (%)');
title('Worst-Case Sidelobe Level for Ideal TTD Beamforming');

cb = colorbar;
ylabel(cb, 'Worst-Case Sidelobe Level (dB)');

clim([-13 0]);


% Worst Case PS Pointing Error for varying BWs & Steering Angles
theta0_values = 0:1:60;
BW_values = [0.05 0.1 0.2 0.3 0.5];
worst_case_PS = zeros(length(BW_values), length(theta0_values));

i = 1;

for BW = BW_values

    j = 1;

    for theta0 = theta0_values

        theta0_rad = deg2rad(theta0);
        fmin = f0 * (1 - BW/2);
        fmax = f0 * (1 + BW/2);
        f_sweep = fmin:0.01*f0:fmax;
        argument = (f0./f_sweep) * sin(theta0_rad);

        if any(abs(argument) > 1)
            worst_case_PS(i,j) = NaN;
        else
            theta_peak = asin(argument);
            pointing_error = abs(rad2deg(theta_peak) - theta0);
            worst_case_PS(i,j) = max(pointing_error);
        end

        j=j+1;
    end

    i = i + 1;
end

figure('Color','k');

h = imagesc(theta0_values, 1:length(BW_values), worst_case_PS);

set(h,'AlphaData',~isnan(worst_case_PS));
set(gca,'YDir','normal','Color',[0.3 0.3 0.3], ...
    'XColor','w','YColor','w');

yticks(1:length(BW_values));
yticklabels(BW_values*100);

xlabel('Steering Angle (degrees)');
ylabel('Fractional Bandwidth (%)');
title('Worst-Case Phase-Shifter Pointing Error');

cb = colorbar;
cb.Color = 'w';
ylabel(cb,'Pointing Error (degrees)');

clim([0 34]);

text(59,0.65,'gray: no real analytical steering solution', ...
    'Color','w','HorizontalAlignment','right');