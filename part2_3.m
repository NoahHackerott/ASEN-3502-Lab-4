%% Part 2.3
clear; close all; clc;

sc = scenario("nominal");

N = 10000;
theta = sample_theta(sc, N);   % N x 4: [m_pay, CD0, e_b, V_w]

% Histograms of the sampled parameters (Claude assisted with formatting)
names  = {'Payload mass m_{pay}', 'Zero-lift drag C_{D0}', ...
          'Battery specific energy e_b', 'Wind speed V_w'};
xlabs  = {'m_{pay} [kg]', 'C_{D0} [-]', 'e_b [J/kg]', 'V_w [m/s]'};

figure(1); clf
for k = 1:4
    subplot(2, 2, k)
    histogram(theta(:, k), 50)
    xlabel(xlabs{k})
    ylabel('Count')
    title(names{k})
    grid on
end
sgtitle(sprintf('Task 2.3: Sampled \\theta (nominal scenario), N = %d', N))