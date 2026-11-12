%% Part 2.4
clear; close all; clc;
addpath('Lab-4');   % evaluate.p

% Deterministic optimum (same setup as Task 2.2)
params = default_params();
% nominal parameters
theta  = [2.0, 0.0325, 648000, 4.0];   

objective = @(x) -uav_range(x, theta, params)/1e5;
% Use deal as a way to essentially create an anonymous function that
% returns two outputs
nonlcon   = @(x) deal(uav_constraints(x, theta, params), []);

% Bounds
lb = [0.3, 4, 0.2, 8];
ub = [3, 20, 10, 40];

options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-12, 'ConstraintTolerance', 1e-10, ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000);

% Multiple starts, keep the best
x0 = [0.5, 12, 4,   20;
      1.0,  8, 2,   15;
      0.4, 18, 5.5, 22;
      2.0, 10, 7,   30];

best_R = -Inf;
for i = 1:size(x0, 1)
    [x_opt, ~, exitflag] = fmincon(objective, x0(i,:), [], [], [], [], lb, ub, nonlcon, options);
    [R, g] = uav_range(x_opt, theta, params);
    if exitflag > 0 && all(g <= 1e-8) && R > best_R
        best_R = R;
        x_star = x_opt;
    end
end

fprintf('Deterministic optimum: S = %.6f, AR = %.6f, m_b = %.6f, V = %.6f\n', x_star);
fprintf('R* = %.3f km\n', best_R/1e3);

%% Challenge submission
evaluate(x_star, 'noah.hackerott@colorado.edu')

%% Local function
function c = uav_constraints(x, theta, params)
    [~, g] = uav_range(x, theta, params);
    c = g(1,:);
end