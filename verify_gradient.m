%% Verify Script 1.3
clc; clear;
params = trim_params();

% Objective (same f used everywhere else)
g = @(x) trim_objective(x(1), x(2), params);

% Test points [alpha; deltaE] in rad, one per column
X_test = [0.1  0.05;
          0.0 -0.10;
          0.2  0.10]';

for j = 1:size(X_test, 2)
    x = X_test(:, j);
    alpha = x(1);
    deltaE = x(2);

    G_analytic = trim_gradient(alpha, deltaE, params);

    % Forward difference with h_k = sqrt(eps)*|x_k|
    G_numeric = zeros(2,1);
    for i = 1:2
        h = sqrt(eps)*abs(x(i));
        if h == 0
            h = sqrt(eps);   % avoid h = 0 when x_k = 0
        end
        dx = zeros(2,1);
        dx(i) = h;
        G_numeric(i) = ...
            (g(x + dx) - g(x))/h;
    end

    fprintf('\nTest point: alpha = %.3f rad, deltaE = %.3f rad\n', alpha, deltaE);
    disp('Analytical gradient:');
    disp(G_analytic);
    disp('Numerical gradient:');
    disp(G_numeric);

    error_value = norm(G_analytic - G_numeric)/norm(G_analytic);
    fprintf('Relative error norm: %.12e\n', error_value);
    if error_value < 1e-5
        disp('PASS: Gradient is verified.');
    else
        disp('FAIL: Check your gradient function.');
    end
end