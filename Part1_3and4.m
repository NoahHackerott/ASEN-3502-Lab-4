%% 1.4  Gradient descent with backtracking line search
params = trim_params();
p = params;
% two-argument form (plot_contours)
f = @(alpha, deltaE) ((p.q*p.S*(p.CL0 + p.CLalpha*alpha + p.CLdelta*deltaE ...
        + p.DCLp0*exp(-(alpha/p.alpha_p).^2)) - p.W)/p.W).^2 ...
        + (p.CM0 + p.CMalpha*alpha + p.CMdelta*deltaE).^2;  
% vector form (methods)
fsearch = @(x) f(x(1), x(2));   
% gradient handle
G = @(x) trim_gradient(x(1), x(2), params);                   

% Reference optimum (same call as Task 1.2), needed for distances
opts  = optimset('TolX', 1e-12, 'TolFun', 1e-16);
x_ref = fminsearch(fsearch, [deg2rad(10); deg2rad(-5)], opts);

tol = 1e-6;
max_iterations = 1000;

% TODO: update gradient_descent.m to this interface:
%   [position, path] = gradient_descent(g, initial_guess, tol, max_iterations, f, l_init)
%   - g = gradient handle (use g(x), delete the finite-difference gradient loop)
%   - f = objective, only for the sufficient-decrease test (rho = 0.5, mu1 = 1e-4)
%   - stop when norm(x_new - x_old) < tol, or after max_iterations
%   - path is 2 x (k+1) and includes initial_guess as its first column

%% 1.4

% Claude helped with formatting
[position_a, path_a] = gradient_descent(G, [0; 0], tol, max_iterations, f, 0.01);
fprintf('\nTask 1.4a: x0 = [0; 0], l_init = 0.01\n');
fprintf('  alpha  = %.4f deg\n', rad2deg(position_a(1)));
fprintf('  deltaE = %.4f deg\n', rad2deg(position_a(2)));
fprintf('  f      = %.4e,  |grad f| = %.3e\n', f(position_a), norm(G(position_a)));
fprintf('  iterations = %d\n', size(path_a,2) - 1);

plot_contours(2, f2); hold on
plot(rad2deg(path_a(1,:)), rad2deg(path_a(2,:)), 'r.-', 'DisplayName', 'GD path (a)');
plot(rad2deg(path_a(1,1)), rad2deg(path_a(2,1)), 'ko', 'DisplayName', 'start');
plot(rad2deg(path_a(1,end)), rad2deg(path_a(2,end)), 'r*', 'MarkerSize', 10, 'DisplayName', 'end');
plot(rad2deg(x_ref(1)), rad2deg(x_ref(2)), 'kp', 'MarkerSize', 12, 'DisplayName', 'fminsearch optimum');
title('Task 1.4a: Gradient descent from [0, 0]'); legend('Location', 'best'); hold off

% TODO (writeup): Is this a minimum? Is it a trim point?
%   Look at |grad f| (minimum?) and at f itself (trim requires f = 0).

%% 1.4b
[position_b, path_b] = gradient_descent(G, [0.2; 0.1], tol, max_iterations, f, 0.045);
fprintf('\nTask 1.4b: x0 = [0.2; 0.1], l_init = 0.045\n');
fprintf('  alpha  = %.4f deg\n', rad2deg(position_b(1)));
fprintf('  deltaE = %.4f deg\n', rad2deg(position_b(2)));
fprintf('  iterations = %d\n', size(path_b,2) - 1);
fprintf('  distance from reference = %.4e rad\n', norm(position_b - x_ref));
% (path_b gets plotted together with Newton's path in 1.5)

%% 1.5  Newton's method with backtracking
% TODO: write newtons_method.m with the same interface:
%   [position, path] = newtons_method(g, initial_guess, tol, max_iterations, f, H, ...)
%   with H = @(x) trim_hessian(x(1), x(2), params)
% TODO: run from [0.2; 0.1], report your line-search parameters,
%       plot path_b and the Newton path on one contour plot,
%       and build the comparison table (alpha, deltaE in deg, iterations, norm(x - x_ref))