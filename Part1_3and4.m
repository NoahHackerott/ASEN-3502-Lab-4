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

%% 1.4a
% Claude helped with formatting
[position_a, path_a] = gradient_descent(G, [0; 0], tol, max_iterations, fsearch, 0.01);
fprintf('\nTask 1.4a: x0 = [0; 0], l_init = 0.01\n');
fprintf('  alpha  = %.4f deg\n', rad2deg(position_a(1)));
fprintf('  deltaE = %.4f deg\n', rad2deg(position_a(2)));
fprintf('  f      = %.4e,  |grad f| = %.3e\n', fsearch(position_a), norm(G(position_a)));
fprintf('  iterations = %d\n', size(path_a,2) - 1);
plot_contours(2, f); hold on
plot(rad2deg(path_a(1,:)), rad2deg(path_a(2,:)), 'r.-', 'DisplayName', 'GD path (a)');
plot(rad2deg(path_a(1,1)), rad2deg(path_a(2,1)), 'ko', 'DisplayName', 'start');
plot(rad2deg(path_a(1,end)), rad2deg(path_a(2,end)), 'r*', 'MarkerSize', 10, 'DisplayName', 'end');
plot(rad2deg(x_ref(1)), rad2deg(x_ref(2)), 'kp', 'MarkerSize', 12, 'DisplayName', 'fminsearch optimum');
title('Task 1.4a: Gradient descent from [0, 0]'); legend('Location', 'best'); hold off


%% 1.4b
% Claude helped with fprintf formatting
[position_b, path_b] = gradient_descent(G, [0.2; 0.1], tol, max_iterations, fsearch, 0.045);
fprintf('\nTask 1.4b: x0 = [0.2; 0.1], l_init = 0.045\n');
fprintf('  alpha  = %.4f deg\n', rad2deg(position_b(1)));
fprintf('  deltaE = %.4f deg\n', rad2deg(position_b(2)));
fprintf('  iterations = %d\n', size(path_b,2) - 1);
fprintf('  distance from reference = %.4e rad\n', norm(position_b - x_ref));


%% 1.5  Newton's method with backtracking
% Claude helped with formatting

H = @(x) trim_hessian(x(1), x(2), params);
[position_n, path_n] = newtons_method(G, [0.2; 0.1], tol, max_iterations, fsearch, H);
fprintf('\nTask 1.5: Newton from x0 = [0.2; 0.1]\n');
fprintf('  alpha  = %.4f deg\n', rad2deg(position_n(1)));
fprintf('  deltaE = %.4f deg\n', rad2deg(position_n(2)));
fprintf('  iterations = %d\n', size(path_n,2) - 1);
fprintf('  distance from reference = %.4e rad\n', norm(position_n - x_ref));

% Combined path plot
plot_contours(3, f); hold on
plot(rad2deg(path_b(1,:)), rad2deg(path_b(2,:)), 'r.-', 'DisplayName', 'Gradient descent (1.4b)');
plot(rad2deg(path_n(1,:)), rad2deg(path_n(2,:)), 'b.-', 'DisplayName', 'Newton');
plot(rad2deg(path_b(1,1)), rad2deg(path_b(2,1)), 'ko', 'DisplayName', 'start');
plot(rad2deg(path_b(1,end)), rad2deg(path_b(2,end)), 'r*', 'MarkerSize', 10, 'DisplayName', 'GD end');
plot(rad2deg(path_n(1,end)), rad2deg(path_n(2,end)), 'b*', 'MarkerSize', 10, 'DisplayName', 'Newton end');
plot(rad2deg(x_ref(1)), rad2deg(x_ref(2)), 'kp', 'MarkerSize', 12, 'DisplayName', 'fminsearch optimum');
title('Task 1.5: Gradient descent vs. Newton from [0.2, 0.1] rad'); legend('Location', 'best'); hold off

% Comparison table
Method     = ["Gradient descent"; "Newton"];
alpha_deg  = rad2deg([position_b(1); position_n(1)]);
deltaE_deg = rad2deg([position_b(2); position_n(2)]);
Iterations = [size(path_b,2) - 1; size(path_n,2) - 1];
Distance_rad = [norm(position_b - x_ref); norm(position_n - x_ref)];
results_15 = table(Method, alpha_deg, deltaE_deg, Iterations, Distance_rad);
disp(results_15)