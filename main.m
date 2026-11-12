%% Part 1.1

params = trim_params();
p = params;


f = @(alpha, deltaE) ((p.q*p.S*(p.CL0 + p.CLalpha*alpha + p.CLdelta*deltaE ...
        + p.DCLp0*exp(-(alpha/p.alpha_p).^2)) - p.W)/p.W).^2 ...
        + (p.CM0 + p.CMalpha*alpha + p.CMdelta*deltaE).^2;

plot_contours(1, f)

%% Part 1.2
opts = optimset('TolX', 1e-12, 'TolFun', 1e-16);

% Choose x0 based off the contour plot
x0 = [deg2rad(10);deg2rad(-5)];

% Redefine f so that it will work with fminsearch
fsearch = @(x) f(x(1),x(2));

[x, fval] = fminsearch(fsearch,x0,opts);

% With help of Claude
fprintf('Task 1.2: Reference optimum (fminsearch)\n');
fprintf('  alpha*   = %.4f deg\n', rad2deg(x(1)));
fprintf('  deltaE*  = %.4f deg\n', rad2deg(x(2)));
fprintf('  f*       = %.4e\n', fval);

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

%% Part 2.1
p = default_params();


% Test points, each row is one case
X  = [0.5 12 4   20;
      1.0  8 2   15;
      0.4 18 5.5 22];
TH = [2.0 0.0325 648000 4;
      2.3 0.035  612000 6;
      1.7 0.028  684000 2];

R_test = zeros(3,1);
g_test = zeros(3,4);
for k = 1:3
    [R_test(k), g_test(k,:)] = uav_range(X(k,:), TH(k,:), p);
end

% With help of Claude
fprintf('Task 2.1: UAV model test points\n');
fprintf('%4s %14s %14s %14s %14s %14s\n', 'Case', 'R [m]', 'g1', 'g2', 'g3', 'g4');
for k = 1:3
    fprintf('%4d %14.6g %14.6g %14.6g %14.6g %14.6g\n', k, R_test(k), g_test(k,:));
end

% Vectorization check, same theta in 5 rows should give 5 identical rows
[R_vec, g_vec] = uav_range(X(1,:), repmat(TH(1,:), 5, 1), p);
assert(isequal(size(R_vec), [5 1]) && isequal(size(g_vec), [5 4]), 'Wrong output size')
assert(all(R_vec == R_vec(1)), 'Rows of R do not match')

%% Part 2.2
clc;
clear;
close all;


params = default_params();


theta = [2.0, 0.0325, 648000, 4.0];


%% 2. Define objective function



objective = @(x) -uav_range(x, theta, params)/(1e5);


%% 3. Define design variable bounds

% x = [S, AR, mb, V]

lb = [0.3, 4, 0.2, 8];    
ub = [3, 20, 10, 40];     



nonlcon = @(x) deal( ...
    uav_range_constraints(x, theta, params), []);


%% 5. Set optimization options

options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-10, ...
    'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 5000);



x0 = [0.5, 12, 4,   20;
      1.0,  8, 2,   15;
      0.4, 18, 5.5, 22;
      2.0, 10, 7,   30 ];



results = zeros(4,5);

for i = 1:4

    % Run fmincon
    [x_opt, fval, exitflag] = fmincon( ...
        objective, x0(i,:), ...
        [], [], [], [], lb, ub, nonlcon, options);

    % Calculate actual range in meters
    [R, g] = uav_range(x_opt, theta, params);

    % Convert range to kilometers
    R_km = R/1000;

    % Save results
    results(i,:) = [x_opt, R_km];

    % Display results
    fprintf('\n====================================\n');
    fprintf('Starting Point %d\n', i);
    fprintf('====================================\n');

    fprintf('Wing Area:       %.8f m^2\n', x_opt(1));
    fprintf('Aspect Ratio:    %.8f\n', x_opt(2));
    fprintf('Battery Mass:    %.8f kg\n', x_opt(3));
    fprintf('Cruise Speed:    %.8f m/s\n', x_opt(4));
    fprintf('Maximum Range:   %.8f km\n', R_km);

    fprintf('Exit Flag:       %d\n', exitflag);

    fprintf('\nConstraint Values:\n');
    disp(g);

    % Identify active constraints
    active = abs(g) < 1e-6;

    fprintf('Active Constraints [g1 g2 g3 g4]:\n');
    disp(active);

end



fprintf('\n====================================\n');
fprintf('FINAL COMPARISON\n');
fprintf('====================================\n');

disp('Columns: [S, AR, mb, V, Range (km)]');
disp(results);




[best_range, best_index] = max(results(:,5));

fprintf('\nBest Starting Point: %d\n', best_index);
fprintf('Best Range: %.8f km\n', best_range);

fprintf('Optimal Design [S, AR, mb, V]:\n');
disp(results(best_index,1:4));

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


%% 3.1 Adam Function implementation
% These are given from teh first week parameters
params = trim_params();


x = [0.2;0.1];
l = 0.02;
beta1 = 0.3;
beta2 = 0.99;
epsilon = 1e-8;
tol = 1e-6;
maxit=1000;
f = @(a,d) trim_objective(a, d, params);
% Function used for trim gradient
g = @(x) trim_gradient(x(1), x(2), params);

[position, path] = Adam(g, x, tol, maxit, l, beta2, beta1, epsilon);

% Reference optimum from Task 1.2 for the distance comparison
% Claude helped with formatting
opts  = optimset('TolX', 1e-12, 'TolFun', 1e-16);
x_ref = fminsearch(@(z) trim_objective(z(1), z(2), params), [deg2rad(10); deg2rad(-5)], opts);

fprintf('Task 3.1: ADAM from x0 = [0.2; 0.1]\n');
fprintf('  alpha  = %.4f deg\n', rad2deg(position(1)));
fprintf('  deltaE = %.4f deg\n', rad2deg(position(2)));
fprintf('  iterations = %d\n', size(path,2) - 1);
fprintf('  distance from reference = %.4e rad\n', norm(position - x_ref));


plot_contours(1,f);

hold on;
plot(path(1,:)*180/pi,path(2,:)*180/pi,'r-o','LineWidth',2);
plot(position(1)*180/pi,position(2)*180/pi,'bo','MarkerFaceColor','b');
legend('Objective Function','ADAM Path','Final Position');
hold off;

%% Part 3.2

params = trim_params();

f2 = @(a,d) ((8930.2*(0.247 + 6*a + 0.43*d + 0.6*exp(-(a/0.05).^2)) - 11346)/11346).^2 + (0.1 - 1.3*a - 1.28*d).^2;
f = @(x) f2(x(1),x(2));
g = @(x) trim_gradient(x(1), x(2), params);

limits = deg2rad([10; -5]);
tol = 1e-6;

[position,path, lambda] = augmented_lagrange(g, [0.15; -0.1], tol, 1000, f, limits);

plot_contours(1,f2);
hold on;
plot(rad2deg(path(1,:)),rad2deg(path(2,:)), ...
    'b.-','DisplayName','Augmented Lagrangian');
xline(10,'r--','DisplayName','Alpha limit');
yline(-5,'m--','DisplayName','Elevator limit');
legend('Location','best');
hold off;

c = position - limits;

fprintf('Alpha: %.6f deg\n',rad2deg(position(1)));
fprintf('Elevator: %.6f deg\n',rad2deg(position(2)));
fprintf('Total iterations: %d\n',size(path,2)-1);
fprintf('Constraints [rad]: %.8e, %.8e\n',c(1),c(2));
fprintf('Active constraints [alpha, elevator]: %d %d\n', lambda > 0);
fprintf('Both satisfied within tolerance: %d\n', all(c <= tol));

outside = find(any(path > limits + tol,1));
fprintf('Number of infeasible iterations: %d\n',numel(outside));
% Chat helped with plot/printing formatting :)


if ~isempty(outside)
    j = outside(1);
    fprintf('First infeasible iteration: %d, at [%.6f, %.6f] deg\n', ...
        j-1,rad2deg(path(1,j)),rad2deg(path(2,j)));
end






