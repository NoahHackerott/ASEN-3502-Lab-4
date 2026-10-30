clc
clear
% Lab 3 Week 1

%% 1.3 
function G = trim_gradient(alpha, deltaE,params)
rho = params.rho; %kg/m^2
CL0 = params.CL0;
CM0 = params.CM0;
DCLp0 = params.DCLp0;
V = params.V; %m/s
S = params.S; %m^2
W = params.W; %N
CLal = params.CLal; %rad^-1
CMal = params.CMal; %rad^-1
alp = params.alp; %rad
CLdel = params.CLdel; %rad^-1
CMdel = params.CMdel; %rad^-1

syms a de real 
f = ((0.5*rho*V^2*S * (CL0 + CLal*a + CLdel*de + CLal*exp(-(a/alp)^2)) - W)/W)^2 + (CM0 - CMal*a - CMdel*de)^2;

Gsym = [diff(f,a); diff(f,de)];


    G = double(subs(Gsym, [a, de], [alpha, deltaE]));
end

%% Verify Script 1.3
clc
clear

params = struct();
params.rho = 1.225; %kg/m^2
params.CL0 = 0.247;
params.CM0 = 0.1;
params.DCLp0 = 0.6;
params.V = 30; %m/s
params.S = 16.3; %m^2
params.W = 11346; %N
params.CLal = 6.0; %rad^-1
params.CMal = -1.3; %rad^-1
params.alp = 0.05; %rad
params.CLdel = 0.43; %rad^-1
params.CMdel = -1.28; %rad^-1

x = [0.1; 0.05];
alpha = x(1);
deltaE = x(2);


G_analytic = trim_gradient(alpha, deltaE, params);

% Scalar 
g = @(x) ((0.5*params.rho*params.V^2*params.S *(params.CL0 + params.CLal*x(1) ...
    + params.CLdel*x(2) + params.CLal*exp(-(x(1)/params.alp)^2)) ...
    - params.W)/params.W)^2 + (params.CM0 - params.CMal*x(1) ...
    - params.CMdel*x(2))^2;


h = sqrt(eps)*norm(x);
G_numeric = zeros(2,1);

for i = 1:2
    dx = zeros(2,1);
    dx(i) = h;

    G_numeric(i) = ...
        (g(x + dx) - g(x - dx))/(2*h);
end


disp('Analytical gradient:');
disp(G_analytic);

disp('Numerical gradient:');
disp(G_numeric);

error_value = norm(G_analytic - G_numeric);
fprintf('Absolute error norm: %.12e\n', error_value);

if error_value < 1e-5
    disp('PASS: Gradient is verified.');
else
    disp('FAIL: Check your gradient function.');
end
%% 1.4

function [position, path] = gradient_descent(g, initial_guess, tol, max_iterations, l_init)

    rho = 0.5;
    mu1 = 1e-4;

    position = initial_guess(:);
    path = position;

    for k = 1:max_iterations

        % Approximate the gradient using finite differences
        n = length(position);
        grad = zeros(n,1);
        eps_fd = 1e-6;

        g0 = g(position);

        for j = 1:n
            ej = zeros(n,1);
            ej(j) = eps_fd;
            grad(j) = (g(position + ej) - g(position - ej)) ...
                      / (2*eps_fd);
        end

        % Steepest descent direction
        p = -grad;

        % Stop if the gradient is essentially zero
        if norm(grad) < 1e-10
            break;
        end

        % Backtracking line search
        l = l_init;

        while l >= tol
            x_trial = position + l*p;

            if g(x_trial) <= g0 + mu1*l*(grad'*p)
                break;
            end

            l = rho*l;
        end

        % Stop if no acceptable step was found
        if l < tol
            break;
        end

        % Accept the step and save the new position
        position = x_trial;
        path(:,end+1) = position;

    end

    fprintf('Final alpha = %.10f rad\n', position(1));
    fprintf('Final delta_e = %.10f rad\n', position(2));
    fprintf('Number of iterations = %d\n', size(path,2)-1);
end

%% 1.4a
g = @(x) ((0.5*params.rho*params.V^2*params.S * ...
    (params.CL0 + params.CLal*x(1) + params.CLdel*x(2) ...
    + params.CLal*exp(-(x(1)/params.alp)^2)) ...
    - params.W)/params.W)^2 ...
    + (params.CM0 - params.CMal*x(1) ...
    - params.CMdel*x(2))^2;

G = @(x) trim_gradient(x(1), x(2), params);

tol = 1e-6;
max_iterations = 1000;
initial_guess = [0; 0];
l_init = 0.01;

[position_a, path_a] = gradient_descent( ...
    g, initial_guess, tol, max_iterations, l_init);

fprintf('Part (a) iterations: %d\n', size(path_a,2)-1);
disp('Part (a) final position [alpha; deltaE]:');
disp(position_a);

%% 1.4 b
initial_guess = [0.2; 0.1];
l_init = 0.045;

[position_b, path_b] = gradient_descent( ...
    g, initial_guess, tol, max_iterations, l_init);

fprintf('Part (b) iterations: %d\n', size(path_b,2)-1);
disp('Part (b) final position [alpha; deltaE]:');
disp(position_b);
