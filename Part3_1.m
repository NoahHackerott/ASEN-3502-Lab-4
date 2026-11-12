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





