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