%% Part 1.1
f = @(alpha, deltaE) ((8.9302e3*(0.247 + 6*alpha + 0.43*deltaE + 0.6*exp(-(alpha/0.05).^2)) - 11346)/11346).^2 ...
                   + (0.1 - 1.3*alpha - 1.28*deltaE).^2;
plot_contours(1, f)

%% Part 1.2
opts = optimset('TolX', 1e-12, 'TolFun', 1e-16);

[x, fval] = fminsearch(f, );