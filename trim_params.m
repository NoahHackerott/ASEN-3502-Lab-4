function P = trim_params()
P = struct;
P.rho     = 1.225;   % kg/m^3
P.V       = 30;      % m/s
P.S       = 16.2;    % m^2
P.W       = 11346;   % N
P.q       = 0.5*P.rho*P.V^2;   % Pa
P.CL0     = 0.247;
P.CLalpha = 6.0;     % 1/rad
P.CLdelta = 0.43;    % 1/rad
P.CM0     = 0.1;
P.CMalpha = -1.3;    % 1/rad
P.CMdelta = -1.28;   % 1/rad
P.DCLp0   = 0.6;
P.alpha_p = 0.05;    % rad
end