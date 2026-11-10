function G = trim_gradient(alpha, deltaE, params)
% Analytic gradient of the trim objective
% Returns G = [df/dalpha; df/ddeltaE] (2x1), angles in radians

% Parameters
q       = params.q;        % Pa
S       = params.S;        % m^2
W       = params.W;        % N
CL0     = params.CL0;
CLalpha = params.CLalpha;  % 1/rad
CLdelta = params.CLdelta;  % 1/rad
DCLp0   = params.DCLp0;
alpha_p = params.alpha_p;  % rad
CM0     = params.CM0;
CMalpha = params.CMalpha;  % 1/rad
CMdelta = params.CMdelta;  % 1/rad

% Model terms at (alpha, deltaE)
CL = CL0 + CLalpha*alpha + CLdelta*deltaE + DCLp0*exp(-(alpha/alpha_p)^2);
CM = CM0 + CMalpha*alpha + CMdelta*deltaE;   % pitching moment coefficient

% Chain rule on each squared term: df/dx = 2*L*dL/dx + 2*CM*dCM/dx
% Lift term contributions: 2*L*(q*S/W)*dCL/dx
lift_dalpha  = 2 * ((q*S/W)*CL - 1) * (q*S/W)*(CLalpha + DCLp0*exp(-(alpha/alpha_p)^2)*(-2*alpha/alpha_p^2));
lift_ddeltaE = 2 * ((q*S/W)*CL - 1) * (q*S/W)*CLdelta;

% Moment term contributions: 2*CM*dCM/dx
moment_dalpha  = 2 * CM * CMalpha;
moment_ddeltaE = 2 * CM * CMdelta;

dfdalpha  = lift_dalpha  + moment_dalpha;
dfddeltaE = lift_ddeltaE + moment_ddeltaE;

G = [dfdalpha; dfddeltaE];
end