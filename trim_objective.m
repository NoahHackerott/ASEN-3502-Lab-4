function f = trim_objective(alpha, deltaE, p)
% Trim objective f(alpha, delta_e)
% Vectorized: alpha/deltaE may be matrices
CL = p.CL0 + p.CLalpha*alpha + p.CLdelta*deltaE + p.DCLp0*exp(-(alpha/p.alpha_p).^2);
CM = p.CM0 + p.CMalpha*alpha + p.CMdelta*deltaE;
f  = ((p.q*p.S*CL - p.W)/p.W).^2 + CM.^2;
end