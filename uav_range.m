function [R,g] = uav_range(x, theta, params)
p = params;

% Defined variables
S = x(1); % m^2
AR = x(2); % dimensionless
m_b = x(3); % kg
V = x(4); % m/s
m_pay = theta(:,1); % kg
Cd0 = theta(:,2); % Dimensionless zero lift drag coefficient
eb = theta(:,3); % J/kg
V_w = theta(:,4); % m/s

m_w = p.k_w * S^(p.a_w) * AR^(p.b_w);

W = p.g * (m_pay + p.m_fix + m_b + m_w);

q = 0.5 * p.rho * V^2;

e = 1.78 * (1 - 0.045 * AR^0.68) - 0.64;

D = q * S * Cd0 + W.^2/(q * S * pi * e * AR);

r = V / p.V_eta;

eta = p.eta_max * 2 * r / (1 + r^2);

E = (1 - p.f_res) * eb * m_b * eta ./ (D * V);

R = E * V .* (1 - V_w.^2/V^2);

g(:,1) = p.k_stall * sqrt(2*W/(p.rho * S * p.C_Lmax)) - V;
g(:,2) = p.RC_min - (eta * p.P_max - D * V)./W;
g(:,3) = sqrt(AR * S) - p.b_max;
g(:,4) = m_pay + p.m_fix + m_b + m_w - p.m_max; 


end