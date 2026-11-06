function theta = sample_theta(sc, N)
% Draw N random parameter sets for the UAV model.

% Payload mass [kg]
m_pay = sc.mpay_mu + sc.mpay_sigma * randn(N,1);

% Zero-lift drag coefficient
Cd0 = sc.CD0_min + (sc.CD0_max - sc.CD0_min) * rand(N,1);

% Battery specific energy [J/kg], normal N(mu, sigma)
eb = sc.eb_mu + sc.eb_sig * randn(N,1);

% Wind speed [m/s], uniform U(0, Vw_max)
V_w = sc.Vw_max * rand(N,1);

% Assemble in the same column order uav_range expects
theta = [m_pay, Cd0, eb, V_w];
end