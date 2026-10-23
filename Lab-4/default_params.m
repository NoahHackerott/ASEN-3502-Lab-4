function P = default_params()
P = struct;
P.g       = 9.81;
P.rho     = 1.225; % air density [kg/m^3]
P.m_fix   = 1.5; % avionics + motor + fuselage mass [kg]
P.k_w     = 0.35; % wing mass m_w = k_w * S^a_w * AR^b_w  [kg]

P.a_w     = 1.2; % exponent of S in m_w equation
P.b_w     = 0.9; % exponent of AR in m_w equation

P.eta_max = 0.70; % peak propulsive efficiency
P.V_eta   = 18.0; % airspeed at which the fixed-pitch prop is most efficient [m/s]
P.f_res   = 0.20; % battery reserve fraction
P.C_Lmax   = 1.3;
P.k_stall = 1.25; % required V/V_stall
P.P_max   = 500; % motor continuous power [W]
P.RC_min  = 2.0; % required rate of climb at cruise speed [m/s]
P.b_max   = 4.0; % span limit [m]
P.m_max   = 15.0; % MTOW [kg]

P.min = [0.3,4,0.2,8]; % bounds on x = [S, AR, m_b, V]
P.max = [3.0,20,10,40];
end
