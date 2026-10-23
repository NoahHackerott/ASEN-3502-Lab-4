function varargout = scenario(action, varargin)
% Distributions of the uncertain parameters theta = [m_pay, CD0, e_b, V_w].

switch action
    % sc = scenario('nominal') - the nominal scenario used in the lab before the challenge
    case 'nominal'
        sc = struct('scenario_id', 0, 'mpay_mu', 2.0, 'mpay_sig', 0.3, 'CD0_min', 0.025, 'CD0_max', 0.040, 'eb_mu', 648000, 'eb_sig', 54000, 'Vw_max', 8.0);
        varargout = {sc};

    % theta = scenario('sample', sc, n) - n pseudo-random rows of theta
    % theta = scenario('sample', sc, n, true) - n scrambled-Sobol rows of theta
    case 'sample'
        sc = varargin{1};
        n = varargin{2};
        use_sobol = numel(varargin) >= 3 && varargin{3};

        if use_sobol
            randVals = net(scramble(sobolset(4, 'Skip', 1e3), 'MatousekAffineOwen'), n);
        else
            randVals = rand(n, 4);
        end
        
        ninv = @(p) -sqrt(2) * erfcinv(2*p); % standard normal inverse CDF
        theta = [sc.mpay_mu + sc.mpay_sig*ninv(randVals(:,1)), ...
                 sc.CD0_min + (sc.CD0_max - sc.CD0_min)*randVals(:,2), ...
                 sc.eb_mu + sc.eb_sig*ninv(randVals(:,3)), ...
                 sc.Vw_max*randVals(:,4)];
        varargout = {theta};

    % theta = scenario('mean', sc) - the mean parameter vector (1 x 4)
    case 'mean'
        sc = varargin{1};
        varargout = {[sc.mpay_mu, 0.5*(sc.CD0_min + sc.CD0_max), sc.eb_mu, 0.5*sc.Vw_max]};
end

end
