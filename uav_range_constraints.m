function c = uav_range_constraints(x, theta, params)

    % Obtain constraints from the existing UAV model
    [~, g] = uav_range(x, theta, params);

    % Return the four inequality constraints
    c = g(1,:);

end