function c = uav_constraints(x, theta, params)
    [~, g] = uav_range(x, theta, params);
    c = g(1,:);
end