params = trim_params();

f2 = @(a,d) ((8930.2*(0.247 + 6*a + 0.43*d + 0.6*exp(-(a/0.05).^2)) - 11346)/11346).^2 + (0.1 - 1.3*a - 1.28*d).^2;
f = @(x) f2(x(1),x(2));
g = @(x) trim_gradient(x(1), x(2), params);

limits = deg2rad([10; -5]);
tol = 1e-6;

[position,path, lambda] = augmented_lagrange(g, [0.15; -0.1], tol, 1000, f, limits);

plot_contours(1,f2);
hold on;
plot(rad2deg(path(1,:)),rad2deg(path(2,:)), ...
    'b.-','DisplayName','Augmented Lagrangian');
xline(10,'r--','DisplayName','Alpha limit');
yline(-5,'m--','DisplayName','Elevator limit');
legend('Location','best');
hold off;

c = position - limits;

fprintf('Alpha: %.6f deg\n',rad2deg(position(1)));
fprintf('Elevator: %.6f deg\n',rad2deg(position(2)));
fprintf('Total iterations: %d\n',size(path,2)-1);
fprintf('Constraints [rad]: %.8e, %.8e\n',c(1),c(2));
fprintf('Active constraints [alpha, elevator]: %d %d\n', lambda > 0);
fprintf('Both satisfied within tolerance: %d\n', all(c <= tol));

outside = find(any(path > limits + tol,1));
fprintf('Number of infeasible iterations: %d\n',numel(outside));
% Chat helped with plot/printing formatting :)


if ~isempty(outside)
    j = outside(1);
    fprintf('First infeasible iteration: %d, at [%.6f, %.6f] deg\n', ...
        j-1,rad2deg(path(1,j)),rad2deg(path(2,j)));
end


