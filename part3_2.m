f2 = @(a,d) ((8930.2*(0.247 + 6*a + 0.43*d + 0.6*exp(-(a/0.05).^2)) - 11346)/11346).^2 + (0.1 - 1.3*a - 1.28*d).^2;
f = @(x) f2(x(1),x(2));
g = @trim_gradient;

limits = deg2rad([10; -5]);
tol = 1e-6;

[position,path] = augmented_lagrange(g, [0.15; -0.1], tol, 1000, f, limits);

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
fprintf('Both satisfied within tolerance: %d\n',all(c <= tol));
fprintf('Active constraints [alpha, elevator]: %d %d\n',abs(c) <= tol);

outside = find(any(path > limits + tol,1));
fprintf('Number of infeasible iterates: %d\n',numel(outside));
% Chat helped with plot/printing formatting :)


if ~isempty(outside)
    j = outside(1);
    fprintf('First infeasible iterate: %d, at [%.6f, %.6f] deg\n', ...
        j-1,rad2deg(path(1,j)),rad2deg(path(2,j)));
end

function [position,path] = augmented_lagrange(g,initial_guess,tol,max_iterations,f,limits)

    position = initial_guess;
    path = position;
    lambda = [0;0];
    mu = 10;

    for k = 1:20
        previous = position;

        phi = @(x) f(x) + sum(max(0,lambda + mu*(x-limits)).^2 - lambda.^2)/(2*mu);
        grad_phi = @(x) g(x) + max(0,lambda + mu*(x-limits));

        % PUT ADAM HERE
        [position,inner_path] = adam(grad_phi,position,tol,max_iterations,0.0005,0.7,0.99,1e-8);

        path = [path,inner_path(:,2:end)];
        c = position - limits;
        lambda = max(0,lambda + mu*c);

        if all(c <= tol) && norm(position-previous) <= tol
            break
        end

        mu = 4*mu;
    end
end