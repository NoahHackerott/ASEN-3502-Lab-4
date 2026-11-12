function [position, path, lambda] = augmented_lagrange(g, initial_guess, tol, max_iterations, f, limits)
    position = initial_guess;
    path = position;
    lambda = [0;0];
    mu = 10;
    outer_tol = 1e-5;

    for k = 1:20
        previous = position;
        phi = @(x) f(x) + sum(max(0,lambda + mu*(x-limits)).^2 - lambda.^2)/(2*mu);
        grad_phi = @(x) g(x) + max(0,lambda + mu*(x-limits));
        [position,inner_path] = Adam(grad_phi,position,tol,max_iterations,0.0005,0.99,0.7,1e-8);
        path = [path,inner_path(:,2:end)];
        c = position - limits;
        lambda = max(0,lambda + mu*c);
        if all(c <= tol) && norm(position-previous) <= outer_tol
            break
        end
        mu = 4*mu;
    end
end