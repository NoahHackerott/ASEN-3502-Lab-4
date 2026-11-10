function [position, path] = gradient_descent(g, initial_guess, tol, max_iterations, l_init)

    rho = 0.5;
    mu1 = 1e-4;

    position = initial_guess(:);
    path = position;

    for k = 1:max_iterations

        % Approximate the gradient using finite differences
        n = length(position);
        grad = zeros(n,1);
        eps_fd = 1e-6;

        g0 = g(position);

        for j = 1:n
            ej = zeros(n,1);
            ej(j) = eps_fd;
            grad(j) = (g(position + ej) - g(position - ej)) ...
                      / (2*eps_fd);
        end

        % Steepest descent direction
        p = -grad;

        % Stop if the gradient is essentially zero
        if norm(grad) < 1e-10
            break;
        end

        % Backtracking line search
        l = l_init;

        while l >= tol
            x_trial = position + l*p;

            if g(x_trial) <= g0 + mu1*l*(grad'*p)
                break;
            end

            l = rho*l;
        end

        % Stop if no acceptable step was found
        if l < tol
            break;
        end

        % Accept the step and save the new position
        position = x_trial;
        path(:,end+1) = position;

    end

    fprintf('Final alpha = %.10f rad\n', position(1));
    fprintf('Final delta_e = %.10f rad\n', position(2));
    fprintf('Number of iterations = %d\n', size(path,2)-1);
end