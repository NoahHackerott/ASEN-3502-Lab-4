function [position, path] = gradient_descent(g, initial_guess, tol, max_iterations, f, l_init)
    rho = 0.5;
    mu1 = 1e-4;
    position = initial_guess(:);
    path = position;

    for k = 1:max_iterations
        % Gradient from the analytic gradient handle
        grad = g(position);
        f0 = f(position);

        % Steepest descent direction
        p = -grad;

        % Backtracking line search (Armijo sufficient decrease on f)
        l = l_init;
        x_trial = position + l*p;
        while f(x_trial) > f0 + mu1*l*(grad'*p)
            l = rho*l;
            x_trial = position + l*p;
            if l < eps        % no acceptable step exists
                break;
            end
        end

        % Step length before accepting the step
        step = norm(x_trial - position);

        % Accept the step and save the new position
        position = x_trial;
        path(:,end+1) = position;

        % Stop when the step length is below tol
        if step < tol
            break;
        end
    end
end