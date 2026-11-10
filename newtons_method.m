function [position, path] = newtons_method(g, initial_guess, tol, max_iterations, f, H)
% Newton's method with backtracking line search (Task 1.5)
%   g - gradient handle (2x1), f - objective handle, H - Hessian handle (2x2)

    rho = 0.5;    % backtracking factor      
    mu1 = 1e-4;   % sufficient decrease factor 

    x = initial_guess(:);   % column vector
    path = x;

    for i = 1:max_iterations
        xOld = x;
        Gx = g(xOld);       % gradient 
        Hx = H(xOld);       % Hessian 

        % Newton step: solve H*dx = -g
        dx = gauss(Hx, -Gx);

        % Backtracking line search on f, starting from the full Newton step
        l = 1;
        while f(xOld + l*dx) > f(xOld) + mu1*l*(Gx'*dx)
            l = rho*l;
            if l < eps      % no acceptable step exists
                break;
            end
        end

        x = xOld + l*dx;    % update
        path(:, end+1) = x;

        % Stop when the step length is below tol
        if norm(x - xOld) < tol
            break;
        end
    end

    position = x;
end