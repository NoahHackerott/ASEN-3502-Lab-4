function [position, path] = Adam(g,x,tol,maxit,l,gamv,gams,epsilon)

    x = x(:);

    m = zeros(2,1);
    v = zeros(2,1);

    k = 0;

    path = x;

    while k < maxit

        k = k + 1;

        grad = g(x);

        m = gams*m + (1-gams)*grad;

        v = gamv*v + (1-gamv)*grad.^2;

        mhat = m/(1-gams^k);

        vhat = v/(1-gamv^k);

        xnew = x - l*mhat./(sqrt(vhat)+epsilon);

        path(:,end+1) = xnew;

        if norm(xnew - x) < tol
           x = xnew;
            fprintf('%d iterations occurred\n',k)
        break;
        end
        
        x = xnew;

    end

    position = x;

end