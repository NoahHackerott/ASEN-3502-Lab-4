function x = gauss(A,b) %Gauss elim

n = size(A,1);
x = zeros(n,1);



    for i = 1:n-1
        [~, idx] = max(abs(A(i:n,i)));
        idx = idx + i -1;
        A([idx,i], :) = A([i, idx], :);
        b([idx,i]) = b([i, idx]);
        for j = i+1:n
        f = A(j, i)/ A(i, i);
        A(j, i:n) =  A(j, i:n) - f* A(i, i:n);
        b(j) = b(j) - f * b(i);
       
        end
    end

    for i = n:-1:1
        x(i) = (b(i) - A(i, i+1:n) * x(i+1:n))/ A(i,i);
    end

end

