function x = BackSub(U,d)
    n = size(U,1);
    x = zeros(n,1);

    x(n) = d(n)/U(n,n);

    for i = n-1:-1:1
        pivot_recip = 1/U(i,i);

        for j = i+1:n
            x(i) = x(i)+U(i,j)*x(j);
        end
        
        x(i) = (d(i)-x(i))*pivot_recip;
    end
end