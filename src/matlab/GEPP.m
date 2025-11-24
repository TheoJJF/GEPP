function [U,d] = GEPP(A,b)
%{

Perform Gaussian Elimination with Partial Pivoting (GEPP).

Parameters: A (Matrix), b (Vector)
Note: A is assume to be nonsingular
    Error is thrown when the matrix is near/singular.

Returns: U (Matrix), d (Vector)

%}

    [n,m] = size(A);

    if n ~= size(b,1)
        error("GEPP:IncompatibleSystem","Incompatible system.");
    end

    if n ~= m
        warning("GEPP:NonsquareMatrix","Matrix is not square.");
    end

    pivot_factor = zeros(n,1);
    
    for i = 1:n-1
        [Aii_max,i_max_subcol] = max(abs(A(i:n,i)));
        
        if abs(Aii_max) < 1e-6
            error("GEPP:SingularMatrix","Matrix is singular.");
        end

        i_max = i + i_max_subcol - 1;
        
        if i ~= i_max
            A([i_max,i],:) = A([i,i_max],:);
            b([i_max,i]) = b([i,i_max]);
        end

        pivot_recip = 1/A(i,i);

        for j = i+1:n
            pivot_factor(j) = A(j,i)*pivot_recip;
            A(j,i) = 0;
            b(j) = b(j)-pivot_factor(j)*b(i);
        end

        for j = i+1:n
            for k = i+1:n
                A(j,k) = A(j,k)-pivot_factor(j)*A(i,k);
            end
        end
    end
    
    if abs(A(n,n)) < 1e-6
        error("GEPP:SingularMatrix","Matrix is singular.");
    end

    U = A; d = b;
end