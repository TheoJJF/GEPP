function [L,U,P,d] = GEPP(A,b,TOL)
%{

Perform Gaussian Elimination with Partial Pivoting (GEPP).

Parameters: A (Matrix), b (Vector)
Note: A is assume to be nonsingular
    Error is thrown when the matrix is near/singular.

Returns: L (Matrix), U (Matrix), P (Matrix) d (Vector)

%}

    arguments
        A;
        b = zeros(size(A,1),1);
        TOL = 1e-12;
    end

    [n,m] = size(A);

    if n ~= m
        warning("GEPP:NonsquareMatrix","Matrix is not square.");
    end

    P = eye(size(A));
    L = zeros(size(A));
    pivot_factor = zeros(n,1);
    
    for i = 1:n-1
        [Aii_max,i_max_subcol] = max(abs(A(i:n,i)));
        
        if abs(Aii_max) < TOL
            error("GEPP:SingularMatrix","Matrix is singular.");
        end

        i_max = i-1+i_max_subcol;
        
        if i ~= i_max
            A([i_max,i],:) = A([i,i_max],:);
            L([i_max,i],:) = L([i,i_max],:);
            P([i_max,i],:) = P([i,i_max],:);
            b([i_max,i]) = b([i,i_max]);
        end

        pivot_recip = 1/A(i,i);

        for j = i+1:n
            pivot_factor(j) = A(j,i)*pivot_recip;
            L(j,i) = pivot_factor(j);
            b(j) = b(j)-pivot_factor(j)*b(i);
            A(j,i) = 0;
        end

        for j = i+1:n
            for k = i+1:n
                A(j,k) = A(j,k)-pivot_factor(j)*A(i,k);
            end
        end
    end
    
    if abs(A(n,n)) < TOL
        error("GEPP:SingularMatrix","Matrix is singular.");
    end

    U = A; L = L + eye(n); d = b;
end