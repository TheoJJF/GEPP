function [U,d] = SGEPP(A,b)
%{

Perform Sparse Gaussian Elimination with Partial Pivoting (SGEPP).

Parameters: A (Struct/Matrix), b (Vector)
Note: A is assume to be nonsingular
    Error is thrown when the matrix is near/singular.

If A is of the form

A = 

  struct with fields:

    sparse_matrix: {4×1 cell}
                n: 4
                m: 4

A.sparse_matrix = 

  4×1 cell array

    {[dictionary (double --> double) with 2 entries]}
    {[dictionary (double --> double) with 3 entries]}
    {[dictionary (double --> double) with 3 entries]}
    {[dictionary (double --> double) with 2 entries]}

Else, perform standard GEPP

A =

     2    -1     0     0
    -1     2    -1     0
     0    -1     2    -1
     0     0    -1     2


Returns: U (Struct/Matrix), d (Vector)

%}

    if isstruct(A) == false
        [U,d] = GEPP(A,b);
        return
    end

    if A.n ~= size(b,1)
        error("SGEPP:IncompatibleSystem","Incompatible system.");
    end

    if A.n ~= A.m
        warning("SGEPP:NonsquareMatrix","Matrix is not square.");
    end
    
    n = size(A.sparse_matrix,1);

    if A.n ~= n
        error("SGEPP:SingularMatrix","Matrix is singular.");
    end

    U = A.sparse_matrix;
    subcol = zeros(n,1);
    pivot_factor = zeros(n,1);

    for i = 1:n-1
        for j = i:n
            subcol(j) = lookup(U{j},i,FallbackValue=0);
        end

        [Aii_max,i_max] = max(abs(subcol));

        if abs(Aii_max) < 1e-8
            error("SGEPP:SingularMatrix","Matrix is singular.");
        end

        if i ~= i_max
            UTemp = U{i};
            U{i} = U{i_max};
            U{i_max} = UTemp;

            b([i_max,i]) = b([i,i_max]);
        end

        pivot_recip = 1/Aii_max;

        for j = i+1:n
            if lookup(U{j},i,FallbackValue=0) ~= 0
                pivot_factor(j) = U{j}(i)*pivot_recip;
            end
            U{j} = remove(U{j},i);

            b(j) = b(j)-pivot_factor(j)*b(i);
        end

        subcol = subcol*0;

        for j = i+1:n
            for k = i+1:n
                if lookup(U{j},k,FallbackValue=0) ~= 0
                    if  lookup(U{i},k,FallbackValue=0) ~= 0
                        updated_entry = U{j}(k)-pivot_factor(j)*U{i}(k);
                    else
                        updated_entry = U{j}(k);
                    end
                else
                    continue
                end

                if abs(updated_entry) < 1e-8
                    U{j} = remove(U{j},k);
                else
                    U{j}(k) = updated_entry;
                end
            end
        end
    end

    d = b;
end