function x = SolveGEPP(A,b,TOL)
    arguments
        A; 
        b;
        TOL = 1e-12;
    end

    [L,U,P] = GEPP(A,TOL);
    if P == eye(size(A,1))
        y = ForSub(L,b);
    else
        y = ForSub(L,P*b);
    end
    x = BackSub(U,y);
end