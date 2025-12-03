function x = SolveGEPP(A,b,TOL)
    [L,U,P] = GEPP(A,TOL);
    if P == eye(size(A,1))
        y = ForSub(L,b);
    else
        y = ForSub(L,P*b);
    end
    x = BackSub(U,y);
end