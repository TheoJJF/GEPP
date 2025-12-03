function x = SolveGEPP(A,b,TOL)
    [U,d] = GEPP(A,b,TOL);
    x = BackSub(U,d);
end