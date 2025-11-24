function x = SolveGEPP(A,b)
    [U,d] = GEPP(A,b);
    x = BackSub(U,d);
end