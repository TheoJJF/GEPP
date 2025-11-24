function x = SolveSGEPP(A,b)
    if isstruct(A) == false
        x = SolveGEPP(A,b);
        return
    end

    [U,d] = SGEPP(A,b);
    x = SBackSub(U,d);
end