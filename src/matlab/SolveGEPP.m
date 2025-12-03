function x = SolveGEPP(A,b,TOL)
    arguments
        A; 
        b;
        TOL = 1e-12;
    end

    [~,U,~,d] = GEPP(A,b,TOL);
    x = BackSub(U,d);
end