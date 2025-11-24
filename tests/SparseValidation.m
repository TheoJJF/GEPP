classdef SparseValidation < matlab.unittest.TestCase

    methods (Test)
        
        function case1(testUpperTriCase)
            data = load(fullfile(pwd,"../data","SparseValidation.mat"));
            verification = load(fullfile(pwd,"../data","expected.mat"));

            A = data.A;
            b = data.b;
            expected = verification.expected;

            actual = SGEPP(A,b);
            testUpperTriCase.verifyEqual(actual,expected);
        end

        function case2(testSolutionCase)
            data = load(fullfile(pwd,"../data","SparseValidation.mat"));
            
            sparse_A = data.A;
            sparse_b = data.b;

            A = [2,-1,0,0;-1,2,-1,0;0,-1,2,-1;0,0,-1,2];
            b = ones(4,1);

            actual = SolveSGEPP(sparse_A,sparse_b);
            expected = A\b;

            testSolutionCase.verifyEqual(actual,expected);
        end

        function case3(testSolutionCase)
            data = load(fullfile(pwd,"../data","SparseValidation.mat"));
            
            sparse_A = data.A;
            sparse_b = data.b;

            A = [2,-1,0,0;-1,2,-1,0;0,-1,2,-1;0,0,-1,2];
            b = ones(4,1);

            actual = SolveSGEPP(sparse_A,sparse_b);
            expected = SolveGEPP(A,b);

            testSolutionCase.verifyEqual(actual,expected);
        end
    end

end