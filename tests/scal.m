classdef scal < matlab.unittest.TestCase
%SCAL Unit tests for the xSCAL class of functions

    properties
        hopts
        dopts

        rf

        tol
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            % Choose the rounding function to use
            if strcmpi( getenv('CHOPBLAS_ROUND_FUNC'), 'cpfloat' )
                testCase.rf = @cpfloat;
            else
                % Default to chop
                testCase.rf = @chop;
            end

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        % Test the chop-based rounding methods
        function chop_no_opts(testCase)
            alpha = 2;
            x = [1, 2, 3, 4, 5, -3, -2];

            % Set the default format for chop to double precision
            testCase.rf( [], testCase.dopts );

            y = chscal( alpha, x, 'Rounding', testCase.rf );
            testCase.verifyEqual( y, 2.*x );

            % Alpha not a scalar
            testCase.verifyError( @() chscal( [2; 2], x ), "chscal:AlphaMustBeScalar" );
        end

        function chop_round_func(testCase)
            alpha = 2;
            x = [1, 2, 3, 4, 5, -3, -2];

            % Test the default rounder using chop
            chop( [], testCase.dopts );

            y = chscal( alpha, x );
            testCase.verifyEqual( y, 2.*x );

            % Test a trivial rounding function
            y = chscal( alpha, x, 'Rounding', @(x, y) zeros(length(x), 1) );
            testCase.verifyEqual( y, zeros(length(x), 1) );
        end

        function chop_opts(testCase)
            alpha = 20000;
            x = [1, 2, -2, 20, 30];

            y = chscal( alpha, x, testCase.hopts, 'Rounding', testCase.rf );

            testCase.verifyEqual( y, [20000, 40000, -40000, Inf, Inf] );
        end
    end
end
