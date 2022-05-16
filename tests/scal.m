classdef scal < matlab.unittest.TestCase
%SCAL Unit tests for the xSCAL class of functions

    properties
        options
        tol
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.options.format = 's';

            % Set the default format for chop to half precision
            chop( [], testCase.options );

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        % Test the chop-based rounding methods
        function chop_no_opts(testCase)
            alpha = 2;
            x = [1, 2, 3, 4, 5, -3, -2];

            y = chscal( alpha, x );

            testCase.verifyEqual( y, 2.*x );
        end

        function chop_opts(testCase)
            alpha = 20000;
            x = [1, 2, -2, 20, 30];

            opts.format = 'h';

            y = chscal( alpha, x, opts );

            testCase.verifyEqual( y, [20000, 40000, -40000, Inf, Inf] );
        end
    end
end
