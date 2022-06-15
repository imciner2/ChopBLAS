classdef axpy < matlab.unittest.TestCase
%AXPY Unit tests for the xAXPY class of functions

    properties
        sopts
        hopts
        dopts
        tol
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.sopts.format = 's';
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            % Set the default format for chop to double precision (mimic no rounding)
            chop( [], testCase.dopts );

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        %% Test the chop-based rounding methods

        % No options specified, use global default of the double mode
        function chop_no_opts(testCase)
            x = [1; 2; 3; 4; 5; -3; -2];
            y = [2; 3; 4; 5; 6;  7;  8];
            alpha = 3;

            chop( [], testCase.dopts );

            z = chaxpy( alpha, x, y );

            testCase.verifyEqual( z, (alpha.*x) + y );
        end

        % Only one mode specified, it is used in both operations
        function chop_one_mode(testCase)
            x = [1.0005; 2.34; 3.24; 4];
            y = [2.00125; 3.00225; 4.00014; 5.0025];
            alpha = 3;

            z = chaxpy( alpha, x, y, testCase.hopts );

            c = zeros(length(x), 1);
            for i=1:length(x)
                c(i) = chop( chop( alpha * x(i), testCase.hopts ) + y(i), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end

        % Two modes specified, each operation has a different mode
        function chop_two_modes(testCase)
            x = [1.0005; 2.34; 3.24; 4];
            y = [2.00125; 3.00225; 4.00014; 5.0025];
            alpha = 3;

            z = chaxpy( alpha, x, y, testCase.sopts, testCase.hopts );

            c = zeros(length(x), 1);
            for i=1:length(x)
                c(i) = chop( chop( alpha * x(i), testCase.sopts ) + y(i), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end
    end
end
