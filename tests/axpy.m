classdef axpy < matlab.unittest.TestCase
%AXPY Unit tests for the xAXPY class of functions

    properties
        sopts
        hopts
        dopts

        rf

        xint
        yint
        xdec
        ydec
        alpha

        tol
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.sopts.format = 's';
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            if strcmpi( getenv('CHOPBLAS_ROUND_FUNC'), 'cpfloat' )
                testCase.rf = @cpfloat;
            else
                % Default to chop
                testCase.rf = @chop;
            end

            testCase.xint = [1; 2; 3; 4; 5; -3; -2];
            testCase.yint = [2; 3; 4; 5; 6;  7;  8];
            testCase.xdec = [1.0005; 2.34; 3.24; 4];
            testCase.ydec = [2.00125; 3.00225; 4.00014; 5.0025];
            testCase.alpha = 3;

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        % No options specified, use global default of the double mode
        function chop_no_opts(testCase)
            testCase.rf( [], testCase.dopts );

            % Use the provided rounding function
            z = chaxpy( testCase.alpha, testCase.xint, testCase.yint, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, (testCase.alpha.*testCase.xint) + testCase.yint );
        end

        % Specify the rounding function
        function chop_round_func(testCase)
            % Use the default rounding function (chop)
            chop( [], testCase.dopts );

            z = chaxpy( testCase.alpha, testCase.xint, testCase.yint );
            testCase.verifyEqual( z, (testCase.alpha.*testCase.xint) + testCase.yint );

            % Test a trivial rounding function to ensure it overrides it
            z = chaxpy( testCase.alpha, testCase.xint, testCase.yint, 'Rounding', @(x, y) zeros(length(x), 1) );
            testCase.verifyEqual( z, zeros(length(testCase.xint), 1) );
        end

        % Only one mode specified, it is used in both operations
        function chop_one_mode(testCase)
            z = chaxpy( testCase.alpha, testCase.xdec, testCase.ydec, testCase.hopts, 'Rounding', testCase.rf );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(testCase.xdec)
                c(i) = testCase.rf( testCase.rf( testCase.alpha * testCase.xdec(i), testCase.hopts ) + testCase.ydec(i), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end

        % Two modes specified, each operation has a different mode
        function chop_two_modes(testCase)
            z = chaxpy( testCase.alpha, testCase.xdec, testCase.ydec, testCase.sopts, testCase.hopts, 'Rounding', testCase.rf );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(testCase.xdec)
                c(i) = testCase.rf( testCase.rf( testCase.alpha * testCase.xdec(i), testCase.sopts ) + testCase.ydec(i), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end
    end
end
