classdef dot < matlab.unittest.TestCase
%DOT Unit tests for the xDOT class of functions

    properties
        sopts
        hopts
        dopts

        rf

        xint
        yint
        xdec
        ydec

        tol
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.sopts.format = 's';
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            % Choose the rounding function to use
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

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        % No options specified, use global default of the double mode
        function chop_no_opts(testCase)
            testCase.rf( [], testCase.dopts );

            res = testCase.xint'*testCase.yint;

            % Test with global rounding options
            z = chdot( testCase.xint, testCase.yint, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, res );

            z = chdot( testCase.xint', testCase.yint, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, res );

            z = chdot( testCase.xint, testCase.yint', 'Rounding', testCase.rf );
            testCase.verifyEqual( z, res );

            % Test argument verification
            testCase.verifyError( @() chdot( [2, 3], [3, 4, 5] ), "chdot:xyMustBeCompatibleSize" );
            testCase.verifyError( @() chdot( [2, 3, 4], [3, 4] ), "chdot:xyMustBeCompatibleSize" );
        end

        % Use custom rounding function
        function chop_round_func(testCase)
            % Test with default of chop
            chop( [], testCase.dopts );

            z = chdot( testCase.xint, testCase.yint );
            testCase.verifyEqual( z, testCase.xint'*testCase.yint );

            % Test with trivial rounding function
            z = chdot( testCase.xint, testCase.yint, 'Rounding', @(x, y) zeros(length(x), 1) );
            testCase.verifyEqual( z, 0 );
        end

        % Only one mode specified, it is used in both operations
        function chop_one_mode(testCase)
            z = chdot( testCase.xdec, testCase.ydec, testCase.hopts, 'Rounding', testCase.rf );

            c = 0;
            for i=1:length(testCase.xdec)
                c = testCase.rf( c + chop( testCase.xdec(i) * testCase.ydec(i), testCase.hopts ), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end

        % Two modes specified, each operation has a different mode
        function chop_two_modes(testCase)
            z = chdot( testCase.xdec, testCase.ydec, testCase.sopts, testCase.hopts, 'Rounding', testCase.rf );

            c = 0;
            for i=1:length(testCase.xdec)
                c = testCase.rf( c + testCase.rf( testCase.xdec(i) * testCase.ydec(i), testCase.sopts ), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end
    end
end
