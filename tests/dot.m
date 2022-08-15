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

        xodd
        xoddunsorted

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

            testCase.xodd = [1.0005; 2.34; 3.24; 4; 5.25678];

            testCase.xoddunsorted = testCase.xodd( randperm( length(testCase.xodd) ) );

            % Set a tolerance for all the tests
            testCase.tol = 1e-7;
        end
    end

    methods(Test)
        % No options specified, use global default of the double mode
        function chop_no_opts(testCase)
            testCase.rf( [], testCase.dopts );

            res = testCase.xint'*testCase.yint;

            % Test with global rounding options
            z = chdot( testCase.xint, testCase.yint, ...
                       'Rounding', testCase.rf );
            testCase.verifyEqual( z, res );

            z = chdot( testCase.xint', testCase.yint, ...
                       'Rounding', testCase.rf );
            testCase.verifyEqual( z, res );

            z = chdot( testCase.xint, testCase.yint', ...
                       'Rounding', testCase.rf );
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
            z = chdot( testCase.xint, testCase.yint, ...
                       'Rounding', @(x, y) zeros(length(x), 1) );
            testCase.verifyEqual( z, 0 );
        end

        % Change the addition algorithm
        function chop_add_algorithm(testCase)
            % Test with default of chop
            chop( [], testCase.dopts );

            z = chdot( testCase.xodd, testCase.xodd, ...
                       'Rounding', testCase.rf, ...
                       'Summation', 'recursive' );
            testCase.verifyEqual( z, testCase.xodd'*testCase.xodd, 'AbsTol', testCase.tol );

            % Test with an unknown algorithm specified
            testCase.verifyError( @() chdot( [2, 3, 5], [3, 4, 5], 'Summation', 'random' ), "chdot:unknownSummationAlgorithm" );

            % Test with recursive algorithm
            x = half( testCase.xodd.*testCase.xodd );
            res = x(1);
            for i=2:1:length(x)
                res = half( res + x(i) );
            end

            z = chdot( testCase.xodd, testCase.xodd, testCase.hopts, ...
                       'Rounding', testCase.rf, ...
                       'Summation', 'recursive' );
            testCase.verifyEqual( z, res );

            % Test with pairwise algorithm
            x = half( testCase.xodd.*testCase.xodd );
            res = half( half( half( x(1) + x(3) ) + half( x(2) + x(4) ) ) + x(5) );

            z = chdot( testCase.xodd, testCase.xodd, testCase.hopts, ...
                       'Rounding', testCase.rf, ...
                       'Summation', 'pairwise' );
            testCase.verifyEqual( z, res );

            % Test sorting with increasing magnitude algorithm
            x = half( testCase.xodd.*testCase.xodd );
            res = x(1);
            for i=2:1:length(x)
                res = half( res + x(i) );
            end

            z = chdot( testCase.xoddunsorted, testCase.xoddunsorted, testCase.hopts, ...
                       'Rounding', testCase.rf, ...
                       'Summation', 'increasing' );
            testCase.verifyEqual( z, res );

            % Test sorting with decreasing magnitude algorithm
            x = half( testCase.xodd.*testCase.xodd );
            res = x(end);
            for i=length(x)-1:-1:1
                res = half( res + x(i) );
            end

            z = chdot( testCase.xoddunsorted, testCase.xoddunsorted, testCase.hopts, ...
                       'Rounding', testCase.rf, ...
                       'Summation', 'decreasing' );
            testCase.verifyEqual( z, res );
        end

        % Only one mode specified, it is used in both operations
        function chop_one_mode(testCase)
            z = chdot( testCase.xdec, testCase.ydec, testCase.hopts, ...
                       'Rounding', testCase.rf );

            c = 0;
            for i=1:length(testCase.xdec)
                c = testCase.rf( c + chop( testCase.xdec(i) * testCase.ydec(i), testCase.hopts ), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end

        % Two modes specified, each operation has a different mode
        function chop_two_modes(testCase)
            z = chdot( testCase.xdec, testCase.ydec, testCase.sopts, testCase.hopts, ...
                       'Rounding', testCase.rf );

            c = 0;
            for i=1:length(testCase.xdec)
                c = testCase.rf( c + testCase.rf( testCase.xdec(i) * testCase.ydec(i), testCase.sopts ), testCase.hopts );
            end

            testCase.verifyEqual( z, c );
        end
    end
end
