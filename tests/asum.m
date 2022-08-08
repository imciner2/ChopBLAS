classdef asum < matlab.unittest.TestCase
%ASUM Unit tests for the xASUM class of functions

    properties
        sopts
        hopts
        dopts

        rf

        xint
        xdec

        xodd;

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
            testCase.xdec = [1.0005; 2.34; 3.24; 4];

            testCase.xodd = [1.0005; 2.34; 3.24; 4; 5.25678];

            % Set a tolerance for all the tests
            testCase.tol = 1e-7;
        end
    end

    methods(Test)
        % No options specified, use global default of the double mode
        function chop_no_opts(testCase)
            testCase.rf( [], testCase.dopts );;

            % Test with global rounding options
            z = chasum( testCase.xint, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, sum( abs( testCase.xint ) ) );
        end

        % Use custom rounding function
        function chop_round_func(testCase)
            % Test with default of chop
            chop( [], testCase.dopts );

            z = chasum( testCase.xint );
            testCase.verifyEqual( z, sum( abs( testCase.xint ) ) );

            % Test with trivial rounding function
            z = chasum( testCase.xint, 'Rounding', @(x, y) zeros(length(x), 1) );
            testCase.verifyEqual( z, 0 );
        end

        % Change the addition algorithm
        function chop_add_algorithm(testCase)
            % Test with default of chop
            chop( [], testCase.dopts );

            z = chasum( testCase.xodd, 'Rounding', testCase.rf, 'Algorithm', 'recursive' );
            testCase.verifyEqual( z, sum( abs( testCase.xodd ) ), 'AbsTol', testCase.tol );

            % Test with an unknown algorithm specified
            testCase.verifyError( @() chasum( [2, 3, 5], 'Algorithm', 'random' ), "chasum:unknownAlgorithm" );

            % Test with recursive algorithm
            x = abs( testCase.xodd );
            res = x(1);
            for i=2:1:length(x)
                res = double( half( res + x(i) ) );
            end

            z = chasum( testCase.xodd, testCase.hopts, 'Rounding', testCase.rf, 'Algorithm', 'recursive' );
            testCase.verifyEqual( z, res );

            % Test with pairwise algorithm
            x = abs( testCase.xodd );
            i1  = double( half( x(1) + x(3) ) );
            i2  = double( half( x(2) + x(4) ) );
            i3  = double( half( i1 + i2 ) );
            res = double( half( i3 + x(5) ) );

            z = chasum( testCase.xodd, testCase.hopts, 'Rounding', testCase.rf, 'Algorithm', 'pairwise' );
            testCase.verifyEqual( z, res );
        end

        % Only one mode specified, it is used in all operations
        function chop_one_mode(testCase)
            z = chasum( testCase.xodd, testCase.hopts, 'Rounding', testCase.rf );

            x = abs( testCase.xodd );
            c = x(1);
            for i=2:1:length(x)
                c = double( half( c + x(i) ) );
            end

            testCase.verifyEqual( z, c );
        end
    end
end
